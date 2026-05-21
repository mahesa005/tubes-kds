"""
Main entry point -- Membrane Transport Simulator.

Interactively asks the user for simulation parameters, then runs all three
transport mechanism simulations, performs sensitivity analysis, and saves plots.

Usage:
    python main.py

Output:
    - Console: simulation summary + sensitivity rankings
    - /outputs/: 4 PNG plots
"""

from simple_diffusion import simulate as sd_sim
from simple_diffusion import D, L, C1_0 as SD_C1_DEF, C2_0 as SD_C2_DEF

from facilitated_diffusion import simulate as fd_sim
from facilitated_diffusion import JMAX, KM, C1_0 as FD_C1_DEF, C2_0 as FD_C2_DEF

from active_transport import simulate as at_sim
from active_transport import (
    NA_IN_0, NA_OUT, K_IN_0, K_OUT,
    ATP_NORMAL, ATP_LOW,
    VMAX_PUMP, KM_NA, KM_K, KM_ATP, KLEAK_NA, KLEAK_K,
)

from sensitivity_analysis import run_all, print_summary
from visualization import generate_all_plots

DIVIDER = '=' * 62


# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------

def _ask_float(label: str, default: float, unit: str = '') -> float:
    """
    Prompt user for a float value. Press Enter to keep the default.
    Repeats until a valid positive number is entered.
    """
    unit_str = f' {unit}' if unit else ''
    while True:
        raw = input(f'    {label} [{default}{unit_str}]: ').strip()
        if raw == '':
            return default
        try:
            val = float(raw)
            if val <= 0:
                print('    Nilai harus lebih dari 0. Coba lagi.')
                continue
            return val
        except ValueError:
            print('    Input tidak valid. Masukkan angka (contoh: 5.0 atau 1e-9).')


def _ask_yes_no(prompt: str, default: bool = True) -> bool:
    """Ask a yes/no question. Returns bool."""
    hint = '[Y/n]' if default else '[y/N]'
    while True:
        raw = input(f'  {prompt} {hint}: ').strip().lower()
        if raw == '':
            return default
        if raw in ('y', 'ya', 'yes'):
            return True
        if raw in ('n', 'tidak', 'no'):
            return False
        print('  Ketik y atau n.')


def _section(title: str) -> None:
    print(f'\n{DIVIDER}')
    print(f'  {title}')
    print(DIVIDER)


# ---------------------------------------------------------------------------
# Configuration input
# ---------------------------------------------------------------------------

def configure_simple_diffusion() -> dict:
    _section('KONFIGURASI: Difusi Sederhana')
    print('  Parameter default (nilai biologis tipikal):')
    print(f'    D     = {D:.2e} m2/s  (koefisien difusi)')
    print(f'    L     = {L:.2e} m     (tebal membran)')
    print(f'    C1_0  = {SD_C1_DEF*1000:.0f} mM          (konsentrasi awal kompartemen 1)')
    print(f'    C2_0  = {SD_C2_DEF*1000:.0f} mM           (konsentrasi awal kompartemen 2)')

    if not _ask_yes_no('Gunakan parameter default?'):
        print('  Masukkan nilai baru (tekan Enter untuk pakai default):')
        d  = _ask_float('D (koef. difusi, m2/s)', D, 'm2/s')
        l  = _ask_float('L (tebal membran, m)', L, 'm')
        c1 = _ask_float('C1_0 (konsentrasi kompartemen 1, mM)', SD_C1_DEF * 1000, 'mM') / 1000
        c2 = _ask_float('C2_0 (konsentrasi kompartemen 2, mM)', SD_C2_DEF * 1000, 'mM') / 1000
        if c1 <= c2:
            print('  Peringatan: C1_0 harus lebih besar dari C2_0 untuk terjadi difusi net.')
    else:
        d, l, c1, c2 = D, L, SD_C1_DEF, SD_C2_DEF

    return {'D': d, 'L': l, 'C1_0': c1, 'C2_0': c2}


def configure_facilitated_diffusion() -> dict:
    _section('KONFIGURASI: Difusi Terfasilitasi')
    print('  Parameter default:')
    print(f'    Jmax  = {JMAX:.2e} mol/m2/s  (kapasitas transpor maksimum)')
    print(f'    Km    = {KM*1000:.0f} mM              (konsentrasi setengah-maksimum)')
    print(f'    C1_0  = {FD_C1_DEF*1000:.0f} mM              (konsentrasi awal kompartemen 1)')
    print(f'    C2_0  = {FD_C2_DEF*1000:.0f} mM               (konsentrasi awal kompartemen 2)')

    if not _ask_yes_no('Gunakan parameter default?'):
        print('  Masukkan nilai baru (tekan Enter untuk pakai default):')
        jmax = _ask_float('Jmax (mol/m2/s)', JMAX, 'mol/m2/s')
        km   = _ask_float('Km (mM)', KM * 1000, 'mM') / 1000
        c1   = _ask_float('C1_0 (mM)', FD_C1_DEF * 1000, 'mM') / 1000
        c2_raw = input(f'    C2_0 (mM) [{FD_C2_DEF*1000:.0f} mM]: ').strip()
        c2   = float(c2_raw) / 1000 if c2_raw else FD_C2_DEF
    else:
        jmax, km, c1, c2 = JMAX, KM, FD_C1_DEF, FD_C2_DEF

    return {'Jmax': jmax, 'Km': km, 'C1_0': c1, 'C2_0': c2}


def configure_active_transport() -> dict:
    _section('KONFIGURASI: Transpor Aktif (Pompa Na+/K+)')
    print('  Parameter default:')
    print(f'    Na+_in  = {NA_IN_0*1000:.0f} mM   (konsentrasi Na+ intrasel awal)')
    print(f'    K+_in   = {K_IN_0*1000:.0f} mM   (konsentrasi K+ intrasel awal)')
    print(f'    ATP normal   = {ATP_NORMAL*1000:.0f} mM  (kondisi sehat)')
    print(f'    ATP abnormal = {ATP_LOW*1000:.1f} mM  (kondisi iskemia / ATP rendah)')
    print(f'    Vmax    = {VMAX_PUMP:.2e} mol/L/s')

    if not _ask_yes_no('Gunakan parameter default?'):
        print('  Masukkan nilai baru (tekan Enter untuk pakai default):')
        na_in    = _ask_float('Na+_in awal (mM)', NA_IN_0 * 1000, 'mM') / 1000
        k_in     = _ask_float('K+_in awal (mM)',  K_IN_0  * 1000, 'mM') / 1000
        atp_norm = _ask_float('ATP normal (mM)',   ATP_NORMAL * 1000, 'mM') / 1000
        atp_low  = _ask_float('ATP abnormal (mM)', ATP_LOW    * 1000, 'mM') / 1000
        vmax     = _ask_float('Vmax (mol/L/s)',    VMAX_PUMP, 'mol/L/s')
        km_na    = _ask_float('Km_Na (mM)',        KM_NA  * 1000, 'mM') / 1000
        km_atp   = _ask_float('Km_ATP (mM)',       KM_ATP * 1000, 'mM') / 1000
    else:
        na_in    = NA_IN_0
        k_in     = K_IN_0
        atp_norm = ATP_NORMAL
        atp_low  = ATP_LOW
        vmax     = VMAX_PUMP
        km_na    = KM_NA
        km_atp   = KM_ATP

    return {
        'Na_in_0':  na_in,
        'Na_out':   NA_OUT,
        'K_in_0':   k_in,
        'K_out':    K_OUT,
        'ATP_normal':   atp_norm,
        'ATP_abnormal': atp_low,
        'Vmax':     vmax,
        'Km_Na':    km_na,
        'Km_K':     KM_K,
        'Km_ATP':   km_atp,
        'k_leak_Na': KLEAK_NA,
        'k_leak_K':  KLEAK_K,
    }


# ---------------------------------------------------------------------------
# Simulation runners (print results)
# ---------------------------------------------------------------------------

def run_simple_diffusion(p: dict) -> tuple:
    _section('HASIL: Difusi Sederhana')
    print(f'  D = {p["D"]:.2e} m2/s | L = {p["L"]:.2e} m')
    print(f'  C_awal: {p["C1_0"]*1000:.1f} mM vs {p["C2_0"]*1000:.1f} mM')

    t, J, C1, C2 = sd_sim(D=p['D'], L=p['L'], C1_0=p['C1_0'], C2_0=p['C2_0'])
    tau_ms = t[-1] / 5 * 1000
    C_eq   = (C1[-1] + C2[-1]) / 2

    print(f'  Konstanta waktu (tau) : {tau_ms:.3f} ms')
    print(f'  Konsentrasi ekuilibrium: {C_eq*1000:.2f} mM')
    print(f'  Fluks awal            : {J[0]:.4e} mol/m2/s')
    print(f'  Fluks akhir           : {J[-1]:.4e} mol/m2/s')
    return t, J, C1, C2


def run_facilitated_diffusion(p: dict) -> tuple:
    _section('HASIL: Difusi Terfasilitasi')
    print(f'  Jmax = {p["Jmax"]:.2e} mol/m2/s | Km = {p["Km"]*1000:.2f} mM')
    print(f'  C_awal: {p["C1_0"]*1000:.1f} mM vs {p["C2_0"]*1000:.1f} mM')

    t, J, C1, C2 = fd_sim(Jmax=p['Jmax'], Km=p['Km'], C1_0=p['C1_0'], C2_0=p['C2_0'])
    C_eq = (C1[-1] + C2[-1]) / 2

    print(f'  Konsentrasi ekuilibrium: {C_eq*1000:.2f} mM')
    print(f'  Fluks awal            : {J[0]:.4e} mol/m2/s')
    print(f'  Fluks akhir           : {J[-1]:.4e} mol/m2/s')
    print(f'  Fluks saat saturasi (C=Km): Jmax/2 = {p["Jmax"]/2:.2e} mol/m2/s')
    return t, J, C1, C2


def run_active_transport(p: dict) -> tuple:
    _section('HASIL: Transpor Aktif (Pompa Na+/K+)')
    print(f'  Vmax = {p["Vmax"]:.2e} mol/L/s | Km_Na = {p["Km_Na"]*1000:.0f} mM')
    print(f'  Na+: {p["Na_in_0"]*1000:.0f} mM (in) vs {p["Na_out"]*1000:.0f} mM (out)')
    print(f'  K+ : {p["K_in_0"]*1000:.0f} mM (in) vs {p["K_out"]*1000:.0f} mM (out)')

    print(f'\n  -- Kondisi Normal (ATP = {p["ATP_normal"]*1000:.1f} mM) --')
    t_n, Jp_n, Na_n, K_n = at_sim(
        Na_in_0=p['Na_in_0'], Na_out=p['Na_out'],
        K_in_0=p['K_in_0'],  K_out=p['K_out'],
        ATP=p['ATP_normal'],  Vmax=p['Vmax'],
        Km_Na=p['Km_Na'],    Km_K=p['Km_K'],
        Km_ATP=p['Km_ATP'],  k_leak_Na=p['k_leak_Na'],
        k_leak_K=p['k_leak_K'],
    )
    print(f'  Na+_in akhir: {Na_n[-1]*1000:.2f} mM  (target fisiologis: ~15 mM)')
    print(f'  K+_in  akhir: {K_n[-1]*1000:.2f} mM  (target fisiologis: ~140 mM)')
    print(f'  Laju pompa SS: {Jp_n[-1]*1e6:.3f} umol/L/s')

    print(f'\n  -- Kondisi Abnormal (ATP = {p["ATP_abnormal"]*1000:.1f} mM, iskemia) --')
    t_a, Jp_a, Na_a, K_a = at_sim(
        Na_in_0=p['Na_in_0'], Na_out=p['Na_out'],
        K_in_0=p['K_in_0'],  K_out=p['K_out'],
        ATP=p['ATP_abnormal'], Vmax=p['Vmax'],
        Km_Na=p['Km_Na'],    Km_K=p['Km_K'],
        Km_ATP=p['Km_ATP'],  k_leak_Na=p['k_leak_Na'],
        k_leak_K=p['k_leak_K'],
    )
    print(f'  Na+_in akhir: {Na_a[-1]*1000:.2f} mM  (meningkat: pompa terganggu)')
    print(f'  K+_in  akhir: {K_a[-1]*1000:.2f} mM  (menurun: pompa terganggu)')
    print(f'  Laju pompa SS: {Jp_a[-1]*1e6:.3f} umol/L/s')
    if Jp_n[-1] > 0:
        print(f'  Penurunan aktivitas pompa: {(1 - Jp_a[-1]/Jp_n[-1])*100:.1f}%')

    return (t_n, Jp_n, Na_n, K_n), (t_a, Jp_a, Na_a, K_a)


def run_sensitivity(sd_cfg: dict, fd_cfg: dict, at_cfg: dict) -> tuple:
    _section('ANALISIS SENSITIVITAS')
    print('  Setiap parameter divariasikan -50%, -20%, +20%, +50% dari baseline.')
    print('  Metrik: mean |elastisitas| = |(dOutput/Output) / (dParam/Param)|')

    sd_sa_params = {k: sd_cfg[k] for k in ('D', 'L', 'C1_0', 'C2_0')}
    fd_sa_params = {k: fd_cfg[k] for k in ('Jmax', 'Km', 'C1_0')}
    at_sa_params = {
        'Vmax':       at_cfg['Vmax'],
        'Km_Na':      at_cfg['Km_Na'],
        'Km_ATP':     at_cfg['Km_ATP'],
        'k_leak_Na':  at_cfg['k_leak_Na'],
        'ATP':        at_cfg['ATP_normal'],
    }

    sd_res, fd_res, at_res = run_all(
        sd_params=sd_sa_params,
        fd_params=fd_sa_params,
        at_params=at_sa_params,
    )
    for res in [sd_res, fd_res, at_res]:
        print_summary(res)

    return sd_res, fd_res, at_res


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print(f'\n{DIVIDER}')
    print('  SIMULASI MEKANISME TRANSPOR MEMBRAN SEL')
    print('  IF3211 Komputasi Domain-Spesifik -- ITB')
    print(DIVIDER)
    print('\n  Program akan meminta konfigurasi parameter untuk setiap model.')
    print('  Tekan Enter untuk menggunakan nilai default biologis.\n')

    # 1. Collect configuration
    sd_cfg = configure_simple_diffusion()
    fd_cfg = configure_facilitated_diffusion()
    at_cfg = configure_active_transport()

    # 2. Run simulations
    sd_data = run_simple_diffusion(sd_cfg)
    fd_data = run_facilitated_diffusion(fd_cfg)
    at_normal_data, at_abnormal_data = run_active_transport(at_cfg)

    # 3. Sensitivity analysis
    sd_res, fd_res, at_res = run_sensitivity(sd_cfg, fd_cfg, at_cfg)

    # 4. Generate plots
    _section('MEMBUAT GRAFIK')
    paths = generate_all_plots(
        sd_data=sd_data,
        fd_data=fd_data,
        at_normal_data=at_normal_data,
        at_abnormal_data=at_abnormal_data,
        sd_params={'D': sd_cfg['D'], 'L': sd_cfg['L']},
        fd_params={'Jmax': fd_cfg['Jmax'], 'Km': fd_cfg['Km']},
        atp_normal=at_cfg['ATP_normal'],
        atp_abnormal=at_cfg['ATP_abnormal'],
        sd_res=sd_res,
        fd_res=fd_res,
        at_res=at_res,
    )
    print('\n  Grafik tersimpan di folder outputs/:')
    for p in paths:
        print(f'    - {p.split("outputs")[-1].lstrip("/\\")}')

    print(f'\n{DIVIDER}')
    print('  Simulasi selesai.')
    print(DIVIDER)


if __name__ == '__main__':
    main()
