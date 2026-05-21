"""
Visualization module for all membrane transport simulations.

Each plot function accepts pre-computed simulation data as arguments,
so parameters are controlled from main.py rather than hardcoded here.

Generates four figures saved to /outputs/:
  1. Concentration vs Time  -- all three mechanisms
  2. Flux vs Concentration  -- simple vs facilitated (saturation effect)
  3. Sensitivity Analysis   -- parameter ranking bar chart
  4. Normal vs Abnormal     -- Na+/K+ pump under healthy vs low-ATP conditions
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from facilitated_diffusion import flux_vs_concentration
from simple_diffusion import D, L

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'outputs')
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLORS = {
    'simple':      '#2196F3',
    'facilitated': '#4CAF50',
    'active':      '#F44336',
    'normal':      '#1565C0',
    'abnormal':    '#B71C1C',
}
plt.rcParams.update({
    'figure.dpi': 120,
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'legend.fontsize': 9,
    'axes.spines.top': False,
    'axes.spines.right': False,
})


def plot_concentration_vs_time(
    sd_data: tuple,
    fd_data: tuple,
    at_data: tuple,
) -> str:
    """
    Plot 1: Concentration vs Time for all three transport mechanisms.

    Parameters
    ----------
    sd_data : (t, J, C1, C2) from simple_diffusion.simulate()
    fd_data : (t, J, C1, C2) from facilitated_diffusion.simulate()
    at_data : (t, J_pump, Na_in, K_in) from active_transport.simulate()

    Returns path to saved PNG.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    fig.suptitle('Membrane Transport: Concentration vs Time', fontweight='bold', y=1.01)

    # Panel 1: Simple Diffusion
    ax = axes[0]
    t, J, C1, C2 = sd_data
    t_ms = t * 1000
    ax.plot(t_ms, C1 * 1000, color=COLORS['simple'], lw=2, label='Compartment 1 (high)')
    ax.plot(t_ms, C2 * 1000, color=COLORS['simple'], lw=2, ls='--', label='Compartment 2 (low)')
    ax.axhline((C1[-1] + C2[-1]) / 2 * 1000, color='gray', lw=0.8, ls=':', label='Equilibrium')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Concentration (mM)')
    ax.set_title("Simple Diffusion\n(Fick's Law)")
    ax.legend()

    # Panel 2: Facilitated Diffusion
    ax = axes[1]
    t, J, C1, C2 = fd_data
    ax.plot(t, C1 * 1000, color=COLORS['facilitated'], lw=2, label='Compartment 1 (high)')
    ax.plot(t, C2 * 1000, color=COLORS['facilitated'], lw=2, ls='--', label='Compartment 2 (low)')
    ax.axhline((C1[-1] + C2[-1]) / 2 * 1000, color='gray', lw=0.8, ls=':', label='Equilibrium')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Concentration (mM)')
    ax.set_title('Facilitated Diffusion\n(Michaelis-Menten)')
    ax.legend()

    # Panel 3: Active Transport
    ax = axes[2]
    t, Jp, Na, K = at_data
    ax.plot(t, Na * 1000, color=COLORS['active'], lw=2, label='[Na+]_in')
    ax.plot(t, K * 1000, color='#9C27B0', lw=2, label='[K+]_in')
    ax.axhline(15, color=COLORS['active'], lw=0.8, ls=':', alpha=0.6, label='Target Na+ (15 mM)')
    ax.axhline(140, color='#9C27B0', lw=0.8, ls=':', alpha=0.6, label='Target K+ (140 mM)')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Concentration (mM)')
    ax.set_title('Active Transport\n(Na+/K+ ATPase Pump)')
    ax.legend(fontsize=8)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'plot1_concentration_vs_time.png')
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


def plot_flux_vs_concentration(
    sd_params: dict,
    fd_params: dict,
) -> str:
    """
    Plot 2: Flux vs Concentration -- simple vs facilitated diffusion.

    Parameters
    ----------
    sd_params : dict with keys 'D', 'L'
    fd_params : dict with keys 'Jmax', 'Km'

    Returns path to saved PNG.
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    Jmax = fd_params['Jmax']
    Km   = fd_params['Km']
    D_   = sd_params['D']
    L_   = sd_params['L']

    C_range, J_fd = flux_vs_concentration(Jmax=Jmax, Km=Km)
    C_mM = C_range * 1000

    J_simple = D_ * C_range * 1000 / L_
    scale = Jmax / J_simple.max() if J_simple.max() > 0 else 1.0
    J_simple_scaled = J_simple * scale

    ax.plot(C_mM, J_fd * 1e5, color=COLORS['facilitated'], lw=2.5,
            label='Facilitated diffusion (Michaelis-Menten)')
    ax.plot(C_mM, J_simple_scaled * 1e5, color=COLORS['simple'], lw=2.5,
            ls='--', label='Simple diffusion (scaled, linear)')
    ax.axhline(Jmax * 1e5, color=COLORS['facilitated'], lw=1, ls=':',
               alpha=0.7, label=f'Jmax = {Jmax:.2e} mol/m2/s')
    ax.axvline(Km * 1000, color='gray', lw=1, ls=':',
               alpha=0.7, label=f'Km = {Km*1000:.1f} mM')

    ax.set_xlabel('Substrate Concentration (mM)')
    ax.set_ylabel('Flux (x10$^{-5}$ mol/m$^2$/s)')
    ax.set_title('Flux vs Concentration: Saturation in Facilitated Diffusion')
    ax.legend()

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'plot2_flux_vs_concentration.png')
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


def plot_sensitivity_analysis(sd_res: dict, fd_res: dict, at_res: dict) -> str:
    """
    Plot 3: Sensitivity analysis -- parameter ranking bar chart.

    Parameters
    ----------
    sd_res, fd_res, at_res : dicts returned by sensitivity_analysis.analyze_*()

    Returns path to saved PNG.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle('Sensitivity Analysis: Parameter Impact Rankings', fontweight='bold', y=1.01)

    model_colors = [COLORS['simple'], COLORS['facilitated'], COLORS['active']]

    for ax, res, color in zip(axes, [sd_res, fd_res, at_res], model_colors):
        params = res['rankings']
        means = [np.mean([abs(s) for _, s in res[p]]) for p in params]

        bars = ax.barh(params[::-1], means[::-1], color=color, alpha=0.8, edgecolor='white')
        for bar, val in zip(bars, means[::-1]):
            ax.text(val + 0.01, bar.get_y() + bar.get_height() / 2,
                    f'{val:.2f}', va='center', fontsize=8)

        metric_label = 'flux' if 'Active' not in res['model'] else '[Na+]_in'
        ax.set_xlabel('Mean |Sensitivity Index|')
        ax.set_title(f"{res['model']}\n(metric: {metric_label})")
        ax.axvline(1.0, color='gray', lw=0.8, ls='--', alpha=0.5)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'plot3_sensitivity_analysis.png')
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


def plot_normal_vs_abnormal(
    at_normal_data: tuple,
    at_abnormal_data: tuple,
    atp_normal: float,
    atp_abnormal: float,
) -> str:
    """
    Plot 4: Normal vs abnormal cell -- Na+/K+ pump with normal vs low ATP.

    Parameters
    ----------
    at_normal_data   : (t, J_pump, Na_in, K_in) with normal ATP
    at_abnormal_data : (t, J_pump, Na_in, K_in) with low ATP
    atp_normal       : normal ATP value [mol/L]
    atp_abnormal     : abnormal ATP value [mol/L]

    Returns path to saved PNG.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Normal vs Abnormal Cell: Na+/K+ Pump Under ATP Depletion',
                 fontweight='bold', y=1.01)

    t_n, _, Na_n, K_n = at_normal_data
    t_a, _, Na_a, K_a = at_abnormal_data

    # Panel 1: Na+
    ax = axes[0]
    ax.plot(t_n, Na_n * 1000, color=COLORS['normal'], lw=2.5,
            label=f'Normal (ATP = {atp_normal*1000:.1f} mM)')
    ax.plot(t_a, Na_a * 1000, color=COLORS['abnormal'], lw=2.5, ls='--',
            label=f'Abnormal (ATP = {atp_abnormal*1000:.1f} mM)')
    ax.axhline(15, color='gray', lw=1, ls=':', alpha=0.6, label='Physiological target (15 mM)')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('[Na+]_intracellular (mM)')
    ax.set_title('Intracellular Na+')
    ax.legend()

    # Panel 2: K+
    ax = axes[1]
    ax.plot(t_n, K_n * 1000, color=COLORS['normal'], lw=2.5,
            label=f'Normal (ATP = {atp_normal*1000:.1f} mM)')
    ax.plot(t_a, K_a * 1000, color=COLORS['abnormal'], lw=2.5, ls='--',
            label=f'Abnormal (ATP = {atp_abnormal*1000:.1f} mM)')
    ax.axhline(140, color='gray', lw=1, ls=':', alpha=0.6, label='Physiological target (140 mM)')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('[K+]_intracellular (mM)')
    ax.set_title('Intracellular K+')
    ax.legend()

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, 'plot4_normal_vs_abnormal.png')
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    return path


def generate_all_plots(
    sd_data, fd_data, at_normal_data, at_abnormal_data,
    sd_params, fd_params,
    atp_normal, atp_abnormal,
    sd_res, fd_res, at_res,
) -> list:
    """Run all four visualization functions and return list of saved file paths."""
    paths = []
    print('  Generating plot 1: Concentration vs Time...')
    paths.append(plot_concentration_vs_time(sd_data, fd_data, at_normal_data))
    print('  Generating plot 2: Flux vs Concentration...')
    paths.append(plot_flux_vs_concentration(sd_params, fd_params))
    print('  Generating plot 3: Sensitivity Analysis...')
    paths.append(plot_sensitivity_analysis(sd_res, fd_res, at_res))
    print('  Generating plot 4: Normal vs Abnormal...')
    paths.append(plot_normal_vs_abnormal(at_normal_data, at_abnormal_data, atp_normal, atp_abnormal))
    return paths
