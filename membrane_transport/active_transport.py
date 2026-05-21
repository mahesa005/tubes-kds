"""
Active Transport Simulator — Na+/K+ ATPase Pump.

Biology: The Na+/K+ pump (Na+/K+-ATPase) uses ATP hydrolysis to pump
3 Na+ OUT and 2 K+ IN per cycle, working AGAINST their electrochemical gradients.
This creates and maintains the resting membrane potential critical for nerve
and muscle function.

ODE system:
    d[Na+]_in/dt = J_leak_Na - 3 * J_pump
    d[K+]_in/dt  = -J_leak_K + 2 * J_pump

J_pump   : pump rate, depends on [Na+]_in, [K+]_out, and [ATP]
J_leak_X : passive leak flux proportional to gradient (always present)

Reference: Skou, J.C. (1957) Biochim. Biophys. Acta 23:394–401 (Nobel Prize 1997)
           Bhatt, D.L. et al. — typical physiological values from Koeppen & Stanton,
           Berne & Levy Physiology, 6th ed., p. 5
"""

import numpy as np
from scipy.integrate import solve_ivp

# ── Baseline parameters ────────────────────────────────────────────────────
# Physiological concentrations [mmol/L = mM]
# Source: Koeppen & Stanton, Berne & Levy Physiology, 6th ed., Table 1-1
NA_IN_0 = 15.0e-3    # [Na+] intracellular, normal  [mol/L]  (15 mM)
NA_OUT = 145.0e-3    # [Na+] extracellular, fixed   [mol/L]  (145 mM)
K_IN_0 = 140.0e-3   # [K+]  intracellular, normal  [mol/L]  (140 mM)
K_OUT = 5.0e-3      # [K+]  extracellular, fixed   [mol/L]  (5 mM)

# ATP concentration [mol/L]
ATP_NORMAL = 5.0e-3   # Normal cell ATP [mol/L]   (5 mM)
ATP_LOW = 0.5e-3      # Low ATP (ischemia/disease) [mol/L]  (0.5 mM)

# Pump kinetics
# Source: Lauger, P. (1991) Electrogenic Ion Pumps, p. 178
# Vmax tuned so that at physiological Na+ (15 mM), K+ (140 mM), ATP (5 mM),
# pump flux exactly balances passive leaks, producing a stable steady state.
VMAX_PUMP = 5.0e-5    # Max pump rate [mol/L/s] at saturating conditions
KM_NA = 10.0e-3       # [Na+]_in affinity [mol/L]  (10 mM)
KM_K = 1.5e-3         # [K+]_out affinity [mol/L]  (1.5 mM)
KM_ATP = 0.5e-3       # ATP affinity [mol/L]        (0.5 mM)

# Passive leak coefficients [s⁻¹]
# Chosen so that at the physiological steady state, leak = pump * stoichiometry:
#   k_leak_Na = 3 * Jp_ss / (Na_out - Na_in_ss)  ~5e-4
#   k_leak_K  = 2 * Jp_ss / (K_in_ss - K_out)   ~3.1e-4
KLEAK_NA = 5.0e-4     # Na+ leak rate constant [s⁻¹]
KLEAK_K = 3.1e-4      # K+ leak rate constant [s⁻¹]

T_END = 2000.0        # Simulation duration [s]  — steady state takes ~30 min
N_POINTS = 1000


def _pump_rate(Na_in: float, K_out: float, ATP: float,
               Vmax: float, Km_Na: float, Km_K: float, Km_ATP: float) -> float:
    """
    Na+/K+ ATPase pump rate using cooperative Michaelis-Menten kinetics.

    The pump requires intracellular Na+ (3 sites, Hill ~2.7), extracellular K+
    (2 sites, Hill ~1.5), and ATP. Simplified here as independent M-M terms.

    J_pump = Vmax * ([Na+]_in / (Km_Na + [Na+]_in))
                  * ([K+]_out / (Km_K  + [K+]_out))
                  * ([ATP]    / (Km_ATP + [ATP]))
    """
    f_Na = Na_in / (Km_Na + Na_in) if Na_in > 0 else 0.0
    f_K = K_out / (Km_K + K_out) if K_out > 0 else 0.0
    f_ATP = ATP / (Km_ATP + ATP) if ATP > 0 else 0.0
    return Vmax * f_Na * f_K * f_ATP


def simulate(
    Na_in_0: float = NA_IN_0,
    Na_out: float = NA_OUT,
    K_in_0: float = K_IN_0,
    K_out: float = K_OUT,
    ATP: float = ATP_NORMAL,
    Vmax: float = VMAX_PUMP,
    Km_Na: float = KM_NA,
    Km_K: float = KM_K,
    Km_ATP: float = KM_ATP,
    k_leak_Na: float = KLEAK_NA,
    k_leak_K: float = KLEAK_K,
    t_end: float = T_END,
    n_points: int = N_POINTS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Simulate Na+/K+ pump dynamics using an ODE system.

    The pump continuously works against passive leak to maintain
    intracellular concentrations far from equilibrium.

    ODEs (intracellular concentrations change; extracellular fixed):
        d[Na+]_in/dt = k_leak_Na * (Na_out - Na_in) - 3 * J_pump
        d[K+]_in/dt  = k_leak_K  * (K_out  - K_in)  + 2 * J_pump

    Parameters
    ----------
    Na_in_0   : initial intracellular [Na+] [mol/L]
    Na_out    : fixed extracellular [Na+] [mol/L]
    K_in_0    : initial intracellular [K+] [mol/L]
    K_out     : fixed extracellular [K+] [mol/L]
    ATP       : ATP concentration [mol/L]  (constant — simplified)
    Vmax      : maximum pump rate [mol/L/s]
    Km_Na     : Na+ intracellular affinity [mol/L]
    Km_K      : K+ extracellular affinity [mol/L]
    Km_ATP    : ATP affinity [mol/L]
    k_leak_Na : Na+ passive leak rate constant [s⁻¹]
    k_leak_K  : K+ passive leak rate constant [s⁻¹]
    t_end     : simulation duration [s]
    n_points  : number of output time points

    Returns
    -------
    t      : time array [s]
    J_pump : pump rate array [mol/L/s]
    Na_in  : intracellular [Na+] over time [mol/L]
    K_in   : intracellular [K+] over time [mol/L]
    """
    def odes(t, y):
        Na_in, K_in = y
        Na_in = max(Na_in, 0.0)
        K_in = max(K_in, 0.0)

        Jp = _pump_rate(Na_in, K_out, ATP, Vmax, Km_Na, Km_K, Km_ATP)

        J_leak_Na = k_leak_Na * (Na_out - Na_in)   # passive Na+ influx
        J_leak_K = k_leak_K * (K_out - K_in)        # passive K+ efflux (K_in > K_out)

        dNa_in = J_leak_Na - 3.0 * Jp
        dK_in = J_leak_K + 2.0 * Jp

        return [dNa_in, dK_in]

    t_eval = np.linspace(0, t_end, n_points)
    sol = solve_ivp(odes, [0, t_end], [Na_in_0, K_in_0],
                    t_eval=t_eval, method='RK45', rtol=1e-9, atol=1e-14)

    Na_in = sol.y[0]
    K_in = sol.y[1]
    J_pump = np.array([
        _pump_rate(max(Na_in[i], 0), K_out, ATP, Vmax, Km_Na, Km_K, Km_ATP)
        for i in range(len(sol.t))
    ])

    return sol.t, J_pump, Na_in, K_in


if __name__ == "__main__":
    print("Active Transport — Na+/K+ ATPase Pump")

    t, Jp, Na, K = simulate(ATP=ATP_NORMAL)
    print(f"\n  [Normal ATP = {ATP_NORMAL*1000:.1f} mM]")
    print(f"  Na+ initial: {NA_IN_0*1000:.0f} mM  ->  final: {Na[-1]*1000:.1f} mM  (target: ~15 mM)")
    print(f"  K+  initial: {K_IN_0*1000:.0f} mM  ->  final: {K[-1]*1000:.1f} mM  (target: ~140 mM)")
    print(f"  Pump rate at SS: {Jp[-1]*1e6:.4f} umol/L/s")

    t2, Jp2, Na2, K2 = simulate(ATP=ATP_LOW)
    print(f"\n  [Low ATP = {ATP_LOW*1000:.1f} mM]")
    print(f"  Na+ initial: {NA_IN_0*1000:.0f} mM  ->  final: {Na2[-1]*1000:.1f} mM")
    print(f"  K+  initial: {K_IN_0*1000:.0f} mM  ->  final: {K2[-1]*1000:.1f} mM")
    print(f"  Pump rate at SS: {Jp2[-1]*1e6:.4f} umol/L/s")
