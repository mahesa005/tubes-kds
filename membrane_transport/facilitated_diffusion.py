"""
Facilitated Diffusion Simulator using Michaelis-Menten kinetics.

Biology: Facilitated diffusion uses membrane protein channels/carriers to move
molecules down their concentration gradient. Unlike simple diffusion, transport
rate saturates at high concentrations because transporter binding sites become
fully occupied (similar to enzyme kinetics).

J = Jmax * C / (Km + C)   [Michaelis-Menten]

Jmax : maximum flux when all transporters are occupied [mol/m²/s]
Km   : concentration at half-maximum flux (affinity constant) [mol/L]
C    : substrate concentration [mol/L]

Reference: Lodish et al., Molecular Cell Biology, 8th ed., p. 470
"""

import numpy as np
from scipy.integrate import solve_ivp

# ── Baseline parameters ────────────────────────────────────────────────────
# Source: Lodish et al., Molecular Cell Biology, 8th ed.
JMAX = 1e-5        # Maximum flux [mol/m²/s] — GLUT1 glucose transporter order of magnitude
KM = 5e-3          # Michaelis constant [mol/L] — affinity of transporter for substrate
C1_0 = 50e-3       # Initial concentration, compartment 1 (extracellular) [mol/L]
C2_0 = 0e-3        # Initial concentration, compartment 2 (intracellular) [mol/L]
V = 1e-15           # Volume of each compartment [m³]
A = 1e-12           # Membrane area [m²]

T_END = 200.0      # Simulation duration [s]
N_POINTS = 500


def _mm_flux(C: float, Jmax: float, Km: float) -> float:
    """
    Michaelis-Menten flux.

    Derived from transporter kinetics: at low C, flux is linear in C (like
    simple diffusion). At high C, flux plateaus at Jmax as all transporters
    become saturated.
    """
    return Jmax * C / (Km + C)


def simulate(
    Jmax: float = JMAX,
    Km: float = KM,
    C1_0: float = C1_0,
    C2_0: float = C2_0,
    V: float = V,
    A: float = A,
    t_end: float = T_END,
    n_points: int = N_POINTS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Simulate facilitated diffusion between two compartments.

    Net flux depends on the concentration difference through M-M kinetics:
        J_net = Jmax * C1/(Km + C1) - Jmax * C2/(Km + C2)

    Compartment dynamics:
        dC1/dt = -J_net * A / V
        dC2/dt = +J_net * A / V

    Transport stops when C1 == C2 (thermodynamic equilibrium, no energy input).

    Parameters
    ----------
    Jmax    : maximum transport flux [mol/m²/s]
    Km      : Michaelis constant (half-saturation concentration) [mol/L]
    C1_0    : initial concentration in compartment 1 [mol/L]
    C2_0    : initial concentration in compartment 2 [mol/L]
    V       : volume of each compartment [m³]
    A       : membrane area [m²]
    t_end   : simulation duration [s]
    n_points: number of output time points

    Returns
    -------
    t  : time array [s]
    J  : net flux array [mol/m²/s]
    C1 : concentration in compartment 1 over time [mol/L]
    C2 : concentration in compartment 2 over time [mol/L]
    """
    def odes(t, y):
        C1, C2 = y
        C1 = max(C1, 0.0)
        C2 = max(C2, 0.0)
        J_net = _mm_flux(C1, Jmax, Km) - _mm_flux(C2, Jmax, Km)
        dC1 = -J_net * A / V
        dC2 = +J_net * A / V
        return [dC1, dC2]

    t_eval = np.linspace(0, t_end, n_points)
    sol = solve_ivp(odes, [0, t_end], [C1_0, C2_0],
                    t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-12)

    C1 = sol.y[0]
    C2 = sol.y[1]
    J = np.array([
        _mm_flux(max(C1[i], 0), Jmax, Km) - _mm_flux(max(C2[i], 0), Jmax, Km)
        for i in range(len(sol.t))
    ])

    return sol.t, J, C1, C2


def flux_vs_concentration(
    Jmax: float = JMAX,
    Km: float = KM,
    c_max: float = 100e-3,
    n_points: int = 200,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute unidirectional flux as a function of substrate concentration.

    Shows the saturation curve: linear at low [C], plateau at Jmax for high [C].
    Useful for comparing the saturation behavior vs simple diffusion (which is
    always linear).

    Parameters
    ----------
    Jmax     : maximum flux [mol/m²/s]
    Km       : Michaelis constant [mol/L]
    c_max    : maximum concentration to evaluate [mol/L]
    n_points : number of concentration points

    Returns
    -------
    C : concentration array [mol/L]
    J : flux array [mol/m²/s]
    """
    C = np.linspace(0, c_max, n_points)
    J = _mm_flux(C, Jmax, Km)
    return C, J


if __name__ == "__main__":
    t, J, C1, C2 = simulate()
    print("Facilitated Diffusion")
    print(f"  Jmax          : {JMAX:.2e} mol/m2/s")
    print(f"  Km            : {KM*1000:.1f} mM")
    print(f"  Equilibrium concentration: {(C1[-1]+C2[-1])/2*1000:.1f} mM")
    print(f"  Initial flux  : {J[0]:.4e} mol/m2/s")
    print(f"  Final flux    : {J[-1]:.4e} mol/m2/s")
    print(f"  C1 final: {C1[-1]*1000:.2f} mM  |  C2 final: {C2[-1]*1000:.2f} mM")
