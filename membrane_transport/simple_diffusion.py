"""
Simple Diffusion Simulator using Fick's First Law.

Biology: Simple diffusion is passive movement of molecules across a membrane
down their concentration gradient, requiring no energy or transporters.
Rate is proportional to the concentration gradient (Fick's First Law).

The two-compartment system has an exact analytical solution — an exponential
relaxation to equilibrium with time constant τ = V·L / (2·D·A).

Reference: Alberts et al., Molecular Biology of the Cell, 6th ed., p. 616
"""

import numpy as np

# ── Baseline parameters ────────────────────────────────────────────────────
# Source: Alberts et al., Molecular Biology of the Cell, 6th ed.
D = 1e-9           # Diffusion coefficient [m²/s] — small polar molecule in lipid bilayer
L = 7e-9           # Membrane thickness [m] — typical phospholipid bilayer
C1_0 = 100e-3      # Initial concentration, compartment 1 [mol/L]
C2_0 = 10e-3       # Initial concentration, compartment 2 [mol/L]
V = 1e-15           # Volume of each compartment [m³] — ~1 femtoliter
A = 1e-12           # Membrane area [m²] — ~1 µm² patch

N_POINTS = 500     # Number of output points


def simulate(
    D: float = D,
    L: float = L,
    C1_0: float = C1_0,
    C2_0: float = C2_0,
    V: float = V,
    A: float = A,
    n_points: int = N_POINTS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Simulate simple diffusion between two compartments (analytical solution).

    Fick's First Law: J = -D * (C1 - C2) / L   [mol/m²/s]

    ODE system (converting mol/L to mol/m³ internally):
        dC1/dt = -D*A/(V*L) * (C1 - C2)
        dC2/dt = +D*A/(V*L) * (C1 - C2)

    Exact solution via sum/difference variables:
        C_eq = (C1_0 + C2_0) / 2
        τ    = V * L / (2 * D * A)
        C1(t) = C_eq + (C1_0 - C_eq) * exp(-t / τ)
        C2(t) = C_eq - (C1_0 - C_eq) * exp(-t / τ)   [wait, check sign]

    Parameters
    ----------
    D        : diffusion coefficient [m²/s]
    L        : membrane thickness [m]
    C1_0     : initial concentration in compartment 1 [mol/L]
    C2_0     : initial concentration in compartment 2 [mol/L]
    V        : volume of each compartment [m³]
    A        : membrane area [m²]
    n_points : number of output time points

    Returns
    -------
    t  : time array [s]
    J  : flux array [mol/m²/s]  (positive = flow from compartment 1 → 2)
    C1 : concentration in compartment 1 over time [mol/L]
    C2 : concentration in compartment 2 over time [mol/L]
    """
    # Rate constant [s⁻¹]; convert concentrations to mol/m³ (×1000) then back
    k = D * A / (V * L)          # s⁻¹
    tau = 1.0 / (2.0 * k)        # equilibration time constant [s]

    # Simulate to 5 τ to see full approach to equilibrium
    t_end = 5.0 * tau
    t = np.linspace(0, t_end, n_points)

    C_eq = (C1_0 + C2_0) / 2.0
    amplitude = C1_0 - C_eq      # = (C1_0 - C2_0) / 2

    C1 = C_eq + amplitude * np.exp(-t / tau)
    C2 = C_eq - amplitude * np.exp(-t / tau)

    # Fick's flux [mol/m²/s]  — convert mol/L → mol/m³ for the gradient
    J = D * (C1 - C2) * 1000.0 / L

    return t, J, C1, C2


def steady_state_flux(D: float = D, L: float = L,
                      C1_0: float = C1_0, C2_0: float = C2_0) -> float:
    """Return initial (maximum) flux at t=0 [mol/m²/s]."""
    return D * (C1_0 - C2_0) * 1000.0 / L


if __name__ == "__main__":
    t, J, C1, C2 = simulate()
    tau = t[-1] / 5.0
    print(f"Simple Diffusion")
    print(f"  Time constant (tau)    : {tau*1000:.4f} ms")
    print(f"  Equilibrium concentration: {(C1[-1]+C2[-1])/2*1000:.1f} mM")
    print(f"  Initial flux           : {J[0]:.4e} mol/m²/s")
    print(f"  Final flux             : {J[-1]:.4e} mol/m²/s")
    print(f"  C1 final: {C1[-1]*1000:.2f} mM  |  C2 final: {C2[-1]*1000:.2f} mM")
