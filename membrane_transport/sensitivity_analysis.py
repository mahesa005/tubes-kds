"""
Sensitivity Analysis for all three transport mechanisms.

For each model, each parameter is varied by -50%, -20%, +20%, +50% from
its baseline. The impact is measured as the normalized change in the
output metric (steady-state flux or final concentration).

Sensitivity score (relative sensitivity index):
    S = (delta_output / output_base) / (delta_param / param_base)

This is the standard local sensitivity / elasticity coefficient.

Reference: Saltelli, A. et al. (2008) Global Sensitivity Analysis: The Primer.
"""

import numpy as np
from simple_diffusion import simulate as sd_simulate, D, L, C1_0, C2_0, V, A
from facilitated_diffusion import simulate as fd_simulate, JMAX, KM
from facilitated_diffusion import C1_0 as FD_C1_0, C2_0 as FD_C2_0
from active_transport import simulate as at_simulate
from active_transport import (
    NA_IN_0, NA_OUT, K_IN_0, K_OUT, ATP_NORMAL,
    VMAX_PUMP, KM_NA, KM_K, KM_ATP, KLEAK_NA, KLEAK_K
)

# Perturbation fractions to apply
PERTURBATIONS = [-0.50, -0.20, +0.20, +0.50]


def _sensitivity_score(base_val: float, perturb_frac: float,
                        output_base: float, output_perturbed: float) -> float:
    """
    Compute normalized sensitivity (elasticity) coefficient.

    S = (dOutput/Output_base) / (dParam/Param_base)
    """
    if abs(output_base) < 1e-30 or abs(perturb_frac) < 1e-15:
        return 0.0
    delta_out = output_perturbed - output_base
    return (delta_out / output_base) / perturb_frac


def analyze_simple_diffusion(base_params: dict = None) -> dict:
    """
    Sensitivity analysis for simple diffusion.

    Output metric: initial steady-state flux J_0 = D * (C1_0 - C2_0) * 1000 / L
    Parameters varied: D, L, C1_0, C2_0

    Parameters
    ----------
    base_params : dict with keys 'D', 'L', 'C1_0', 'C2_0'. If None, uses module defaults.

    Returns a dict mapping parameter name -> list of (perturbation, score) tuples,
    plus 'rankings' key with parameters sorted by mean |sensitivity|.
    """
    def metric(params):
        t, J, c1, c2 = sd_simulate(
            D=params['D'], L=params['L'], C1_0=params['C1_0'], C2_0=params['C2_0']
        )
        return J[0]   # initial flux

    if base_params is None:
        base_params = {'D': D, 'L': L, 'C1_0': C1_0, 'C2_0': C2_0}
    base_metric = metric(base_params)

    results = {}
    for param_name, base_val in base_params.items():
        scores = []
        for frac in PERTURBATIONS:
            perturbed = {k: v for k, v in base_params.items()}
            perturbed[param_name] = base_val * (1 + frac)
            m = metric(perturbed)
            scores.append((frac, _sensitivity_score(base_val, frac, base_metric, m)))
        results[param_name] = scores

    rankings = sorted(results.keys(),
                      key=lambda p: np.mean([abs(s) for _, s in results[p]]),
                      reverse=True)
    results['rankings'] = rankings
    results['model'] = 'Simple Diffusion'
    return results


def analyze_facilitated_diffusion(base_params: dict = None) -> dict:
    """
    Sensitivity analysis for facilitated diffusion.

    Output metric: initial net flux J_0 at t=0 (= Jmax * C1_0 / (Km + C1_0))
    Parameters varied: Jmax, Km, C1_0

    Parameters
    ----------
    base_params : dict with keys 'Jmax', 'Km', 'C1_0'. If None, uses module defaults.

    Returns same structure as analyze_simple_diffusion().
    """
    def metric(params):
        t, J, c1, c2 = fd_simulate(
            Jmax=params['Jmax'], Km=params['Km'], C1_0=params['C1_0']
        )
        return J[0]

    if base_params is None:
        base_params = {'Jmax': JMAX, 'Km': KM, 'C1_0': FD_C1_0}
    base_metric = metric(base_params)

    results = {}
    for param_name, base_val in base_params.items():
        scores = []
        for frac in PERTURBATIONS:
            perturbed = {k: v for k, v in base_params.items()}
            perturbed[param_name] = base_val * (1 + frac)
            m = metric(perturbed)
            scores.append((frac, _sensitivity_score(base_val, frac, base_metric, m)))
        results[param_name] = scores

    rankings = sorted(results.keys(),
                      key=lambda p: np.mean([abs(s) for _, s in results[p]]),
                      reverse=True)
    results['rankings'] = rankings
    results['model'] = 'Facilitated Diffusion'
    return results


def analyze_active_transport(base_params: dict = None) -> dict:
    """
    Sensitivity analysis for active transport.

    Output metric: steady-state intracellular [Na+] (final value of simulation).
    Lower Na+_in means the pump is working better.
    Parameters varied: Vmax, Km_Na, Km_ATP, k_leak_Na, ATP

    Parameters
    ----------
    base_params : dict with keys 'Vmax', 'Km_Na', 'Km_ATP', 'k_leak_Na', 'ATP'.
                  If None, uses module defaults.

    Returns same structure as analyze_simple_diffusion().
    """
    def metric(params):
        t, Jp, Na, K = at_simulate(
            Vmax=params['Vmax'], Km_Na=params['Km_Na'], Km_ATP=params['Km_ATP'],
            k_leak_Na=params['k_leak_Na'], ATP=params['ATP']
        )
        return Na[-1]   # final [Na+]_in -- lower = healthier pump

    if base_params is None:
        base_params = {
            'Vmax': VMAX_PUMP,
            'Km_Na': KM_NA,
            'Km_ATP': KM_ATP,
            'k_leak_Na': KLEAK_NA,
            'ATP': ATP_NORMAL,
        }
    base_metric = metric(base_params)

    results = {}
    for param_name, base_val in base_params.items():
        scores = []
        for frac in PERTURBATIONS:
            perturbed = {k: v for k, v in base_params.items()}
            perturbed[param_name] = base_val * (1 + frac)
            m = metric(perturbed)
            scores.append((frac, _sensitivity_score(base_val, frac, base_metric, m)))
        results[param_name] = scores

    rankings = sorted(results.keys(),
                      key=lambda p: np.mean([abs(s) for _, s in results[p]]),
                      reverse=True)
    results['rankings'] = rankings
    results['model'] = 'Active Transport (Na+/K+ Pump)'
    return results


def run_all(
    sd_params: dict = None,
    fd_params: dict = None,
    at_params: dict = None,
) -> tuple[dict, dict, dict]:
    """
    Run sensitivity analysis for all three models and return results.

    Parameters
    ----------
    sd_params : base params for simple diffusion (optional, uses defaults if None)
    fd_params : base params for facilitated diffusion (optional)
    at_params : base params for active transport (optional)
    """
    return (
        analyze_simple_diffusion(sd_params),
        analyze_facilitated_diffusion(fd_params),
        analyze_active_transport(at_params),
    )


def print_summary(results: dict) -> None:
    """Print a human-readable sensitivity summary for one model."""
    print(f"\n{'='*55}")
    print(f"  {results['model']}")
    print(f"{'='*55}")
    print(f"  {'Parameter':<15} {'Mean |Sensitivity|':>20}  Rank")
    print(f"  {'-'*45}")
    for rank, param in enumerate(results['rankings'], 1):
        mean_s = np.mean([abs(s) for _, s in results[param]])
        print(f"  {param:<15} {mean_s:>20.4f}  #{rank}")


if __name__ == "__main__":
    sd, fd, at = run_all()
    for r in [sd, fd, at]:
        print_summary(r)
