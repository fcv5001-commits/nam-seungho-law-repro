"""Diagonal-error Pantheon+ comparison for the proposed no-Lambda NS expansion law.

This is an open-data reproduction check, not a full Pantheon+ likelihood.
The published statistical+systematic covariance is not used in v1.0.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import differential_evolution, minimize

from data_loader import SOURCES, fetch

C_KM_S = 299792.458
H0 = 70.0  # absorbed by the fitted magnitude offset


def distance_modulus(z: np.ndarray, expansion) -> np.ndarray:
    grid = np.unique(np.concatenate(([0.0], np.asarray(z, dtype=float))))
    inv_e = 1.0 / expansion(grid)
    dc = (C_KM_S / H0) * cumulative_trapezoid(inv_e, grid, initial=0.0)
    dl = (1.0 + z) * np.interp(z, grid, dc)
    return 5.0 * np.log10(dl) + 25.0


def profile_offset(mu_obs, mu_model, sigma):
    weight = 1.0 / np.square(sigma)
    return float(np.sum(weight * (mu_obs - mu_model)) / np.sum(weight))


def chi2_profiled(mu_obs, mu_model, sigma):
    offset = profile_offset(mu_obs, mu_model, sigma)
    residual = (mu_obs - mu_model - offset) / sigma
    return float(residual @ residual), offset


def load_sample(path: Path):
    table = np.genfromtxt(path, names=True, dtype=None, encoding="utf-8")
    mask = (table["IS_CALIBRATOR"] == 0) & (table["zHD"] > 0.01)
    return (
        table["zHD"][mask].astype(float),
        table["m_b_corr"][mask].astype(float),
        table["m_b_corr_err_DIAG"][mask].astype(float),
    )


def main():
    data_path, checksum = fetch("pantheon_plus")
    z, mu_obs, sigma = load_sample(data_path)

    def base_objective(x):
        omega_m = float(x[0])
        expansion = lambda zz: np.sqrt(1.0 + omega_m * ((1.0 + zz) ** 3 - 1.0))
        return chi2_profiled(mu_obs, distance_modulus(z, expansion), sigma)[0]

    base = minimize(base_objective, x0=[0.3], bounds=[(0.001, 1.5)])
    omega_m = float(base.x[0])
    base_expansion = lambda zz: np.sqrt(1.0 + omega_m * ((1.0 + zz) ** 3 - 1.0))
    base_chi2, base_offset = chi2_profiled(
        mu_obs, distance_modulus(z, base_expansion), sigma
    )

    def ns_objective(x):
        alpha, nu = map(float, x)
        expansion = lambda zz: np.exp(0.5 * alpha * ((1.0 + zz) ** nu - 1.0))
        return chi2_profiled(mu_obs, distance_modulus(z, expansion), sigma)[0]

    seed = differential_evolution(ns_objective, bounds=[(0.01, 4.0), (0.05, 2.0)], seed=20260912)
    ns = minimize(ns_objective, x0=seed.x, bounds=[(0.01, 4.0), (0.05, 2.0)])
    alpha, nu = map(float, ns.x)
    ns_expansion = lambda zz: np.exp(0.5 * alpha * ((1.0 + zz) ** nu - 1.0))
    ns_chi2, ns_offset = chi2_profiled(mu_obs, distance_modulus(z, ns_expansion), sigma)

    result = {
        "project": "Nam Seung-ho Law open-data reproduction",
        "author": "Nam Seung-ho",
        "dataset": SOURCES["pantheon_plus"],
        "dataset_sha256": checksum,
        "sample_rule": "IS_CALIBRATOR == 0 and zHD > 0.01",
        "n": int(z.size),
        "likelihood_scope": "diagonal errors only; full covariance not included",
        "baseline_no_lambda": {
            "formula": "E^2=1+Omega_m*((1+z)^3-1)",
            "omega_m": omega_m,
            "offset": base_offset,
            "chi2": base_chi2,
        },
        "ns_candidate": {
            "formula": "E^2=exp(alpha*((1+z)^nu-1))",
            "alpha": alpha,
            "nu": nu,
            "offset": ns_offset,
            "chi2": ns_chi2,
        },
        "delta_chi2_ns_minus_baseline": ns_chi2 - base_chi2,
        "status": "OPEN",
        "reason": "exploratory in-sample diagonal fit; not a held-out or full-covariance validation",
    }
    Path("results").mkdir(exist_ok=True)
    output = Path("results/pantheon_fit.json")
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
