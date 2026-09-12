"""Reproduce Planck 2018 TT binned residual summary.

The released BestFit column is Planck's baseline reference. No independent
Nam Seung-ho C_ell spectrum is claimed or substituted here.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from data_loader import SOURCES, fetch


def main():
    data_path, checksum = fetch("planck_tt_binned")
    table = np.loadtxt(data_path, comments="#")
    ell, observed, err_lo, err_hi, best_fit = table.T
    sigma = 0.5 * (err_lo + err_hi)
    normalized = (observed - best_fit) / sigma
    result = {
        "project": "Nam Seung-ho Law open-data reproduction",
        "author": "Nam Seung-ho",
        "dataset": SOURCES["planck_tt_binned"],
        "dataset_sha256": checksum,
        "n_bins": int(ell.size),
        "reference": "Planck-released BestFit column",
        "diagnostic_chi2_diagonal": float(normalized @ normalized),
        "mean_normalized_residual": float(np.mean(normalized)),
        "rms_normalized_residual": float(np.sqrt(np.mean(normalized**2))),
        "max_abs_normalized_residual": float(np.max(np.abs(normalized))),
        "status": "OPEN",
        "reason": "diagnostic only; covariance, foreground nuisance likelihood, and independent NS TT/TE/EE spectra are not implemented",
    }
    Path("results").mkdir(exist_ok=True)
    output = Path("results/cmb_residual.json")
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
