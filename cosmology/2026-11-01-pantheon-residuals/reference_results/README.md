# Pantheon+ residual comparison — 2026-11-01

This is an in-sample comparison, not a blind prediction or full Gamma_NS derivation.

Selection: noncalibrators with zHD > 0.01; 1580 light curves, 1466 distinct CID. All repeated observations are retained with the supplied full covariance. Observed magnitude is m_b_corr. H0 and absolute magnitude are absorbed into one profiled additive offset for each model. The luminosity-distance prefactor uses zHEL and the radial integral uses zHD.

NS: E^2=exp(alpha*((1+z)^nu-1)), alpha=1.4835435161693789, nu=0.8263146605293827; no refit of alpha/nu. These coefficients have historical Pantheon links.

GR comparator: E^2=Omega_m*(1+z)^3+(1-Omega_m)*(1+z)^2, Lambda=0, curvature distance S_k retained. Omega_m=6.58796983525e-11; optimized on [0,3].

| Model | chi-square, full covariance | ordinary RMS, mag |
|---|---:|---:|
| NS fixed | 1386.816992579 | 0.153399804 |
| GR no Lambda | 1434.628931388 | 0.157957914 |

Delta chi-square (NS minus GR): -47.811938808. RMS is descriptive and does not replace the full-covariance statistic.

This 1,580-row run is separate from the historic A-path 1,701-row saved table. It does not replace or exactly reproduce that table. No new AIC/BIC claim is made.

Binned points are joint generalized least-squares means at fixed display bins. Their covariance includes removal of the fitted additive offset but conditions on model shape; uncertainty in fitted GR Omega_m is not propagated. The number of display bins is not a physical periodicity claim. Top-panel extreme residuals are clipped only for display.

No period, permanent memory, or physical origin is inferred from a residual plot alone. Full Gamma_NS-to-cosmology and third-party reproduction remain outside this calculation.

Reproduce: place Pantheon+SH0ES.dat and Pantheon+SH0ES_STAT+SYS.cov in data/ next to compare_residuals.py, then run python compare_residuals.py. Dependencies: numpy, pandas, scipy, matplotlib. Input SHA256 hashes and numerical checks are in results/summary.json.

Sources: stored project whole-gravity revision 2026-09-27, section 7; official release https://github.com/PantheonPlusSH0ES/DataRelease ; methodology https://arxiv.org/abs/2202.04077 .

The GR optimum reaches the lower boundary Omega_m=0 (empty open limit); it is not a well-determined positive matter-density measurement. This result must not be generalized to all GR cosmologies or LambdaCDM.
