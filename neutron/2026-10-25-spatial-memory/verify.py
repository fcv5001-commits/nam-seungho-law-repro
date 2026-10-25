"""Verify the published NS spatial-memory response and calibrated readback.

This checks reproduction and linear stability, not a full microscopic neutron
derivation or an independent experimental prediction of the supplied inputs.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import subprocess
import sys
import tempfile

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import eigh, expm
import run as model

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='verification.json')
    args = parser.parse_args()
    checks = {}
    data = model.load()
    with tempfile.TemporaryDirectory(prefix='ns-neutron-') as tmp:
        outputs = []
        for lapse in (1.0, 1.000001):
            out = Path(tmp) / ('baseline.json' if lapse == 1 else 'lapse.json')
            subprocess.run([sys.executable, str(ROOT / 'run.py'), '--lapse',
                            str(lapse), '--output', str(out)], check=True)
            outputs.append(json.loads(out.read_text()))
        baseline, lapse = outputs
    archive = json.loads((ROOT / 'reference_results.json').read_text())
    delta = max(float(np.max(abs(np.asarray(x['poles'])-np.asarray(y['poles']))))
                for x, y in zip(baseline['spatial_samples'], archive['spatial_samples']))
    assert delta < 1e-12
    checks['archived_pole_max_difference'] = delta
    checks['source_sha256_matches_runtime'] = all(
        hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
        for name, digest in baseline['sha256'].items())
    assert checks['source_sha256_matches_runtime']
    pub = data['public_inputs']
    cal = baseline['public_input_effective_neutron']
    lcal = lapse['public_input_effective_neutron']
    tau = -pub['hbar_MeV_s'] / (2 * cal['rest_pole_relative_energy_MeV'][1])
    assert abs(tau-pub['lifetime_s']) < 1e-10
    assert abs(lcal['coordinate_lifetime_s']-tau/lcal['external_lapse']) < 1e-10
    assert abs(cal['rest_mass_energy_readback_MeV']-pub['rest_energy_MeV']) < 1e-12
    assert np.max(abs(np.asarray(cal['static_gravity_response_MeV'])
                      -pub['rest_energy_MeV']*np.eye(2))) < 1e-12
    checks['calibration_and_constant_lapse_identities'] = True
    minima = {}
    for name in ('Kd', 'H', 'Kx'):
        matrix = data[name]
        assert np.max(abs(matrix-matrix.T)) < 1e-12
        minima[name] = float(np.linalg.eigvalsh(matrix)[0])
        assert minima[name] > 0
    assert data['eta'] > 0 and data['tau'] > 0
    checks['positive_matrix_minima'] = minima
    checks['routh_hurwitz_margin_eta'] = data['eta']
    # For every real k: H+k^2 Kx is positive definite. With positive Kd,
    # each generalized eigenvalue m is positive. The modal polynomial is
    # tau*s^3+s^2+(tau*m+eta)*s+m, whose Routh determinant is eta>0.
    checks['all_real_k_linear_stability_condition'] = True
    ode = []
    eye = np.eye(4)
    zero = np.zeros((4, 4))
    y0 = np.r_[np.array([.1, .2, -.1, .05]), np.zeros(8)]
    for k in (0., .1, 1., 5.):
        A = np.block([[zero, eye, zero],
                      [-np.linalg.solve(data['Kd'], data['H']+k*k*data['Kx']),
                       zero, -data['eta']*eye],
                      [zero, eye/data['tau'], -eye/data['tau']]])
        sol = solve_ivp(lambda t, y: A@y, (0, 2), y0,
                        t_eval=[0., 1., 2.], rtol=1e-10, atol=1e-12)
        exact = np.column_stack([expm(A*t)@y0 for t in sol.t])
        err = float(np.max(abs(sol.y-exact)))
        assert sol.success and err < 1e-8
        ode.append({'k': k, 'maximum_ode_exponential_difference': err})
    checks['time_evolution_crosscheck'] = ode
    result = {
        'status': 'PASS_NUMERICAL_REPRODUCTION_AND_LINEAR_STABILITY',
        'checks': checks,
        'input_roles': 'Public mass and lifetime calibrate the energy map; not independent test targets.',
        'scope': 'Supplied four-mode linear spatial-memory response and seven-coordinate embedding.',
        'limitations': baseline['not_implemented'] + [
            'Finite probe initial condition is numerical, not measured neutron state.',
            'Constant lapse is externally specified, not a measured or solved gravitational field.',
            'Matrix positivity and stability are numerical evidence plus an analytic conditional argument; no interval certificate.',
            'No nonlinear all-field stability or independent experimental validation is established.'],
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('run.py', 'inputs.json', 'verify.py', 'reference_results.json')},
    }
    out = Path(args.output)
    if not out.is_absolute():
        out = ROOT/out
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(result['status'])
    print('Verification:', out.name)


if __name__ == '__main__':
    main()
