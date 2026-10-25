"""Nam Seung-ho law: stored spatial-memory response and explicit data calibration.

This is a numerical implementation of available response blocks, not a recovered
universal neutron functional. Mass/lifetime below are calibration inputs.
"""
from pathlib import Path
import argparse, datetime, hashlib, json, platform
import numpy as np
import scipy
from scipy.linalg import eigh, eigvals

ROOT = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load():
    data = json.loads((ROOT / 'inputs.json').read_text())
    for name in ['Kd', 'H', 'Kx', 'V', 'A', 'Q']:
        data[name] = np.asarray(data[name], dtype=float)
    return data

def response(data, w, k):
    """Dimensionless retarded response in the supplied four-coordinate basis."""
    eta, tau = data['eta'], data['tau']
    return (data['H'] + k*k*data['Kx'] - w*w*data['Kd']
            + (-1j*eta*w/(1-1j*tau*w))*data['Kd'])

def poles(data, k):
    """Direct augmented ODE: z'=v, v'=-K^-1(H+k²X)z-eta*r,
    tau*r'=v-r. Fourier convention exp(-i*w*t); w=i*lambda."""
    z, eye = np.zeros((4, 4)), np.eye(4)
    eta, tau = data['eta'], data['tau']
    dyn = np.block([[z, eye, z],
                    [-np.linalg.solve(data['Kd'], data['H']+k*k*data['Kx']), z, -eta*eye],
                    [z, eye/tau, -eye/tau]])
    return 1j*eigvals(dyn)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lapse', type=float, default=1.0,
                        help='External constant static lapse N>0, not a solved metric')
    parser.add_argument('--output', default='results.json')
    args = parser.parse_args()
    if not np.isfinite(args.lapse) or args.lapse <= 0:
        parser.error('--lapse must be finite and positive')
    d = load(); K,H,X,V,A,Q = (d[n] for n in ['Kd','H','Kx','V','A','Q'])
    checks = {
        'K_positive_min_eigenvalue': float(eigh(K,eigvals_only=True)[0]),
        'X_positive_min_generalized_eigenvalue': float(eigh(X,K,eigvals_only=True)[0]),
        'normalization_max_error': float(np.max(abs(V.T@K@V-np.eye(4)))),
        'H_diagonalization_max_error': float(np.max(abs(V.T@H@V-np.diag(d['mu_squared'])))),
        'Q_embedding_max_error': float(np.max(abs(Q-np.vstack([V,A@V])))),
        'Q_rank': int(np.linalg.matrix_rank(Q)),
    }
    assert checks['K_positive_min_eigenvalue'] > 0
    assert checks['X_positive_min_generalized_eigenvalue'] > 0
    assert checks['normalization_max_error'] < 1e-12
    assert checks['H_diagonalization_max_error'] < 1e-11
    assert checks['Q_embedding_max_error'] < 1e-12
    samples=[]
    for k in [0., .1, 1., 5.]:
        direct=poles(d,k)
        mus=eigh(H+k*k*X,K,eigvals_only=True)
        reference=np.concatenate([np.roots([1j*d['tau'],-1,
            -1j*(d['tau']*m+d['eta']),m]) for m in mus])
        error=float(max(min(abs(w-r) for r in reference) for w in direct))
        assert error<1e-10
        assert max(w.imag for w in direct)<0
        # Relative smallest singular value of the original matrix at each pole.
        residuals=[]
        for w in direct:
            D=response(d,w,k)
            residuals.append(float(np.linalg.svd(D,compute_uv=False)[-1]/max(1,np.linalg.norm(D,2))))
        # Check the actual Schur response in the stored coordinate basis.
        D=response(d,1+.2j,k)
        S=D[:1,:1]-D[:1,1:]@np.linalg.solve(D[1:,1:],D[1:,:1])
        schur=float(abs(np.linalg.det(D)-np.linalg.det(D[1:,1:])*S[0,0]))
        assert schur<1e-9
        samples.append({'k_dimensionless':k,'poles':[[float(w.real),float(w.imag)]
            for w in sorted(direct,key=lambda w:(w.real,w.imag))],
            'matrix_vs_cubic_max_difference':error,
            'matrix_pole_max_relative_residual':max(residuals),
            'Schur_identity_absolute_error':schur})
    # New phenomenological identification: smallest positive-frequency k=0 mode.
    # Neither the neutron quantum numbers nor this identification are derived.
    physical_candidates=[w for w in poles(d,0.) if w.real>1e-8]
    w0=min(physical_candidates,key=lambda w:w.real)
    pub=d['public_inputs']; hb=pub['hbar_MeV_s']; mass=pub['rest_energy_MeV']; lifetime=pub['lifetime_s']
    scale=hb/(2*lifetime*(-w0.imag))
    # Keep tiny relative complex energy separate from the rest energy.
    delta=scale*(w0-w0.real)
    N=args.lapse
    calibration={
        'classification':'PUBLIC_INPUT_CALIBRATION_NOT_INDEPENDENT_PREDICTION',
        'new_assumptions':['identify lowest positive k=0 response mode with a phenomenological neutron pole',
                           'affine energy map fixed by observed mass and lifetime',
                           'constant externally supplied lapse scales the complete complex pole',
                           'static whole-state coupling Hg=mass_energy*phi'],
        'frequency_origin':float(w0.real), 'energy_scale_MeV':float(scale),
        'rest_mass_energy_readback_MeV':float(mass+delta.real),
        'rest_pole_relative_energy_MeV':[float(delta.real),float(delta.imag)],
        'proper_lifetime_readback_s':float(-hb/(2*delta.imag)),
        'external_lapse':N,'coordinate_lifetime_s':float(-hb/(2*N*delta.imag)),
        'static_gravity_response_MeV':(mass*np.eye(2)).tolist(),
        'survival_at_878_4_coordinate_seconds':float(np.exp(-N*878.4/lifetime)),
        'mass_lifetime_independently_validated':False,
        'metric_solved_from_Gamma_NS':False,
        'spin_charge_mapping_derived':False,
    }
    assert np.isclose(calibration['proper_lifetime_readback_s'],lifetime,rtol=1e-13)
    result={
        'name':'남승호법칙 중성자 공간 기억응답 수치구현',
        'time_KST':datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        'scope':d['scope'],'checks':checks,'spatial_samples':samples,
        'public_input_effective_neutron':calibration,
        'not_implemented':['full 3D topological and binding spatial operator',
                           'weak-transition vertex from same Gamma_NS',
                           'neutron spin and charge sector derivation',
                           'universal metric variation and full gravity field'],
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'sha256':{'run.py':digest(Path(__file__)),'inputs.json':digest(ROOT/'inputs.json')},
        'classification':'IMPLEMENTATION_AND_NUMERICAL_CHECKS_COMPLETE_WITHIN_DECLARED_SCOPE',
        'whole_neutron_theory_validation_complete':False,
    }
    output=Path(args.output)
    if not output.is_absolute(): output=ROOT/output
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('공간 기억응답 검산: PASS')
    print('12개 극 × 4개 k; 성장극: 0')
    print('공개입력 기반 질량에너지:',calibration['rest_mass_energy_readback_MeV'],'MeV')
    print('공개입력 기반 고유수명:',calibration['proper_lifetime_readback_s'],'s')
    print('중력장은 외부입력; 전체 중성자 이론검증 완료 아님')
    print('결과:',output.name)

if __name__=='__main__': main()
