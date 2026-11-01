"""Independent numerical cross-check of the fixed 1580-row residual comparison.

Adaptive quadrature replaces the primary Gauss-Legendre integration; an analytic
empty-open GR expression checks the fitted boundary. Same environment, not a
third-party reproduction or validation of the physical source of residuals.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import cho_factor,cho_solve

ROOT=Path(__file__).parent
s=json.loads((ROOT/'reference_results/summary.json').read_text())
for name,digest in s['input_sha256'].items():
    if hashlib.sha256((ROOT/'data'/name).read_bytes()).hexdigest()!=digest:
        raise RuntimeError(f'Input checksum mismatch: {name}')
df=pd.read_csv(ROOT/'data/Pantheon+SH0ES.dat',sep=r'\s+')
keep=np.flatnonzero((df.IS_CALIBRATOR.to_numpy()==0)&(df.zHD.to_numpy()>.01))
d=df.iloc[keep];assert len(d)==1580
raw=np.loadtxt(ROOT/'data/Pantheon+SH0ES_STAT+SYS.cov')
N=int(raw[0]);C=raw[1:].reshape(N,N)[np.ix_(keep,keep)];C=(C+C.T)/2
cf=cho_factor(C,lower=True);one=np.ones(len(d));ci1=cho_solve(cf,one);den=one@ci1
a=s['NS']['alpha'];v=s['NS']['nu'];z=d.zHD.to_numpy();zh=d.zHEL.to_numpy();y=d.m_b_corr.to_numpy()
dc=np.array([quad(lambda t:np.exp(-a/2*((1+t)**v-1)),0,float(zi),epsabs=2e-12,epsrel=2e-12)[0] for zi in z])
milne=((1+z)-1/(1+z))/2
audit={}
reference=pd.read_csv(ROOT/'reference_results/row_residuals.csv')
assert np.array_equal(keep,reference.original_row_index.to_numpy())
for name,dist in [('NS_fixed',dc),('GR_noLambda_curvature',milne)]:
    shape=5*np.log10((1+zh)*dist);offset=float(ci1@(y-shape)/den);res=y-shape-offset
    chi=float(res@cho_solve(cf,res));delta=chi-s['models'][name]['chi2']
    row_difference=float(np.max(np.abs(res-reference[name+'_residual_mag'].to_numpy())))
    if abs(delta)>1e-6 or row_difference>1e-8:raise RuntimeError(f'Reproduction mismatch: {name}')
    audit[name]={'chi2':chi,'difference_from_reference':delta,'max_row_residual_difference_mag':row_difference}
    produced=ROOT/'results/row_residuals.csv'
    if produced.exists():
        current=pd.read_csv(produced)
        assert np.array_equal(keep,current.original_row_index.to_numpy())
        assert np.max(np.abs(res-current[name+'_residual_mag'].to_numpy()))<1e-8
nodes,weights=np.polynomial.legendre.leggauss(64);u=z[:,None]*(nodes+1)/2
grid=np.linspace(0,3,61);scores=[]
for om in grid:
    radial=z/2*((1/np.sqrt(om*(1+u)**3+(1-om)*(1+u)**2))@weights)
    ok=1-om
    dist=radial if abs(ok)<1e-12 else (np.sinh(np.sqrt(ok)*radial)/np.sqrt(ok) if ok>0 else np.sin(np.sqrt(-ok)*radial)/np.sqrt(-ok))
    shape=5*np.log10((1+zh)*dist);offset=ci1@(y-shape)/den;res=y-shape-offset
    scores.append(float(res@cho_solve(cf,res)))
assert int(np.argmin(scores))==0
audit['GR_grid_minimum_Omega_m']=float(grid[np.argmin(scores)])
audit['scope']='same-environment numerical cross-check; not third-party physical validation'
(ROOT/'results').mkdir(exist_ok=True)
(ROOT/'results/verification.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2));print('NUMERICAL REPRODUCTION PASS')
