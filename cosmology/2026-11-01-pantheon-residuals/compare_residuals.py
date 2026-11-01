import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,hashlib,platform,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve
from scipy.special import roots_legendre
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
df=pd.read_csv(ROOT/'data/Pantheon+SH0ES.dat',sep=r'\s+')
raw=np.loadtxt(ROOT/'data/Pantheon+SH0ES_STAT+SYS.cov')
N=int(raw[0]); full=raw[1:].reshape(N,N)
mask=(df.zHD.to_numpy()>.01)&(df.IS_CALIBRATOR.to_numpy()==0)
indices=np.flatnonzero(mask); d=df.loc[mask].copy()
C0=full[np.ix_(indices,indices)]
asym=float(np.max(np.abs(C0-C0.T)))
C=(C0+C0.T)/2
cf=cho_factor(C,lower=True)
one=np.ones(len(d)); ci1=cho_solve(cf,one); denom=one@ci1
z=d.zHD.to_numpy(); zh=d.zHEL.to_numpy(); y=d.m_b_corr.to_numpy()
alpha=1.4835435161693789; nu=.8263146605293827

def distances(kind,om=None,order=96):
    nodes,weights=roots_legendre(order)
    u=z[:,None]*(nodes[None,:]+1)/2
    if kind=='NS': invE=np.exp(-alpha/2*((1+u)**nu-1))
    else: invE=1/np.sqrt(om*(1+u)**3+(1-om)*(1+u)**2)
    dc=z/2*(invE@weights)
    if kind!='NS':
        ok=1-om
        if ok>1e-10: dc=np.sinh(np.sqrt(ok)*dc)/np.sqrt(ok)
        elif ok< -1e-10: dc=np.sin(np.sqrt(-ok)*dc)/np.sqrt(-ok)
    return (1+zh)*dc

def evaluate(kind,om=None,order=96):
    shape=5*np.log10(distances(kind,om,order))
    offset=float(ci1@(y-shape)/denom)
    pred=shape+offset; r=y-pred
    chi=float(r@cho_solve(cf,r))
    return dict(shape=shape,offset=offset,pred=pred,r=r,chi2=chi,rms=float(np.sqrt(np.mean(r*r))))

fit=minimize_scalar(lambda om:evaluate('GR',om)['chi2'],bounds=(0,3),method='bounded',options={'xatol':1e-10})
models={'NS_fixed':evaluate('NS'),'GR_noLambda_curvature':evaluate('GR',fit.x)}
checks={}
for key,kind,om in [('NS_fixed','NS',None),('GR_noLambda_curvature','GR',fit.x)]:
    hi=evaluate(kind,om,192)
    checks[key]=dict(max_magnitude_change_96_to_192=float(np.max(np.abs(hi['pred']-models[key]['pred']))),chi2_change=hi['chi2']-models[key]['chi2'],weighted_offset_residual=float(ci1@models[key]['r']))

bins=np.array([.01,.025,.05,.1,.2,.35,.5,.7,1.,1.4,2.3])
which=np.digitize(z,bins)-1
X=np.column_stack([(which==j).astype(float) for j in range(len(bins)-1)])
assert np.all(X.sum(axis=1)==1)
ciX=cho_solve(cf,X); Bcov=np.linalg.inv(X.T@ciX)
W=Bcov@ciX.T
bin_cov=Bcov-np.ones_like(Bcov)/denom
be=np.sqrt(np.maximum(np.diag(bin_cov),0))
binned=pd.DataFrame({'z_low':bins[:-1],'z_high':bins[1:],'n':X.sum(axis=0).astype(int),'z_plot':np.sqrt(bins[:-1]*bins[1:]),'conditional_sigma_mag':be})
for key,m in models.items():
    binned[key+'_residual_mag']=W@m['r']
    d[key+'_predicted_m_b']=m['pred'];d[key+'_residual_mag']=m['r']
d['original_row_index']=indices
d['sigma_marginal_mag']=np.sqrt(np.diag(C))
d[['original_row_index','CID','IDSURVEY','zHD','zHEL','m_b_corr','sigma_marginal_mag','NS_fixed_predicted_m_b','NS_fixed_residual_mag','GR_noLambda_curvature_predicted_m_b','GR_noLambda_curvature_residual_mag']].to_csv(OUT/'row_residuals.csv',index=False)
binned.to_csv(OUT/'binned_residuals.csv',index=False)
np.savetxt(OUT/'binned_residual_covariance.txt',bin_cov)

summary={'run_date_KST':'2026-11-01','classification':'in_sample_fixed_NS_background_residual_comparison_not_independent_prediction','selection':'IS_CALIBRATOR == 0 and zHD > 0.01; retain all selected light-curve rows and corresponding covariance','input_rows':N,'selected_rows':len(d),'selected_unique_CID':int(d.CID.nunique()),'covariance':'full STAT+SYS selected submatrix; symmetrized as (C+C.T)/2','maximum_input_asymmetry':asym,'observable':'m_b_corr; one additive magnitude offset profiled separately for every model; no H0 determination','distance':'dL_dimensionless=(1+zHEL)*S_k(integral_0^zHD du/E(u))','NS':{'alpha':alpha,'nu':nu,'Lambda':0,'coefficients_refit_this_run':False},'GR_noLambda':{'Omega_m':float(fit.x),'Omega_k':float(1-fit.x),'bounds_Omega_m':[0,3],'optimizer_success':bool(fit.success)},'models':{k:{field:m[field] for field in ['offset','chi2','rms']} for k,m in models.items()},'delta_chi2_NS_minus_GR':models['NS_fixed']['chi2']-models['GR_noLambda_curvature']['chi2'],'numerical_checks':checks,'binning':'fixed display edges; joint generalized least squares bin means with full covariance; conditional covariance after intercept removal; GR shape uncertainty not propagated','not_established':['A-path N=1701 saved-table exact reproduction','permanent physical residual or residual-memory causal identification','blind target-excluded independence','third-party reproduction'],'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'data').iterdir() if p.is_file()},'environment':{'python':platform.python_version(),'numpy':np.__version__}}
(OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))

plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True,gridspec_kw={'height_ratios':[1.4,1]},layout='constrained')
colors={'NS_fixed':'#146b80','GR_noLambda_curvature':'#b95328'}
labels={'NS_fixed':'NS: fixed alpha, nu; no Lambda','GR_noLambda_curvature':'GR: fitted matter + curvature; no Lambda'}
for key,m in models.items():
    axes[0].scatter(z,m['r'],s=5,alpha=.12,color=colors[key],rasterized=True)
    axes[1].errorbar(binned.z_plot,binned[key+'_residual_mag'],yerr=be,marker='o',capsize=3,lw=1.6,color=colors[key],label=labels[key])
for ax in axes:
    ax.axhline(0,color='#555555',lw=1,ls='--');ax.set_xscale('log');ax.set_ylabel('Observed - model (mag)');ax.grid(alpha=.18)
axes[0].set_ylim(-.85,.85)
axes[0].set_title('Pantheon+ residual comparison | 1,580 selected light curves',loc='left',fontweight='bold')
axes[0].text(.02,.04,'Full covariance used for fitting; scatter is a visual display.\nTop panel clips extreme residuals; every row remains in the fit and CSV.',transform=axes[0].transAxes,fontsize=9)
axes[1].legend(loc='best',fontsize=9);axes[1].set_xlabel('Redshift zHD')
fig.savefig(OUT/'residual_comparison.png',dpi=170)

report=f'''# Pantheon+ residual comparison — 2026-11-01\n\nThis is an in-sample comparison, not a blind prediction or full Gamma_NS derivation.\n\nSelection: noncalibrators with zHD > 0.01; {len(d)} light curves, {d.CID.nunique()} distinct CID. All repeated observations are retained with the supplied full covariance. Observed magnitude is m_b_corr. H0 and absolute magnitude are absorbed into one profiled additive offset for each model. The luminosity-distance prefactor uses zHEL and the radial integral uses zHD.\n\nNS: E^2=exp(alpha*((1+z)^nu-1)), alpha={alpha}, nu={nu}; no refit of alpha/nu. These coefficients have historical Pantheon links.\n\nGR comparator: E^2=Omega_m*(1+z)^3+(1-Omega_m)*(1+z)^2, Lambda=0, curvature distance S_k retained. Omega_m={fit.x:.12g}; optimized on [0,3].\n\n| Model | chi-square, full covariance | ordinary RMS, mag |\n|---|---:|---:|\n| NS fixed | {models['NS_fixed']['chi2']:.9f} | {models['NS_fixed']['rms']:.9f} |\n| GR no Lambda | {models['GR_noLambda_curvature']['chi2']:.9f} | {models['GR_noLambda_curvature']['rms']:.9f} |\n\nDelta chi-square (NS minus GR): {summary['delta_chi2_NS_minus_GR']:.9f}. RMS is descriptive and does not replace the full-covariance statistic.\n\nThis 1,580-row run is separate from the historic A-path 1,701-row saved table. It does not replace or exactly reproduce that table. No new AIC/BIC claim is made.\n\nBinned points are joint generalized least-squares means at fixed display bins. Their covariance includes removal of the fitted additive offset but conditions on model shape; uncertainty in fitted GR Omega_m is not propagated. The number of display bins is not a physical periodicity claim. Top-panel extreme residuals are clipped only for display.\n\nNo period, permanent memory, or physical origin is inferred from a residual plot alone. Full Gamma_NS-to-cosmology and third-party reproduction remain outside this calculation.\n\nReproduce: place Pantheon+SH0ES.dat and Pantheon+SH0ES_STAT+SYS.cov in data/ next to compare_residuals.py, then run python compare_residuals.py. Dependencies: numpy, pandas, scipy, matplotlib. Input SHA256 hashes and numerical checks are in results/summary.json.\n\nSources: stored project whole-gravity revision 2026-09-27, section 7; official release https://github.com/PantheonPlusSH0ES/DataRelease ; methodology https://arxiv.org/abs/2202.04077 .\n'''
(OUT/'README.md').write_text(report)
print(json.dumps(summary,ensure_ascii=False,indent=2))
