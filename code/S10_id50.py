#!/usr/bin/env python3
"""S10 — ID50: the %identity at which transfer fidelity falls to half of its near-plateau value.
Fits a constrained logistic to (Wang vs %identity) per track x aspect, with bootstrap CI and an
identifiability rule (FIX 1, pre-manuscript hardening). When the logistic is not identifiable
(flat curve / boundary-pinned / wide CI / poor fit) we fall back to the robust interpolation ID50
and flag it. Representative value = ID50_rep.

Inputs: each track's metrics. Outputs: results/crossclade/{id50.tsv, fig_id50.svg/.pdf}
"""
import os
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = "/var2/lsg/Claude_Code/Cross-species-GeneOntology"
CLADE = {"fish":"fish","mammal":"mammal","plant_rice":"plant","plant_arabidopsis":"plant",
         "fungi":"fungi","insect":"insect"}
OUTD = f"{ROOT}/results/crossclade"; os.makedirs(OUTD, exist_ok=True)
NBOOT = 1000
K_UPPER = 5.0           # relaxed (was hard 2.0); boundary flagged via identifiability instead
CI_WIDTH_MAX = 25.0     # %p; wider => not identifiable
R2_MIN = 0.5
RNG = np.random.default_rng(20260616)

def logistic(x, lo, hi, id50, k):
    return lo + (hi-lo)/(1+np.exp(-k*(x-id50)))

def interp_id50(x, y):
    o=np.argsort(x); x,y=np.array(x)[o],np.array(y)[o]
    hi=np.mean(np.sort(y)[-2:]); lo=np.min(y); tgt=(hi+lo)/2
    for i in range(len(x)-1):
        if (y[i]-tgt)*(y[i+1]-tgt)<=0 and y[i+1]!=y[i]:
            return x[i]+(tgt-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]), lo, hi
    return np.nan, lo, hi

def fit_logistic(x, y):
    """constrained logistic fit; returns (lo,hi,id50,k,R2) or None."""
    lo0, hi0 = float(np.min(y)), float(np.mean(np.sort(y)[-2:]))
    xmn, xmx = float(np.min(x)), float(np.max(x))
    # constraints: floor>=0, plateau in [floor, observed max], id50 within observed range, k>0
    bounds = ([0.0, lo0-1e-6, xmn, 1e-3], [max(hi0,lo0)+1e-6, 1.0, xmx, K_UPPER])
    p0 = [lo0, hi0, np.clip(interp_id50(x,y)[0] if not np.isnan(interp_id50(x,y)[0]) else np.median(x), xmn, xmx), 0.2]
    try:
        popt,_ = curve_fit(logistic, x, y, p0=p0, bounds=bounds, maxfev=30000)
    except Exception:
        return None
    yhat = logistic(x, *popt)
    ss_res = float(np.sum((y-yhat)**2)); ss_tot = float(np.sum((y-np.mean(y))**2))
    r2 = 1 - ss_res/ss_tot if ss_tot>0 else np.nan
    return (*popt, r2)

rows=[]; curves={}
for track,clade in CLADE.items():
    f=f"{ROOT}/results/{track}/metrics/semantic_sim.tsv"
    if not os.path.exists(f): continue
    sem=pd.read_csv(f,sep="\t")
    m=pd.read_csv(f"{ROOT}/results/{track}/metrics/setlevel_metrics.tsv",sep="\t")
    pid=m[m.ic_bin=="all"].groupby("ref").median_pident.first()
    sem=sem[sem.list_id.str.startswith("SYN_")].copy()
    sem["ref"]=sem.label.str.split(".").str[0]; sem["method"]=sem.label.str.split(".").str[1]; sem["evset"]=sem.label.str.split(".").str[2]
    sem=sem[(sem.method=="rbh")&(sem.evset=="all")]
    sem["pid"]=sem.ref.map(pid)
    for asp in ["MF","BP","CC"]:
        a=sem[sem.aspect==asp].groupby("ref").agg(wang=("wang_bma","mean"),pid=("pid","first")).dropna()
        if len(a)<5: continue
        x=a.pid.values.astype(float); y=a.wang.values.astype(float)
        xmn,xmx=float(x.min()),float(x.max())
        id50_i, lo, hi = interp_id50(x,y)
        fit = fit_logistic(x,y)
        id50_f = fit[2] if fit else np.nan
        k_f    = fit[3] if fit else np.nan
        r2     = fit[4] if fit else np.nan
        # bootstrap CI on logistic ID50 (resample reference points)
        boot=[]
        if fit:
            n=len(x)
            for _ in range(NBOOT):
                idx=RNG.integers(0,n,n)
                if len(np.unique(x[idx]))<4: continue
                fb=fit_logistic(x[idx],y[idx])
                if fb: boot.append(fb[2])
        ci_lo=ci_hi=np.nan
        if len(boot)>=50:
            ci_lo,ci_hi=np.percentile(boot,[2.5,97.5])
        # identifiability rule
        k_boundary  = (not np.isnan(k_f)) and (k_f>=0.99*K_UPPER or k_f<=1.5e-3)
        id50_pinned = (not np.isnan(id50_f)) and (id50_f<=xmn+1 or id50_f>=xmx-1)  # half-point at/outside data range
        ci_outside  = (not np.isnan(ci_lo)) and (ci_lo<xmn-1 or ci_hi>xmx+1)
        ci_wide     = (not np.isnan(ci_lo)) and ((ci_hi-ci_lo)>CI_WIDTH_MAX)
        poor_fit    = (not np.isnan(r2)) and (r2<R2_MIN)
        no_ci       = np.isnan(ci_lo)
        identifiable = bool(fit) and not (k_boundary or id50_pinned or ci_outside or ci_wide or poor_fit or no_ci)
        if identifiable:
            id50_rep=id50_f; method_used="logistic"
        else:
            id50_rep=id50_i; method_used="interp"
        rows.append(dict(track=track,clade=clade,aspect=asp,n=len(a),
                         plateau=round(hi,3),floor=round(lo,3),
                         ID50_interp=round(id50_i,1) if not np.isnan(id50_i) else np.nan,
                         ID50_logistic=round(id50_f,1) if not np.isnan(id50_f) else np.nan,
                         slope_k=round(k_f,3) if not np.isnan(k_f) else np.nan,
                         fit_R2=round(r2,3) if not np.isnan(r2) else np.nan,
                         ID50_CI_low=round(ci_lo,1) if not np.isnan(ci_lo) else np.nan,
                         ID50_CI_high=round(ci_hi,1) if not np.isnan(ci_hi) else np.nan,
                         identifiable=identifiable, method_used=method_used,
                         ID50_rep=round(id50_rep,1) if not np.isnan(id50_rep) else np.nan))
        curves[(track,asp)]=(x,y,id50_f,k_f,lo,hi,identifiable,ci_lo,ci_hi)
res=pd.DataFrame(rows); res.to_csv(f"{OUTD}/id50.tsv",sep="\t",index=False)
pd.set_option("display.width",220,"display.max_columns",30)
print(res.to_string(index=False))
print(f"\nidentifiable logistic fits: {res.identifiable.sum()}/{len(res)}; "
      f"rest fall back to interpolated ID50 (method_used).")

# figure: BP fits per track, CI band, non-identifiable dashed/grey
fig,ax=plt.subplots(figsize=(9.5,5.8))
col={"fish":"#1b9e77","mammal":"#d95f02","plant":"#7570b3","fungi":"#e7298a","insect":"#66a61e"}
for (track,asp),(x,y,id50,k,lo,hi,ident,clo,chi) in curves.items():
    if asp!="BP": continue
    c=col[CLADE[track]]; ax.scatter(x,y,color=c,s=24,alpha=.75)
    if not np.isnan(id50):
        xs=np.linspace(min(x)-3,max(x)+3,100)
        ax.plot(xs,logistic(xs,lo,hi,id50,k),color=c,alpha=.55,
                ls="-" if ident else "--", lw=1.8 if ident else 1.2)
        if ident and not np.isnan(clo):
            ax.axvspan(clo,chi,color=c,alpha=.08)          # CI band
            ax.axvline(id50,color=c,ls=":",alpha=.5)
    ax.annotate(track+("" if ident else "*"),(x[np.argmin(x)],y[np.argmin(x)]),fontsize=7,color=c)
ax.invert_xaxis(); ax.set_xlabel("median ortholog % identity"); ax.set_ylabel("Wang semantic similarity (BP)")
ax.set_title("ID50 — half-collapse %identity per track (BP, RBH; solid=identifiable, dashed*=NE→interp; band=95% CI)")
ax.grid(alpha=.3); fig.tight_layout()
fig.savefig(f"{OUTD}/fig_id50.svg"); fig.savefig(f"{OUTD}/fig_id50.pdf")
print(f"\n[S10] wrote id50.tsv + fig_id50.svg/.pdf -> {OUTD}")
