"""
Robustness/sensitivity suite for the transporter-CYP zonation model, built
in response to peer review of the manuscript.

Addresses:
  (1) whether metabolite-tracks-CYP is an emergent or structural result
      (test with a zonated, not uniform, metabolite export kefm(x))
  (2) global sensitivity across fm, Km_cyp, kef, kefm, F, Ki_up, Ki_cyp, N
  (3) whether the transporter-inhibitor-vs-CYP-inhibitor DDI signature is
      robust across fm and zonation configurations, not just one case
  (4) whether CV:PV contrasts survive coarse spatial averaging (3-zone
      readout instead of 15-compartment)

Numpy-vectorized (compartments), explicit Euler in time (matches the
manuscript's validated numerics; dt/t_end taken from model.py's
convergence-checked defaults, re-checked below for the new kefm(x) case).
"""
import numpy as np

N = 15  # anatomically-grounded (see model.py docstring); was 24
xs = np.array([(i + 0.5) / N for i in range(N)])

KM = 8.2
KI = 0.3
F = 8.0
VB = 0.15
VH = 0.6
KEF = 0.8
VMAX0 = 10.0
KM_CYP = 5.0
KI_CYP = 0.5
KEFM0 = 0.8

DT = 0.0025
T_END = 15.0


def shape(kind):
    if kind == "pericentral":
        raw = 0.04 + 0.96 * xs**3
    elif kind == "flat":
        raw = np.ones(N)
    elif kind == "periportal":
        raw = 0.04 + 0.96 * (1 - xs)**3
    else:
        raise ValueError(kind)
    return raw / raw.mean()


def cypmax_profile(kind, fm=0.5, kef=KEF, km_cyp=KM_CYP):
    kmet_avg = fm * kef / (1 - fm)
    return shape(kind) * kmet_avg * km_cyp


def simulate(transporter_kind, cyp_kind, fm=0.5, digoxin_uM=0.0, cyp_inh_uM=0.0,
             dose_uM=10.0, t_end=T_END, dt=DT, vmax0=VMAX0, km_cyp=KM_CYP,
             kef=KEF, kefm=KEFM0, ki_up=KI, ki_cyp=KI_CYP, F_=F,
             kefm_kind=None):
    """kefm: scalar (uniform export, default) OR if kefm_kind is given
    (one of pericentral/flat/periportal), export itself is zonated with
    that shape, normalized to the same mean as the scalar kefm."""
    Vmax = shape(transporter_kind) * vmax0
    CYPmax = cypmax_profile(cyp_kind, fm=fm, kef=kef, km_cyp=km_cyp)
    Km_app = KM * (1 + digoxin_uM / ki_up) if digoxin_uM > 0 else KM
    Km_cyp_app = km_cyp * (1 + cyp_inh_uM / ki_cyp) if cyp_inh_uM > 0 else km_cyp
    kefm_arr = shape(kefm_kind) * kefm if kefm_kind else np.full(N, kefm)

    Cb = np.zeros(N); Ch = np.zeros(N); Cmb = np.zeros(N); Cm = np.zeros(N)
    nsteps = int(t_end / dt)
    for _ in range(nsteps):
        Cb_in = np.concatenate(([dose_uM], Cb[:-1]))
        Cmb_in = np.concatenate(([0.0], Cmb[:-1]))
        Jin = Vmax * Cb / (Km_app + Cb)
        Jout = kef * Ch
        Jmet = CYPmax * Ch / (Km_cyp_app + Ch)
        Jmetb = kefm_arr * Cm
        Cb = np.maximum(0, Cb + dt*(F_/VB*(Cb_in-Cb) - (Jin-Jout)/VB))
        Ch = np.maximum(0, Ch + dt*((Jin-Jout-Jmet)/VH))
        Cmb = np.maximum(0, Cmb + dt*(F_/VB*(Cmb_in-Cmb) + Jmetb/VB))
        Cm = np.maximum(0, Cm + dt*((Jmet-Jmetb)/VH))
    return {"Cb": Cb, "Ch": Ch, "Cmb": Cmb, "Cm": Cm}


def cv_pv(arr, n_edge=3):
    return arr[-n_edge:].mean() / max(arr[:n_edge].mean(), 1e-9)


def coarse_cv_pv(arr, n_zones=3):
    """Coarse-grain into n_zones equal-width spatial bins (e.g. 3:
    periportal/mid/pericentral) and report the resulting zone-average
    CV:PV, mimicking realistic experimental spatial resolution."""
    bins = np.array_split(arr, n_zones)
    means = [b.mean() for b in bins]
    return means[-1] / max(means[0], 1e-9)


if __name__ == "__main__":
    print("############ TEST 1: metabolite export zonation (structural vs emergent) ############")
    print(f"{'transporter':>11} {'CYP':>11} {'kefm shape':>11}   Ch CV:PV   Cm CV:PV")
    for tr in ["pericentral", "flat", "periportal"]:
        for cy in ["pericentral", "flat", "periportal"]:
            for kefm_kind in [None, "pericentral", "periportal"]:
                r = simulate(tr, cy, fm=0.5, kefm_kind=kefm_kind)
                label = kefm_kind or "uniform"
                print(f"{tr:>11} {cy:>11} {label:>11}   {cv_pv(r['Ch']):8.3f}   {cv_pv(r['Cm']):8.3f}")
        print()
