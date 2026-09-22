"""
Extended lobule model: transporter zonation + CYP metabolic zonation, for a
joint transporter-CYP hepatic DDI model.

N compartments PV(1) -> CV(N), a sinusoidal blood-flow/transport backbone
with an added metabolic sink in each hepatocyte compartment producing a
tracked metabolite pool.

N=15: Ruijter et al. (2004, Hepatology 39:343-352) directly report 13-15
cells along the centroportal (PV->CV) axis in rat liver, via stereological
mapping (compare 16-18 in human, 9-10 in mouse); N=15 sits at the top of
that reported range. Supersedes an earlier N=24, an uncited round-number
discretization used in an earlier draft of this model.

Pure Python, no dependencies (fast, dt~1e-3, t~20 with explicit Euler).
"""

N = 15
xs = [(i + 0.5) / N for i in range(N)]

# ---- shared / transporter parameters ----
KM = 8.2        # uM, Oatp1a4-type transporter Km (Akanuma et al. 2019)
KI = 0.3        # uM, digoxin (transporter inhibitor) Ki
F = 8.0         # blood flow, a.u./min
VB = 0.15       # fractional blood volume / compartment
VH = 0.6        # fractional hepatocyte volume / compartment
KEF = 0.8       # /min, reversible transporter efflux (hepatocyte -> blood)
VMAX0 = 10.0    # total transporter capacity, a.u.

# ---- CYP / metabolism parameters ----
KM_CYP = 5.0    # uM, illustrative CYP Km
KI_CYP = 0.5    # uM, illustrative CYP inhibitor (e.g. ketoconazole-like) Ki
KEFM = 0.8      # /min, metabolite hepatocyte -> blood permeation (passive, one-way)


def shape(kind):
    """Unnormalized spatial shape, PV->CV, shared functional forms for both
    transporter and CYP zonation hypotheses."""
    if kind == "pericentral":     # sharp CV-biased (Akanuma Oatp1a4 IHC; rat Cyp3a IHC)
        raw = [0.04 + 0.96 * x**3 for x in xs]
    elif kind == "flat":          # non-zonated
        raw = [1.0 for _ in xs]
    elif kind == "periportal":    # mirror-image, PV-biased
        raw = [0.04 + 0.96 * (1 - x)**3 for x in xs]
    else:
        raise ValueError(kind)
    m = sum(raw) / len(raw)
    return [r / m for r in raw]


def vmax_transporter(kind, vmax0=VMAX0):
    return [s * vmax0 for s in shape(kind)]


def cypmax_profile(kind, fm=0.5, kef=KEF, km_cyp=KM_CYP):
    """fm-parameterized total CYP capacity: fm = kmet_avg / (kmet_avg + kef)
    is the fraction of hepatocyte elimination going through metabolism
    (vs. reversible transporter efflux) at the average (non-saturating)
    condition, matching real fm values reportable for a drug without
    needing a fully-fitted kinetic dataset."""
    kmet_avg = fm * kef / (1 - fm)
    s = shape(kind)
    return [si * kmet_avg * km_cyp for si in s]  # since mean(s)=1, CYPmax_i = s_i * kmet_avg * Km_cyp


def simulate(transporter_kind, cyp_kind, fm=0.5, digoxin_uM=0.0, cyp_inh_uM=0.0,
             dose_uM=10.0, t_end_min=10.0, dt=0.0025, vmax0=VMAX0):
    """Integrate parent (Cb, Ch) and metabolite (Cmb, Cm) compartments forward
    by explicit Euler. Returns dict of the four length-N arrays at t_end."""
    Vmax = vmax_transporter(transporter_kind, vmax0)
    CYPmax = cypmax_profile(cyp_kind, fm=fm)
    Km_app = KM * (1 + digoxin_uM / KI) if digoxin_uM > 0 else KM
    Km_cyp_app = KM_CYP * (1 + cyp_inh_uM / KI_CYP) if cyp_inh_uM > 0 else KM_CYP

    Cb = [0.0] * N
    Ch = [0.0] * N
    Cmb = [0.0] * N   # metabolite, blood
    Cm = [0.0] * N    # metabolite, hepatocyte
    nsteps = int(t_end_min / dt)

    for _ in range(nsteps):
        Cb_in = [dose_uM] + Cb[:-1]
        Cmb_in = [0.0] + Cmb[:-1]
        nCb, nCh, nCmb, nCm = [0.0]*N, [0.0]*N, [0.0]*N, [0.0]*N
        for i in range(N):
            Jin = Vmax[i] * Cb[i] / (Km_app + Cb[i])
            Jout = KEF * Ch[i]
            Jmet = CYPmax[i] * Ch[i] / (Km_cyp_app + Ch[i])
            Jmetb = KEFM * Cm[i]

            dCb = F/VB*(Cb_in[i] - Cb[i]) - (Jin - Jout)/VB
            dCh = (Jin - Jout - Jmet)/VH
            dCmb = F/VB*(Cmb_in[i] - Cmb[i]) + Jmetb/VB
            dCm = (Jmet - Jmetb)/VH

            nCb[i] = max(0.0, Cb[i] + dt*dCb)
            nCh[i] = max(0.0, Ch[i] + dt*dCh)
            nCmb[i] = max(0.0, Cmb[i] + dt*dCmb)
            nCm[i] = max(0.0, Cm[i] + dt*dCm)
        Cb, Ch, Cmb, Cm = nCb, nCh, nCmb, nCm

    return {"Cb": Cb, "Ch": Ch, "Cmb": Cmb, "Cm": Cm}


def cv_pv_ratio(arr, n_edge=3):
    pv = sum(arr[:n_edge]) / n_edge
    cv = sum(arr[-n_edge:]) / n_edge
    return cv / max(pv, 1e-9)


if __name__ == "__main__":
    transporters = ["pericentral", "flat", "periportal"]
    cyps = ["pericentral", "flat", "periportal"]

    print("=== Baseline grid: parent (Ch) and metabolite (Cm) CV:PV, t=10 min, fm=0.5 ===")
    print(f"{'transporter':>12} {'CYP':>12}   Ch CV:PV   Cm CV:PV")
    for tr in transporters:
        for cy in cyps:
            r = simulate(tr, cy, fm=0.5, t_end_min=10.0)
            print(f"{tr:>12} {cy:>12}   {cv_pv_ratio(r['Ch']):7.3f}   {cv_pv_ratio(r['Cm']):7.3f}")
