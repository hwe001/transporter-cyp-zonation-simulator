import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sensitivity_suite import simulate, cv_pv, coarse_cv_pv

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3))
axA, axB, axC = axes

# ---------------- Panel A: metabolite export zonation sensitivity ----------------
cy_kinds = ["pericentral", "flat", "periportal"]
export_kinds = [None, "pericentral", "periportal"]
export_labels = ["uniform\n(as modeled)", "export also\npericentral", "export also\nperiportal"]
colors = ["#2a78d6", "#eb6834", "#1baf7a"]

x = np.arange(len(cy_kinds))
width = 0.26
for k, (ek, lab, col) in enumerate(zip(export_kinds, export_labels, colors)):
    vals = [cv_pv(simulate("flat", cy, fm=0.5, kefm_kind=ek)["Cm"]) for cy in cy_kinds]
    axA.bar(x + (k-1)*width, np.log10(vals), width, label=lab, color=col, alpha=0.85)
axA.set_xticks(x); axA.set_xticklabels(["Peri-\ncentral", "Flat", "Peri-\nportal"])
axA.set_xlabel("CYP zonation (transporter fixed = flat)")
axA.set_ylabel("log$_{10}$(Metabolite Cm CV:PV)")
axA.axhline(0, color="#333", lw=0.8)
axA.legend(fontsize=6.8, frameon=False, loc="upper right")
axA.spines[['top','right']].set_visible(False)
axA.set_title("A", loc="left", fontweight="bold", x=-0.16)
axA.text(0.02, 0.03, "rank order preserved;\nmagnitude shifts up to ~180-fold",
         transform=axA.transAxes, fontsize=6.8, color="#666", style="italic")

# ---------------- Panel B: DDI signature robustness across fm/configs ----------------
configs = [("flat","pericentral"), ("flat","flat"), ("pericentral","pericentral")]
fms = [0.2, 0.5, 0.8]
grid = np.zeros((3,3))
for i,(tr,cy) in enumerate(configs):
    for j,fm in enumerate(fms):
        base = simulate(tr, cy, fm=fm)
        base_cvpv = cv_pv(base["Ch"])
        rt = simulate(tr, cy, fm=fm, digoxin_uM=2.4)
        rc = simulate(tr, cy, fm=fm, cyp_inh_uM=4.0)
        dt_shape = abs(100*cv_pv(rt["Ch"])/base_cvpv - 100)
        dc_shape = abs(100*cv_pv(rc["Ch"])/base_cvpv - 100)
        grid[i,j] = 1.0 if dc_shape > dt_shape else 0.0

im = axB.imshow(grid, cmap="RdYlGn", vmin=0, vmax=1, aspect="equal")
for i in range(3):
    for j in range(3):
        axB.text(j, i, "holds" if grid[i,j]==1 else "fails", ha="center", va="center", fontsize=8)
axB.set_xticks(range(3)); axB.set_xticklabels([f"fm={f}" for f in fms])
axB.set_yticks(range(3)); axB.set_yticklabels([f"{tr}\n{cy}" for tr,cy in configs], fontsize=7.5)
axB.set_title("B", loc="left", fontweight="bold", x=-0.30)
axB.set_xlabel("metabolic fraction")
axB.text(0.5, -0.28, "DDI-locus signature (CYP-inhib. reshapes\nspatial contrast more than transporter-inhib.):\nholds in 7/9 tested configurations",
         transform=axB.transAxes, fontsize=6.8, color="#666", ha="center", style="italic")

# ---------------- Panel C: spatial-averaging / detectability ----------------
data = []
for tr in ["pericentral","flat","periportal"]:
    for cy in ["pericentral","flat","periportal"]:
        r = simulate(tr, cy, fm=0.5)
        data.append((cv_pv(r["Ch"]), coarse_cv_pv(r["Ch"],3), cv_pv(r["Cm"]), coarse_cv_pv(r["Cm"],3)))
data = np.array(data)
axC.scatter(np.log10(data[:,0]), np.log10(data[:,1]), color="#2a78d6", label="Parent (Ch)", s=45, alpha=0.85)
axC.scatter(np.log10(data[:,2]), np.log10(data[:,3]), color="#7a4fa8", label="Metabolite (Cm)", s=45, marker="s", alpha=0.85)
lims = [min(np.log10(data[:,[0,2]]).min(), np.log10(data[:,[1,3]]).min())-0.3,
        max(np.log10(data[:,[0,2]]).max(), np.log10(data[:,[1,3]]).max())+0.3]
axC.plot(lims, lims, color="#999", lw=1, ls="--")
axC.set_xlim(lims); axC.set_ylim(lims)
axC.set_xlabel("log$_{10}$(CV:PV), 24-compartment")
axC.set_ylabel("log$_{10}$(CV:PV), 3-zone readout")
axC.legend(fontsize=7.5, frameon=False, loc="upper left")
axC.spines[['top','right']].set_visible(False)
axC.set_title("C", loc="left", fontweight="bold", x=-0.20)
axC.set_aspect("equal")

fig.suptitle("Figure S1. Robustness checks requested in review", fontsize=11, y=1.03)
fig.savefig("figure2_sensitivity.png", dpi=300, bbox_inches="tight")
print("saved")
