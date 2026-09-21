import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from model import simulate, cv_pv_ratio, N, xs

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})

fig = plt.figure(figsize=(14.0, 8.6))
gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.15], wspace=0.9, hspace=0.42)
axA = fig.add_subplot(gs[0, 0:2])
axB = fig.add_subplot(gs[0, 2:4])
axC = fig.add_subplot(gs[0, 4:6])
axD = fig.add_subplot(gs[1, 0:3])
axE = fig.add_subplot(gs[1, 3:6])

# ---------------- Panel A: schematic ----------------
axA.set_xlim(0, 10); axA.set_ylim(-0.3, 4.3); axA.axis("off")
axA.set_title("A", loc="left", fontweight="bold", x=-0.06)
n_boxes = 5
box_w, box_h = 1.55, 1.35
gap = 0.32
x0 = 0.3
for k in range(n_boxes):
    x = x0 + k*(box_w+gap)
    axA.add_patch(Rectangle((x, 1.55), box_w, box_h, facecolor="#f2ded3", edgecolor="#8a5a3c", lw=1.2, zorder=2))
    axA.text(x+box_w/2, 2.22, f"Ch$_{{{k+1}}}$", ha="center", va="center", fontsize=8.5, color="#5a3a28")
    axA.add_patch(Rectangle((x, 0.75), box_w, 0.55, facecolor="#dbe8f7", edgecolor="#2a5c94", lw=1.0, zorder=2))
    axA.text(x+box_w/2, 1.025, f"Cb$_{{{k+1}}}$", ha="center", va="center", fontsize=8, color="#1c3a5e")
    axA.add_patch(Rectangle((x, -0.05), box_w, 0.5, facecolor="#e3ddf2", edgecolor="#5a3c8a", lw=1.0, zorder=2))
    axA.text(x+box_w/2, 0.2, f"Cm$_{{{k+1}}}$", ha="center", va="center", fontsize=7.5, color="#3a2860")
    axA.annotate("", xy=(x+box_w/2, 1.55), xytext=(x+box_w/2, 1.3),
                 arrowprops=dict(arrowstyle="-|>", color="#555", lw=1.1))
    axA.annotate("", xy=(x+box_w*0.7, -0.05), xytext=(x+box_w*0.7, 0.35),
                 arrowprops=dict(arrowstyle="-|>", color="#7a4fa8", lw=1.2))
    if k < n_boxes-1:
        axA.annotate("", xy=(x+box_w+gap, 1.025), xytext=(x+box_w, 1.025),
                     arrowprops=dict(arrowstyle="-|>", color="#2a5c94", lw=1.4))
        axA.annotate("", xy=(x+box_w+gap, 0.2), xytext=(x+box_w, 0.2),
                     arrowprops=dict(arrowstyle="-|>", color="#7a4fa8", lw=1.3))
last_x = x0 + (n_boxes-1)*(box_w+gap)
axA.annotate("", xy=(last_x+box_w+0.45, 1.025), xytext=(last_x+box_w, 1.025),
             arrowprops=dict(arrowstyle="-|>", color="#2a5c94", lw=1.4))
axA.annotate("", xy=(last_x+box_w+0.45, 0.2), xytext=(last_x+box_w, 0.2),
             arrowprops=dict(arrowstyle="-|>", color="#7a4fa8", lw=1.3))
axA.text(0.3, 3.55, "PV", fontsize=10, fontweight="bold", color="#4c7a52")
axA.text(last_x+box_w-0.3, 3.55, "CV", fontsize=10, fontweight="bold", color="#a8503f")
axA.text(x0+box_w/2, 3.0, "transporter\nuptake/efflux", fontsize=6.8, ha="center", color="#5a3a28")
axA.text(x0+box_w/2, -0.28, "CYP\nmetabolism (sink)", fontsize=6.8, ha="center", color="#5a3c8a")
axA.set_xlim(0, last_x+box_w+1.0)

# ---------------- Panel B: baseline 3x3 grid, parent CV:PV ----------------
transporters = ["pericentral", "flat", "periportal"]
cyps = ["pericentral", "flat", "periportal"]
tlabels = ["Peri-\ncentral", "Flat", "Peri-\nportal"]

grid_Ch = np.zeros((3, 3))
grid_Cm = np.zeros((3, 3))
for i, tr in enumerate(transporters):
    for j, cy in enumerate(cyps):
        r = simulate(tr, cy, fm=0.5, t_end_min=15.0)
        grid_Ch[i, j] = cv_pv_ratio(r["Ch"])
        grid_Cm[i, j] = cv_pv_ratio(r["Cm"])

def plot_grid(ax, grid, title, cmap):
    im = ax.imshow(np.log10(grid), cmap=cmap, vmin=-2.5, vmax=2.5, aspect="equal")
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{grid[i,j]:.2f}", ha="center", va="center", fontsize=8.5,
                     color="white" if abs(np.log10(grid[i,j])) > 1.3 else "black")
    ax.set_xticks(range(3)); ax.set_xticklabels(tlabels, fontsize=7.5)
    ax.set_yticks(range(3)); ax.set_yticklabels(tlabels, fontsize=7.5)
    ax.set_xlabel("CYP zonation", fontsize=8)
    ax.set_ylabel("Transporter zonation", fontsize=8)
    ax.set_title(title, fontsize=9)
    return im

axB.set_title("B", loc="left", fontweight="bold", x=-0.24)
plot_grid(axB, grid_Ch, "Parent (hepatocyte Ch) CV:PV", "RdBu_r")

axC.set_title("C", loc="left", fontweight="bold", x=-0.24)
plot_grid(axC, grid_Cm, "Metabolite (hepatocyte Cm) CV:PV", "PRGn")

# ---------------- Panels D/E: DDI discrimination, split into two single-axis panels ----------------
dig_doses = [0, 0.3, 0.6, 1.0, 1.5, 2.4, 4.0]
inh_doses = [0, 0.6, 1.2, 2.0, 3.0, 5.0, 8.0]

ch_dig, cm_dig = [], []
for d in dig_doses:
    r = simulate("flat", "pericentral", fm=0.5, digoxin_uM=d, t_end_min=15.0)
    ch_dig.append(cv_pv_ratio(r["Ch"])); cm_dig.append(cv_pv_ratio(r["Cm"]))

ch_inh, cm_inh = [], []
for d in inh_doses:
    r = simulate("flat", "pericentral", fm=0.5, cyp_inh_uM=d, t_end_min=15.0)
    ch_inh.append(cv_pv_ratio(r["Ch"])); cm_inh.append(cv_pv_ratio(r["Cm"]))

ch_dig_pct = [100*c/ch_dig[0] for c in ch_dig]
cm_dig_pct = [100*c/cm_dig[0] for c in cm_dig]
ch_inh_pct = [100*c/ch_inh[0] for c in ch_inh]
cm_inh_pct = [100*c/cm_inh[0] for c in cm_inh]
y_all = ch_dig_pct + cm_dig_pct + ch_inh_pct + cm_inh_pct
ylim = (0, max(y_all)*1.15)

def plot_ddi_panel(ax, doses, ch_pct, cm_pct, title, panel_letter, xlabel):
    ax.plot(range(len(doses)), ch_pct, "o-", color="#2a78d6", label="Parent (Ch) CV:PV", ms=6, mfc="white", mew=1.6)
    ax.plot(range(len(doses)), cm_pct, "s--", color="#7a4fa8", label="Metabolite (Cm) CV:PV", ms=6, mfc="white", mew=1.6)
    ax.annotate(f"{ch_pct[-1]:.0f}%", (len(doses)-1, ch_pct[-1]), textcoords="offset points",
                xytext=(6, 4), fontsize=7, color="#2a78d6")
    ax.annotate(f"{cm_pct[-1]:.0f}%", (len(doses)-1, cm_pct[-1]), textcoords="offset points",
                xytext=(6, -12), fontsize=7, color="#7a4fa8")
    ax.set_xticks(range(len(doses)))
    ax.set_xticklabels([f"{d:g}" for d in doses], fontsize=7.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("% of no-perpetrator\nbaseline CV:PV")
    ax.set_ylim(ylim)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    ax.set_title(title, fontsize=9)
    ax.set_title(panel_letter, loc="left", fontweight="bold", x=-0.22)

plot_ddi_panel(axD, dig_doses, ch_dig_pct, cm_dig_pct,
                "Transporter inhibitor (digoxin-like)", "D", "Digoxin (µM)")
plot_ddi_panel(axE, inh_doses, ch_inh_pct, cm_inh_pct,
                "CYP inhibitor (ketoconazole-like)", "E", "CYP inhibitor (µM)")

fig.text(0.5, 0.005,
    "D, E share the same y-axis scale. Victim drug: transporter = flat, CYP = pericentral, fm = 0.5. "
    "Signature (E reshapes spatial contrast more than D) held in 7/9 tested fm x zonation configurations (Figure 2B).",
    ha="center", fontsize=7.5, color="#555", style="italic")

fig.savefig(__file__.replace("make_figure1.py", "figure1.png"), dpi=300, bbox_inches="tight")
print("saved")
print("grid_Ch=", grid_Ch.tolist())
print("grid_Cm=", grid_Cm.tolist())
