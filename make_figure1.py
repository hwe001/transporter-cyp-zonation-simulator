import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from model import simulate, cv_pv_ratio, N, xs

plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})

fig = plt.figure(figsize=(13.6, 8.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], wspace=0.32, hspace=0.42)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

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

# ---------------- Panel D: DDI discrimination, flat transporter + pericentral CYP ----------------
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

axD.set_title("D", loc="left", fontweight="bold", x=-0.16)
ax2 = axD.twinx()
l1, = axD.plot(range(len(dig_doses)), [c/ch_dig[0] for c in ch_dig], "o-", color="#2a78d6",
                label="Parent CV:PV (transporter-inhib.)", ms=6, mfc="white", mew=1.6)
l2, = axD.plot(range(len(dig_doses)), [c/ch_inh[0] for c in ch_inh], "s--", color="#eb6834",
                label="Parent CV:PV (CYP-inhib.)", ms=6, mfc="white", mew=1.6)
l3, = ax2.plot(range(len(dig_doses)), [c/cm_dig[0] for c in cm_dig], "^-", color="#2a78d6",
                alpha=0.5, label="Metabolite CV:PV (transporter-inhib.)", ms=6)
l4, = ax2.plot(range(len(dig_doses)), [c/cm_inh[0] for c in cm_inh], "v--", color="#eb6834",
                alpha=0.5, label="Metabolite CV:PV (CYP-inhib.)", ms=6)
axD.set_xticks(range(len(dig_doses)))
axD.set_xticklabels([f"dose {k}" for k in range(len(dig_doses))], fontsize=7)
axD.set_xlabel("Perpetrator dose level (arbitrary steps)")
axD.set_ylabel("Parent CV:PV\n(fraction of no-perpetrator)")
ax2.set_ylabel("Metabolite CV:PV\n(fraction of no-perpetrator)")
axD.spines[['top']].set_visible(False)
lines = [l1, l2, l3, l4]
axD.legend(lines, [l.get_label() for l in lines], fontsize=6.3, frameon=False, loc="upper left")
axD.text(0.98, 0.03, "transporter=flat, CYP=pericentral", transform=axD.transAxes,
          ha="right", fontsize=6.5, color="#666", style="italic")

fig.savefig(__file__.replace("make_figure1.py", "figure1.png"), dpi=300, bbox_inches="tight")
print("saved")
print("grid_Ch=", grid_Ch.tolist())
print("grid_Cm=", grid_Cm.tolist())
