import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent


def loadtxt(name):
    return np.loadtxt((BASE_DIR / name).resolve())

# =========================
# Load data
# =========================
data1 = loadtxt("Pollinator_basin_counts_f0.1.txt")
data2 = loadtxt("Pollinator_basin_counts_f0.2.txt")
data3 = loadtxt("Pollinator_basin_counts_f0.3.txt")

data4 = loadtxt("Pollinator_basin_counts_xu0.1.txt")
data5 = loadtxt("Pollinator_basin_counts_xu0.2.txt")
data6 = loadtxt("Pollinator_basin_counts_xu0.3.txt")

# Panel 1: BS vs x_u for f=0.1,0.2,0.3
datasets_upper = [
    (data1, r"$f=0.1$", 0.1),
    (data2, r"$f=0.2$", 0.2),
    (data3, r"$f=0.3$", 0.3),
]

# Panel 2: BS vs f for x_u=0.1,0.2,0.3
datasets_lower = [
    (data4, r"$x_u=0.1$", 0.1),
    (data5, r"$x_u=0.2$", 0.2),
    (data6, r"$x_u=0.3$", 0.3),
]

# =========================
# Columns (adjust ONLY if your files differ)
# =========================
XCOL = 0      # x-axis column
YCOL = 3      # "BS" column

# Marker map (your request)
MARKER = {0.1: "o", 0.2: "s", 0.3: "^"}

# =========================
# Styling helpers
# =========================
def style_axis(ax, xlabel):
    ax.set_xlabel(xlabel, fontsize=18, fontweight="bold")
    ax.set_ylabel("BS", fontsize=18, fontweight="bold")
    ax.set_ylim(-0.03, 1.03)
    ax.grid(False)
    ax.tick_params(axis="both", labelsize=13, width=2, right=False)
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_fontweight("bold")
    for s in ax.spines.values():
        s.set_linewidth(2.2)

def plot_panel(ax, datasets, xlabel_math):
    for data, lbl, key in datasets:
        x = data[:, XCOL]
        bs = data[:, YCOL]
        bs_non = 1.0 - bs

        m = MARKER.get(key, "o")

        # connecting lines
        ax.plot(x, bs, "-", color="#C43C39", lw=2.2, alpha=0.5)
        ax.plot(x, bs_non, "--", color="#2A9D8F", lw=2.2, alpha=0.5)
        
        
        # markers
        ax.scatter(x, bs, s=34, marker=m, c="#C43C39",
                   edgecolors="k", linewidths=0.05, zorder=3, alpha=0.5)
        ax.scatter(x, bs_non, s=40, marker=m, c="#2A9D8F",
                   edgecolors="k", linewidths=0.05, zorder=3, alpha=0.5)

        # legend handle (black)
        ax.plot([], [], color="w", marker=m, linestyle="-", lw=0.01,
                markersize=5, label=lbl)

    style_axis(ax, xlabel=xlabel_math)
#    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22),
#              ncol=3, frameon=False, fontsize=11,
#              handlelength=1.6, handletextpad=0.5)

# =========================
# Make INDIVIDUAL figures
# =========================

# ---- Figure 1 (upper): BS vs x_u ----
fig1, ax1 = plt.subplots(figsize=(6.0, 4))
plot_panel(ax1, datasets_upper, xlabel_math=r"$\mathbf{x_u}$")
ax1.text(-0.04, 1.2, "(b)", fontsize=18, fontweight="bold", ha="left", va="top")
plt.tight_layout()
plt.savefig(BASE_DIR / "Pollinator_basin_upper.png", dpi=600, bbox_inches="tight")
plt.show()

# ---- Figure 2 (lower): BS vs f ----
fig2, ax2 = plt.subplots(figsize=(6.0, 4))
plot_panel(ax2, datasets_lower, xlabel_math=r"$\mathbf{f}$")
ax2.text(-0.04, 1.2, "(c)", fontsize=18, fontweight="bold", ha="left", va="top")
plt.tight_layout()
plt.savefig(BASE_DIR / "Pollinator_basin_lower.png", dpi=600, bbox_inches="tight")
plt.show()
