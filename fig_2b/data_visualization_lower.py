import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent


def loadtxt_first(*relative_paths):
    for rel in relative_paths:
        candidate = (BASE_DIR / rel).resolve()
        if candidate.exists():
            return np.loadtxt(candidate)
    joined = ", ".join(relative_paths)
    raise FileNotFoundError(f"None of these files were found: {joined}")

# === Load data ===
data1 = loadtxt_first("Pollinator_bif_fwd_xu0.1_0.2.txt")
data2 = loadtxt_first("Pollinator_bif_fwd_xu0.1_0.5.txt")
data3 = loadtxt_first("Pollinator_bif_fwd_xu0.1_0.8.txt")

datasets = [data1, data2, data3]
titles = [r"$\mathbf{IC_1}$", r"$\mathbf{IC_2}$", r"$\mathbf{IC_3}$"]

# === Columns & degrees ===
Y_cols = [3, 5, 7, 11, 12, 17, 18, 19, 21]
Y_deg  = [1, 3, 2, 1, 1, 3, 7, 1, 2]

# === Custom color map by degree ===
color_map = {1: "y", 2: "teal", 3: "tomato", 7: "navy"}

# === Create 3 side-by-side subplots ===
fig, axes = plt.subplots(1, 3, figsize=(24, 6.5), sharey=True)

for ax, data, title in zip(axes, datasets, titles):
    X = data[:, 0]

    # Plot each column with color based on degree
    for col, deg in zip(Y_cols, Y_deg):
        color = color_map.get(deg, "gray")
        ax.plot(X, data[:, col], 'o-', color=color,
                label=f"deg {deg}", linewidth=3.5, markersize=6)

    # Axis labels & titles
    ax.set_xlabel(r"$\mathbf{f}$", fontsize=32, fontweight='bold')
    ax.set_title(title, fontsize=32, fontweight='bold', pad =15)
    ax.grid(False)
    ax.tick_params(axis='both', labelsize=28, width=3, length=7)

    # Bold tick labels
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontweight('bold')

    # Thicken frame (spines)
    for spine in ax.spines.values():
        spine.set_linewidth(4)

    ax.set_ylim(-0.03, 0.6)
# Shared y-label
axes[0].set_ylabel(r"$\boldsymbol{A}$", fontsize=32, fontweight='bold')

# === Unified legend (only deg 1, 2, 3, 7) ===
# === Unified legend (only deg 1, 2, 3, 7) ===
legend_degrees = [1, 2, 3, 7]
handles = [
    plt.Line2D([], [], color=color_map[deg], marker='o', markersize=10,
               linewidth=3.5, label=f"Degree {deg}")
    for deg in legend_degrees
]
fig.legend(handles=handles,
           loc='lower center',
           ncol=4,
           fontsize=32,
           frameon=False,
           handletextpad=0.8,
           columnspacing=1.2,
           )

# === Layout adjustment (gap between plots and legend) ===
#plt.subplots_adjust(bottom=0.35, wspace=0.25)

# === Layout adjustment ===
plt.subplots_adjust(bottom=0.35, wspace=0.1, top=0.88)

# === Save high-resolution image ===
plt.savefig(BASE_DIR / "Pollinator_bifurcation_panels_xu_fixed.png", dpi=600, bbox_inches='tight')

plt.show()
