import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

# ==== USER SETTINGS ==========================================================
paths = [
    [
        "Pollinator_basin_A_xu_0.10_f_0.50_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.20_f_0.50_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.30_f_0.50_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.40_f_0.50_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.50_f_0.50_gamma0_1.8.txt",
    ],
    [
        "Pollinator_basin_A_xu_0.10_f_0.40_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.20_f_0.40_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.30_f_0.40_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.40_f_0.40_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.50_f_0.40_gamma0_1.8.txt",
    ],
    [
        "Pollinator_basin_A_xu_0.10_f_0.30_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.20_f_0.30_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.30_f_0.30_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.40_f_0.30_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.50_f_0.30_gamma0_1.8.txt",
    ],
    [
        "Pollinator_basin_A_xu_0.10_f_0.20_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.20_f_0.20_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.30_f_0.20_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.40_f_0.20_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.50_f_0.20_gamma0_1.8.txt",
    ],
    [
        "Pollinator_basin_A_xu_0.10_f_0.10_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.20_f_0.10_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.30_f_0.10_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.40_f_0.10_gamma0_1.8.txt",
        "Pollinator_basin_A_xu_0.50_f_0.10_gamma0_1.8.txt",
    ],
]

xu_vals = [0.1, 0.2, 0.3, 0.4, 0.5]
f_vals  = [0.5, 0.4, 0.3, 0.2, 0.1]

# Column indices
X_COL, Y_COL, Z_COL = 0, 1, 2
OUTPUT_START = 3
K_AMONG_OUTPUTS = 18
CLASS_COL = OUTPUT_START + (K_AMONG_OUTPUTS - 1)

MARKER_SIZE = 6
ALPHA = 0.01


# ============================================================================

def style_3d_axes(ax, fontsize=11, linewidth=2):
    """Uniform styling for all 3D subplots."""
    ax.grid(False)
#    ax.tick_params(axis="both", labelsize=fontsize, width=linewidth)
#    ax.tick_params(axis="z", labelsize=fontsize, width=linewidth)

    # Thicker axis lines
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.line.set_linewidth(linewidth)

    # Hide tick labels
    ax.set_xticklabels([])
    ax.set_yticklabels([]) 
    ax.set_zticklabels([])
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)

    ax.set_box_aspect([1, 1, 0.8])

def add_axis_arrows(ax):
    """
    Draw arrows along each axis from 0 -> 1 in data space.
    Assumes your data live in [0,1]^3; sets limits accordingly.
    """
    # Ensure limits match the requested arrow domain
#    ax.set_xlim(0, 1)
#    ax.set_ylim(0, 1)
#    ax.set_zlim(0, 1)

    # Draw 3 arrows from the origin along +x, +y, +z
    # arrow_length_ratio controls the head size
#    ax.quiver(0, 0, 0, 1, 0, 0, arrow_length_ratio=0.08, linewidths=2, color="k")
#    ax.quiver(0, 0, 0, 0, 1, 0, arrow_length_ratio=0.08, linewidths=2, color="k")
#    ax.quiver(0, 0, 0, 0, 0, 1, arrow_length_ratio=0.08, linewidths=2, color="k")

    # Optional: tiny labels at the arrow tips (comment out if you don't want them)
#    ax.text(1.02, 0.00, 0.00, "x", ha="left", va="center", fontsize=9)
#    ax.text(0.00, 1.02, 0.00, "y", ha="center", va="bottom", fontsize=9)
#    ax.text(0.00, 0.00, 1.02, "z", ha="center", va="bottom", fontsize=9)

def plot_scatter(ax, df):
    x = df.iloc[:, X_COL].to_numpy()
    y = df.iloc[:, Y_COL].to_numpy()
    z = df.iloc[:, Z_COL].to_numpy()
    c_vals = df.iloc[:, CLASS_COL].astype(float).to_numpy()
#    colors = np.where(c_vals < 0.01, "crimson", "mediumaquamarine")
    colors = np.where(c_vals < 0.01, "#C43C39", "#2A9D8F")
    ax.scatter(x, y, z, c=colors, s=MARKER_SIZE, alpha=ALPHA)


def load_txt(path):
    candidate = (BASE_DIR / path).resolve()
    if not candidate.exists():
        print(f"[warn] File not found: {candidate}")
        return None
    return pd.read_csv(candidate, sep=r"\s+", header=None, comment="#", engine="python")

def ratio_nonzero_over_zero(df):
    """Compute (#non-zero) / (#zero) for the class column used for coloring."""
    c = df.iloc[:, CLASS_COL].astype(float).to_numpy()
    zeros = np.sum(c < 0.01)
    nonzeros = np.sum(c >= 0.01)
    if zeros == 0 and nonzeros == 0:
        return None  # no data
    if zeros == 0:
        return np.inf
    return nonzeros / 27000  # (kept as in your script)

def main(save_path="pollinator_basin_5x5.png"):
    save_path = Path(save_path)
    fig = plt.figure(figsize=(12, 12))

    # big background for arrows + labels
    bg = fig.add_axes([0, 0, 1, 1], frameon=False, zorder=5)
    bg.set_xlim(0, 1); bg.set_ylim(0, 1); bg.axis("off")

    arrow_style = dict(arrowstyle='-|>', color='black', lw=4, mutation_scale=15)

    # outer arrows
    bg.annotate("", xy=(0.96, 0.08), xytext=(0.08, 0.08),
                arrowprops=arrow_style, xycoords="axes fraction")
    bg.annotate("", xy=(0.08, 0.92), xytext=(0.08, 0.08),
                arrowprops=arrow_style, xycoords="axes fraction")

    bg.text(0.538, -0.01, r"$\mathbf{x_u}$", fontsize=34)#, ha="center", va="center")
    bg.text(-0.015, 0.51, r"$\mathbf{f}$",   fontsize=34)#, va="center", ha="center")

    bg.text(0.01, 0.96, "(a)", fontsize=30, fontweight="bold", ha="left", va="top")

    nrows, ncols = 5, 5
    for i in range(nrows):
        for j in range(ncols):
            idx = i * ncols + j + 1
            ax = fig.add_subplot(nrows, ncols, idx, projection="3d")

            df = load_txt(paths[i][j])
            if df is not None:
                plot_scatter(ax, df)
                r = ratio_nonzero_over_zero(df)
                if r is None:
                    label = "N/A"
                elif np.isinf(r):
                    label = "∞"
                else:
                    label = f"{r:.3f}"
            else:
                ax.text2D(0.4, 0.5, "Missing file",
                          transform=ax.transAxes, color="red", fontsize=8)
                label = "N/A"

            style_3d_axes(ax, fontsize=10, linewidth=3)

            # >>> add per-subplot axis arrows 0->1 on x,y,z <<<
            add_axis_arrows(ax)

            # --- write bold value above each cube ---
            ax.text2D(
                0.5, 0.93, label,
                transform=ax.transAxes,
                ha="center", va="bottom",
                fontsize=16, fontweight="bold",
                bbox=dict(facecolor="white", alpha=0.2, edgecolor="none", pad=1.5),
                zorder=10
            )

    # numeric labels along the outer axes
    for j, xu in enumerate(xu_vals):
        x_pos = 0.20 + j * 0.18
        fig.text(x_pos, 0.03, f"{xu:.1f}", ha="center", fontsize=24, fontweight="bold")
    for i, f in enumerate(f_vals):
        y_pos = 0.86 - i * 0.17
        fig.text(0.009, y_pos, f"{f:.1f}", va="center", fontsize=24, fontweight="bold")

    plt.subplots_adjust(left=0.10, right=0.98, bottom=0.10, top=0.95, wspace=0.05, hspace=0.05)
    plt.savefig(save_path, dpi=600, bbox_inches="tight", pad_inches=0.05)
    print(f"[info] saved: {save_path}")
    plt.show()

if __name__ == "__main__":    
    main(BASE_DIR / "pollinator_basin_5x5.png")

