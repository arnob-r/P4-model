#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib as mpl

BASE_DIR = Path(__file__).resolve().parent


def loadtxt_first(*relative_paths):
    for rel in relative_paths:
        candidate = (BASE_DIR / rel).resolve()
        if candidate.exists():
            return np.loadtxt(candidate)
    joined = ", ".join(relative_paths)
    raise FileNotFoundError(f"None of these files were found: {joined}")


# =============================================================================
# Load datasets
# =============================================================================
data1 = loadtxt_first("Pollinator_transition_g_xi_f0.1.txt",
                      "../Pollinator_transition_g_xi_f_0.1.txt")
data2 = loadtxt_first("Pollinator_transition_g_xi_f0.2.txt",
                      "../Pollinator_transition_g_xi_f_0.2.txt")
data3 = loadtxt_first("Pollinator_transition_g_xi_f0.3.txt",
                      "../Pollinator_transition_g_xi_f_0.3.txt")
data4 = loadtxt_first("Pollinator_transition_g_xi_xu0.1.txt")
data5 = loadtxt_first("Pollinator_transition_g_xi_xu0.2.txt")
data6 = loadtxt_first("Pollinator_transition_g_xi_xu0.3.txt")

datasets = [
    (data1, r"$\mathbf{f=0.1}$"),
    (data2, r"$\mathbf{f=0.2}$"),
    (data3, r"$\mathbf{f=0.3}$"),
    (data4, r"$\mathbf{x_u=0.1}$"),
    (data5, r"$\mathbf{x_u=0.2}$"),
    (data6, r"$\mathbf{x_u=0.3}$"),
]

panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]


# =============================================================================
# Colormap
# =============================================================================
cmap = mpl.colormaps["cividis"].copy()
cmap.set_bad(color="#FFFFFF")


# =============================================================================
# Global plot settings
# =============================================================================
plt.rcParams.update({
    "font.size": 28,
    "font.weight": "bold",
    "axes.labelsize": 28,
    "axes.labelweight": "bold",
    "axes.titlesize": 26,
    "axes.titleweight": "bold",
    "axes.linewidth": 3,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "xtick.major.width": 3,
    "ytick.major.width": 3,
    "xtick.minor.width": 3,
    "ytick.minor.width": 3,
    "hatch.linewidth": 2,
})


# =============================================================================
# Create figure
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(24, 16), constrained_layout=True)
axes = axes.flatten()


# =============================================================================
# Global color scale (using only values > 0.09)
# =============================================================================
z_all = np.concatenate([
    d[d[:, 2] > 0.09, 2] for d, _ in datasets if np.any(d[:, 2] > 0.09)
])

if z_all.size == 0:
    raise RuntimeError("No z-values > 0.09 found across datasets; cannot set vmin/vmax.")

vmin = float(z_all.min())
vmax = float(z_all.max())


# =============================================================================
# Plot each dataset
# =============================================================================
ims = []

for idx, (ax, (data, label)) in enumerate(zip(axes, datasets)):
    x, y, z = data[:, 0], data[:, 1], data[:, 2]

    mask = z > 0.09
    if not np.any(mask):
        ax.set_title(label + "\n(No points > 0.09)", fontsize=18, fontweight='bold')
        ax.axis('off')
        ims.append(None)
        continue

    # -------------------------------------------------------------------------
    # Build rectangular grid
    # -------------------------------------------------------------------------
    xi = np.unique(x)
    yi = np.unique(y)
    X, Y = np.meshgrid(xi, yi)
    Z = np.full((len(yi), len(xi)), np.nan, dtype=float)

    for i, xv in enumerate(xi):
        for j, yv in enumerate(yi):
            m = (x == xv) & (y == yv)
            if np.any(m):
                val = z[m][0]
                if val > 0.09:
                    Z[j, i] = val
                else:
                    Z[j, i] = np.nan

    # -------------------------------------------------------------------------
    # Main color plot
    # -------------------------------------------------------------------------
    im = ax.pcolormesh(
        X, Y, Z,
        cmap=cmap,
        shading='auto',
        vmin=vmin,
        vmax=vmax,
        linewidth=0.7,
        edgecolors='face'
    )
    ims.append(im)

    # -------------------------------------------------------------------------
    # Identify right-side NaN region only
    # -------------------------------------------------------------------------
    nan_mask = np.isnan(Z)
    right_nan_mask = np.zeros_like(Z, dtype=float)

    for j in range(Z.shape[0]):
        row = nan_mask[j, :]
        non_nan_indices = np.where(~row)[0]

        if len(non_nan_indices) > 0:
            last_valid = non_nan_indices[-1]
            right_nan_mask[j, last_valid + 1:] = 1.0
        else:
            right_nan_mask[j, :] = 1.0

    # -------------------------------------------------------------------------
    # Apply sparse hatch only on right NaN region
    # -------------------------------------------------------------------------
    if np.any(right_nan_mask > 0):
        hatch_overlay = ax.contourf(
            X, Y, right_nan_mask,
            levels=[0.5, 1.5],
            hatches=['/'],
            colors='none'
        )

        if hasattr(hatch_overlay, "collections"):
            for c in hatch_overlay.collections:
                c.set_edgecolor('black')
                c.set_linewidth(0.4)
        else:
            hatch_overlay.set_edgecolor('black')
            hatch_overlay.set_linewidth(0.4)

    # -------------------------------------------------------------------------
    # Axis formatting
    # -------------------------------------------------------------------------
    ax.set_xlabel(r"$\mathbf{\gamma_0}$", fontsize=28, fontweight='bold')
    ax.set_ylabel(r"$\mathbf{\xi_0}$", fontsize=28, fontweight='bold')
    ax.set_title(label, fontsize=28, fontweight='bold', pad=18)

    ax.tick_params(width=3, length=10, direction='inout')
    ax.set_xlim(0.5, 2.1)
    ax.set_ylim(0.1, 3.0)

    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight('bold')

    # -------------------------------------------------------------------------
    # Panel label: (a), (b), ..., (f)
    # -------------------------------------------------------------------------
    ax.text(
        -0.12, 1.05, panel_labels[idx],
        transform=ax.transAxes,
        fontsize=28,
        fontweight='bold',
        va='bottom',
        ha='left'
    )


# =============================================================================
# Colorbars
# =============================================================================
cbar1 = fig.colorbar(ims[2], ax=axes[2], shrink=0.85, pad=0.02)
cbar1.set_label(r"$\mathbf{x_u^*}$", fontsize=28, fontweight="bold")
cbar1.ax.tick_params(labelsize=24, width=3)
cbar1.outline.set_linewidth(3)
for t in cbar1.ax.get_yticklabels():
    t.set_fontweight("bold")
cbar1.solids.set_edgecolor("face")

cbar2 = fig.colorbar(ims[5], ax=axes[5], shrink=0.85, pad=0.02)
cbar2.set_label(r"$\mathbf{f^*}$", fontsize=28, fontweight="bold")
cbar2.ax.tick_params(labelsize=24, width=3)
cbar2.outline.set_linewidth(3)
for t in cbar2.ax.get_yticklabels():
    t.set_fontweight("bold")
cbar2.solids.set_edgecolor("face")


# =============================================================================
# Save and show
# =============================================================================
fig.savefig(
    BASE_DIR / "Pollinator_Transition_Maps_with_panel_labels.png",
    dpi=600,
    bbox_inches="tight"
)
plt.show()