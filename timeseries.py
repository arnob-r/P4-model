#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Time-series simulation for the P4 model, Eq. 1.

This script solves the pollinator-plant-pest-pesticide model using RK4
and plots the time evolution of pollinator A18 for three management cases:

(a) x_u = 0.1, f = 0.1
(b) x_u = 0.1, f = 0.5
(c) x_u = 0.5, f = 0.1

Initial condition:
IC1: A_i(0) = P_i(0) = S_i(0) = 0.2
"""

import os
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. MODEL CONFIGURATION
# ============================================================

N_A = 23
N_P = 36
N_S = 51

alpha = 0.1
kappa = 0.0001

gamma0 = 1.8
xi0 = 1.0

mA = 0.5
mS = 0.9

h_holling = 0.7
mu = 0.0

theta = 0.5

uH = 1.0
uL = 0.5

biiA, bijA = 1.0, 0.1
biiP, bijP = 1.0, 0.1
biiS, bijS = 1.0, 0.1

dt = 0.01
tmax = 1000.0
save_every = 100

IC1_VALUE = 0.2

# A18 means index 18 in one-based indexing.
# Python index is therefore 17.
A18_INDEX = 17


# ============================================================
# 2. LOAD ADJACENCY MATRICES
# ============================================================

def load_adjacencies(
    pollinator_plant_file="PPo_adj_A.txt",
    plant_pest_file="PPe_adj_A.txt"
):
    if not os.path.exists(pollinator_plant_file):
        raise FileNotFoundError(
            f"Cannot find {pollinator_plant_file}. "
            "Place this file in the same folder as this script."
        )

    if not os.path.exists(plant_pest_file):
        raise FileNotFoundError(
            f"Cannot find {plant_pest_file}. "
            "Place this file in the same folder as this script."
        )

    Adj1 = np.loadtxt(pollinator_plant_file)
    Adj2 = np.loadtxt(plant_pest_file)

    if Adj1.shape != (N_A, N_P):
        raise ValueError(
            f"PPo_adj_A.txt must have shape {(N_A, N_P)}, "
            f"but got {Adj1.shape}."
        )

    if Adj2.shape != (N_P, N_S):
        raise ValueError(
            f"PPe_adj_A.txt must have shape {(N_P, N_S)}, "
            f"but got {Adj2.shape}."
        )

    return Adj1.astype(float), Adj2.astype(float)


# ============================================================
# 3. COMPETITION MATRICES
# ============================================================

def competition_matrix(n, bii=1.0, bij=0.1):
    B = np.full((n, n), bij, dtype=float)
    np.fill_diagonal(B, bii)
    return B


# ============================================================
# 4. DEGREE-SCALED INTERACTION MATRICES
# ============================================================

def safe_degree_power(degree, exponent):
    degree = np.asarray(degree, dtype=float)
    out = np.zeros_like(degree, dtype=float)

    mask = degree > 0
    out[mask] = degree[mask] ** exponent
    out[~mask] = np.inf

    return out


def build_couplings(Adj1, Adj2):
    D_A = Adj1.sum(axis=1)
    D_P_from_A = Adj1.sum(axis=0)
    D_P_to_S = Adj2.sum(axis=1)
    D_S = Adj2.sum(axis=0)

    denom_A = safe_degree_power(D_A, 1.0 - theta)
    denom_P_from_A = safe_degree_power(D_P_from_A, 1.0 - theta)
    denom_P_to_S = safe_degree_power(D_P_to_S, 1.0 - theta)
    denom_S = safe_degree_power(D_S, 1.0 - theta)

    gamma_AP = gamma0 * Adj1 / denom_A[:, None]
    gamma_PA = gamma0 * Adj1.T / denom_P_from_A[:, None]

    xi_PS = xi0 * Adj2 / denom_P_to_S[:, None]
    xi_SP = xi0 * Adj2.T / denom_S[:, None]

    return gamma_AP, gamma_PA, xi_PS, xi_SP


# ============================================================
# 5. EFFECTIVE PESTICIDE PRESSURE
# ============================================================

def effective_pesticide_pressure(xu, f):
    return (1.0 - f) * ((1.0 - xu) * uH + xu * uL)


# ============================================================
# 6. RIGHT-HAND SIDE OF EQ. 1
# ============================================================

def rhs(y, params):
    betaA = params["betaA"]
    betaP = params["betaP"]
    betaS = params["betaS"]

    gamma_AP = params["gamma_AP"]
    gamma_PA = params["gamma_PA"]

    xi_PS = params["xi_PS"]
    xi_SP = params["xi_SP"]

    xu = params["xu"]
    f = params["f"]

    A = y[:N_A]
    P = y[N_A:N_A + N_P]
    S = y[N_A + N_P:]

    ueff = effective_pesticide_pressure(xu, f)

    # Pollinator equation
    comp_A = betaA @ A
    mutual_A_raw = gamma_AP @ P
    mutual_A = mutual_A_raw / (1.0 + h_holling * mutual_A_raw)

    dA = A * (
        alpha
        - kappa
        - comp_A
        + mutual_A
        - mA * ueff
    ) + mu

    # Plant equation
    comp_P = betaP @ P
    mutual_P_raw = gamma_PA @ A
    mutual_P = mutual_P_raw / (1.0 + h_holling * mutual_P_raw)

    pest_loss_raw = xi_PS @ S
    pest_loss = pest_loss_raw / (1.0 + h_holling * pest_loss_raw)

    dP = P * (
        alpha
        - comp_P
        + mutual_P
        - pest_loss
    ) + mu

    # Pest equation
    comp_S = betaS @ S
    plant_gain_raw = xi_SP @ P
    plant_gain = plant_gain_raw / (1.0 + h_holling * plant_gain_raw)

    dS = S * (
        alpha
        - kappa
        - comp_S
        + plant_gain
        - mS * ueff
    ) + mu

    return np.concatenate([dA, dP, dS])


# ============================================================
# 7. RK4 INTEGRATOR
# ============================================================

def rk4_integrate(params, y0, dt=0.01, tmax=1000.0, save_every=100):
    n_steps = int(tmax / dt)
    n_save = n_steps // save_every + 1

    T = np.zeros(n_save)
    Y = np.zeros((n_save, y0.size))

    y = y0.copy()

    save_id = 0
    T[save_id] = 0.0
    Y[save_id] = y

    for step in range(1, n_steps + 1):

        k1 = rhs(y, params)
        k2 = rhs(y + 0.5 * dt * k1, params)
        k3 = rhs(y + 0.5 * dt * k2, params)
        k4 = rhs(y + dt * k3, params)

        y = y + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)

        # Avoid small negative numerical undershoots
        y = np.maximum(y, 0.0)

        if step % save_every == 0:
            save_id += 1
            T[save_id] = step * dt
            Y[save_id] = y

    return T, Y


# ============================================================
# 8. RUN ONE SCENARIO
# ============================================================

def run_scenario(xu, f, common_params):
    A0 = np.full(N_A, IC1_VALUE)
    P0 = np.full(N_P, IC1_VALUE)
    S0 = np.full(N_S, IC1_VALUE)

    y0 = np.concatenate([A0, P0, S0])

    params = dict(common_params)
    params["xu"] = xu
    params["f"] = f

    T, Y = rk4_integrate(
        params=params,
        y0=y0,
        dt=dt,
        tmax=tmax,
        save_every=save_every
    )

    A = Y[:, :N_A]
    P = Y[:, N_A:N_A + N_P]
    S = Y[:, N_A + N_P:]

    result = {
        "T": T,
        "A": A,
        "P": P,
        "S": S,
        "A18": A[:, A18_INDEX],
        "mean_A": A.mean(axis=1),
        "mean_P": P.mean(axis=1),
        "mean_S": S.mean(axis=1),
        "xu": xu,
        "f": f,
    }

    return result


# ============================================================
# 9. PLOTTING FUNCTION: THREE SUBFIGURES
# ============================================================

def plot_three_subfigures_A18(
    results,
    output_name="P4_A18_three_subfigures_IC1.png"
):
    fig, axes = plt.subplots(
        3, 1,
        figsize=(10, 16),
        sharex=True
    )

    panel_labels = ["(a)", "(b)", "(c)"]

    title_fontsize = 34
    panel_title_fontsize = 30
    axis_label_fontsize = 30
    tick_fontsize = 30

    timeseries_lw = 5.0
    box_lw = 3.5
    tick_lw = 3.0
    tick_length = 9

    for i, (ax, res) in enumerate(zip(axes, results)):

        ax.plot(
            res["T"],
            res["A18"],
            lw=timeseries_lw
        )

        # Fixed y-axis limit
        ax.set_ylim(-0.05, 0.55)
        ax.set_yticks(np.arange(0.0, 0.51, 0.1))

        ax.set_title(
            panel_labels[i] + " " +
            rf"$x_u={res['xu']:.1f},\ f={res['f']:.1f}$",
            fontsize=panel_title_fontsize,
            fontweight="bold",
            loc="left",
            pad=12
        )

        ax.set_ylabel(
            r"$\mathbf{A_{18}(t)}$",
            fontsize=axis_label_fontsize,
            fontweight="bold"
        )

        ax.tick_params(
            axis="both",
            labelsize=tick_fontsize,
            width=tick_lw,
            length=tick_length,
            direction="out"
        )

        for tick in ax.get_xticklabels():
            tick.set_fontweight("bold")

        for tick in ax.get_yticklabels():
            tick.set_fontweight("bold")

        for spine in ax.spines.values():
            spine.set_linewidth(box_lw)

        ax.grid(False)

    axes[-1].set_xlabel(
        "Time",
        fontsize=axis_label_fontsize,
        fontweight="bold",
        labelpad=10
    )

    fig.suptitle(
        r"Persistence and extinction dynamics of pollinator $\mathbf{A_{18}}$"
        "\n"
        r"under management parameters",
        fontsize=title_fontsize,
        fontweight="bold",
        y=0.985
    )

    plt.tight_layout(rect=[0, 0, 1, 0.945])
    plt.savefig(output_name, dpi=600, bbox_inches="tight")
    plt.show()


# ============================================================
# 10. SAVE OUTPUT AS CSV
# ============================================================

def save_results_to_csv(results):
    for i, res in enumerate(results, start=1):

        out = np.column_stack([
            res["T"],
            res["A18"],
            res["mean_A"],
            res["mean_P"],
            res["mean_S"]
        ])

        fname = (
            f"P4_IC1_case_{i}_"
            f"xu_{res['xu']:.1f}_"
            f"f_{res['f']:.1f}.csv"
        )

        np.savetxt(
            fname,
            out,
            delimiter=",",
            header="time,A18,mean_A,mean_P,mean_S",
            comments=""
        )

        print(f"Saved: {fname}")


# ============================================================
# 11. MAIN SCRIPT
# ============================================================

Adj1, Adj2 = load_adjacencies(
    pollinator_plant_file="PPo_adj_A.txt",
    plant_pest_file="PPe_adj_A.txt"
)

betaA = competition_matrix(N_A, bii=biiA, bij=bijA)
betaP = competition_matrix(N_P, bii=biiP, bij=bijP)
betaS = competition_matrix(N_S, bii=biiS, bij=bijS)

gamma_AP, gamma_PA, xi_PS, xi_SP = build_couplings(Adj1, Adj2)

common_params = {
    "betaA": betaA,
    "betaP": betaP,
    "betaS": betaS,
    "gamma_AP": gamma_AP,
    "gamma_PA": gamma_PA,
    "xi_PS": xi_PS,
    "xi_SP": xi_SP,
}

scenarios = [
    {"xu": 0.1, "f": 0.1},
    {"xu": 0.1, "f": 0.5},
    {"xu": 0.5, "f": 0.1},
]

results = []

for sc in scenarios:
    print(f"Running xu={sc['xu']}, f={sc['f']} ...")

    res = run_scenario(
        xu=sc["xu"],
        f=sc["f"],
        common_params=common_params
    )

    results.append(res)

save_results_to_csv(results)

plot_three_subfigures_A18(
    results,
    output_name="P4_A18_three_subfigures_IC1.png"
)

print("Done.")