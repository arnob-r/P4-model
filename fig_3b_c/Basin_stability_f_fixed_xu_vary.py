#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code is for calculating basin stability f fixed amd xu is varied. Parallel JIT-accelerated sweep over multiple xu values.
Output: xu, count(A18<0.01), total, fraction
"""

import numpy as np
from pathlib import Path
from itertools import product
from time import perf_counter
import os, math
from multiprocessing import get_context
from numba import njit
from tqdm import tqdm

# -------------------- helpers --------------------
def load_adjacency(fname: str, shape: tuple, allow_transpose: bool = True) -> np.ndarray:
    p = Path(fname)
    if not p.exists():
        p = Path("/mnt/data") / Path(fname).name
    M = np.loadtxt(p, dtype=np.float64)
    if M.ndim == 1:
        M = M[None, :]
    if M.shape == shape:
        return M
    if allow_transpose and M.T.shape == shape:
        print(f"[note] {fname} was {M.shape}; using its transpose {M.T.shape}")
        return M.T
    raise ValueError(f"{fname}: expected {shape} or transpose, got {M.shape}")

def build_beta(n: int, bii=1.0, bij=0.1) -> np.ndarray:
    B = np.full((n, n), bij, dtype=np.float64)
    np.fill_diagonal(B, bii)
    return B

# -------------------- JIT core --------------------
@njit(fastmath=True, cache=True)
def rhs_jit(A, P, S, params, beta1, beta2, beta3, Adj1, Adj2, D1A, D1P, D2P, D2S):
    alpha, kappa, gamma0, xi0, mA, mS, f, ubar, mu, tt, hh = params
    sum1 = beta1 @ A
    sum2 = Adj1 @ P
    sum3 = beta2 @ P
    sum4 = Adj1.T @ A
    sum5 = Adj2 @ S
    sum6 = beta3 @ S
    sum7 = Adj2.T @ P
    powA  = np.power(D1A, (1.0 - tt))
    powP1 = np.power(D1P, (1.0 - tt))
    powP2 = np.power(D2P, (1.0 - tt))
    powS  = np.power(D2S, (1.0 - tt))
    gainA = np.zeros_like(A)
    for i in range(A.size):
        if D1A[i] > 0:
            gainA[i] = (gamma0 / powA[i]) * sum2[i]
    gainP1 = np.zeros_like(P)
    for j in range(P.size):
        if D1P[j] > 0:
            gainP1[j] = (gamma0 / powP1[j]) * sum4[j]
    gainP2 = np.zeros_like(P)
    for j in range(P.size):
        if D2P[j] > 0:
            gainP2[j] = (xi0 / powP2[j]) * sum5[j]
    gainS = np.zeros_like(S)
    for k in range(S.size):
        if D2S[k] > 0:
            gainS[k] = (xi0 / powS[k]) * sum7[k]
    def sat(g): return g / (1.0 + hh * g)
    dA = A * (alpha - kappa - sum1 + sat(gainA) - mA * (1.0 - f) * ubar) + mu
    dP = P * (alpha - sum3 + sat(gainP1) - sat(gainP2)) + mu
    dS = S * (alpha - kappa - sum6 + sat(gainS) - mS * (1.0 - f) * ubar) + mu
    return dA, dP, dS

@njit(fastmath=True, cache=True)
def integrate_once_jit(A0, P0, S0, params,
                       beta1, beta2, beta3, Adj1, Adj2, D1A, D1P, D2P, D2S,
                       h=0.01, mm=200000):
    A = A0.copy(); P = P0.copy(); S = S0.copy()
    for _ in range(mm):
        k1A, k1P, k1S = rhs_jit(A, P, S, params, beta1, beta2, beta3, Adj1, Adj2, D1A, D1P, D2P, D2S)
        A2 = A + 0.5 * h * k1A; P2 = P + 0.5 * h * k1P; S2 = S + 0.5 * h * k1S
        k2A, k2P, k2S = rhs_jit(A2, P2, S2, params, beta1, beta2, beta3, Adj1, Adj2, D1A, D1P, D2P, D2S)
        A3 = A + 0.5 * h * k2A; P3 = P + 0.5 * h * k2P; S3 = S + 0.5 * h * k2S
        k3A, k3P, k3S = rhs_jit(A3, P3, S3, params, beta1, beta2, beta3, Adj1, Adj2, D1A, D1P, D2P, D2S)
        A4 = A + h * k3A; P4 = P + h * k3P; S4 = S + h * k3S
        k4A, k4P, k4S = rhs_jit(A4, P4, S4, params, beta1, beta2, beta3, Adj1, Adj2, D1A, D1P, D2P, D2S)
        A += (h/6.0)*(k1A + 2*(k2A+k3A) + k4A)
        P += (h/6.0)*(k1P + 2*(k2P+k3P) + k4P)
        S += (h/6.0)*(k1S + 2*(k2S+k3S) + k4S)
    return A, P, S

# -------------------- Globals for workers --------------------
_G_PARAMS = None
_G_CACHES = None
_G_H = None
_G_MM = None
_G_N1 = None

def _init_worker(params, caches, h, mm, n1):
    global _G_PARAMS, _G_CACHES, _G_H, _G_MM, _G_N1
    _G_PARAMS = params; _G_CACHES = caches
    _G_H = h; _G_MM = mm; _G_N1 = n1

def _run_combo(combo):
    Ai0, Pj0, Sk0 = combo
    n1 = _G_N1; n2 = _G_CACHES["Adj1"].shape[1]; n3 = _G_CACHES["Adj2"].shape[1]
    A0 = np.full(n1, Ai0); P0 = np.full(n2, Pj0); S0 = np.full(n3, Sk0)
    A_fin, _, _ = integrate_once_jit(A0, P0, S0, _G_PARAMS,
                                     _G_CACHES["beta1"], _G_CACHES["beta2"], _G_CACHES["beta3"],
                                     _G_CACHES["Adj1"], _G_CACHES["Adj2"],
                                     _G_CACHES["D1A"], _G_CACHES["D1P"], _G_CACHES["D2P"], _G_CACHES["D2S"],
                                     h=_G_H, mm=_G_MM)
    return A_fin[17] < 0.01   # A(18)

# -------------------- Sweep for one xu --------------------
def sweep_one_xu(xu, out_path, levels=10, mm=200000, h=0.01, n_workers=None, show_progress=True):
    n1, n2, n3 = 23, 36, 51
    tt=0.5; alpha=0.1; kappa=1e-4
    gamma0=1.8; xi0=1.0
    mA=0.5; mS=0.9
    uh=0.5; ul=1.0
    mu=0.0001; hh=0.7
    f=0.1
    ubar = (1.0 - xu) * uh + xu * ul
    beta1=build_beta(n1); beta2=build_beta(n2); beta3=build_beta(n3)
    Adj1=load_adjacency("PPo_adj_A.txt",(n1,n2))
    Adj2=load_adjacency("PPe_adj_A.txt",(n2,n3))
    D1A=Adj1.sum(axis=1); D1P=Adj1.sum(axis=0)
    D2P=Adj2.sum(axis=1); D2S=Adj2.sum(axis=0)
    params=np.array([alpha,kappa,gamma0,xi0,mA,mS,f,ubar,mu,tt,hh])
    caches=dict(beta1=beta1,beta2=beta2,beta3=beta3,
                Adj1=Adj1,Adj2=Adj2,D1A=D1A,D1P=D1P,D2P=D2P,D2S=D2S)
    grid_vals=np.linspace(0.001,1.0,levels)
    combos=list(product(grid_vals,grid_vals,grid_vals))
    num_runs=len(combos)
    if n_workers is None:
        n_workers=int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    n_workers=min(n_workers,48); n_workers=max(1,n_workers)
    ctx=get_context("fork")
    chunksize=max(1, min(1000, math.ceil(num_runs/(n_workers*4))))
    t_start=perf_counter()
    with ctx.Pool(processes=n_workers,
                  initializer=_init_worker,
                  initargs=(params,caches,h,mm,n1)) as pool:
        iterator=pool.imap_unordered(_run_combo, combos, chunksize=chunksize)
        if show_progress:
            iterator=tqdm(iterator,total=num_runs,desc=f"xu={xu:.2f}",dynamic_ncols=True)
        results=list(iterator)
    count=sum(results)
    t_end=perf_counter()
    frac=count/num_runs
    with open(out_path,"a") as fout:
        fout.write(f"{xu:.4f} {count} {num_runs} {frac:.6f}\n")
    print(f"[done] xu={xu:.4f}, count={count}, total={num_runs}, fraction={frac:.6f}, time={t_end-t_start:.1f}s")

# -------------------- Main loop over xu --------------------
if __name__=="__main__":
    out_file="Pollinator_basin_counts_f0.1.txt"
    if os.path.exists(out_file):
        os.remove(out_file)
    # 🔹 loop from 0.01 to 0.1 inclusive with 20 values
    for xu in np.linspace(0.01,1.0,20):
        sweep_one_xu(xu, out_file, levels=30, mm=100000, h=0.01, n_workers=48)
