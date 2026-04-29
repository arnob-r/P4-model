# P4 Model Datasets and Codes

This repository contains the codes and datasets used to simulate the Pollinator–Plant–Pest–Pesticide (P4) model. The model studies how pesticide regulation and farmer-level Integrated Pest Management (IPM) adoption affect pollinator recovery in a tripartite ecological network.

## Dataset

The `data/` folder contains the empirical tripartite network used in the study. It includes:

- `PPo_adj_A.txt`: pollinator–plant adjacency matrix.
- `PPe_adj_A.txt`: plant–pest adjacency matrix.
- `Sinohara_tripartite_adjacency_A.csv`: combined tripartite network data.
- `Contrasting effects of land‐use changes on herbivory.pdf`: reference paper for the empirical network dataset.

These datasets define the ecological interaction structure among pollinators, plants, and pests.

## Codes

The repository contains separate codes for generating the main numerical results:

- `fig_2b/`: computes pollinator extinction-to-persistence transitions under varying IPM adoption \(x_u\) and policy strength \(f\).
- `fig_3a/`: computes the basin of attraction for the focal pollinator species \(A_{18}\).
- `fig_3b_c/`: computes basin stability curves as management parameters change.
- `fig_4/`: analyzes how mutualistic strength \(\gamma_0\) and antagonistic strength \(\xi_0\) affect recovery thresholds.
- `timeseries.py`: generates representative time-series plots for pollinator persistence and extinction dynamics.

## What to Do with the Codes

1. Use the files in `data/` as input ecological networks.
2. Run the Fortran files to generate numerical simulation outputs.
3. Run the Python scripts to visualize the outputs and reproduce the manuscript figures.
4. Modify \(x_u\), \(f\), \(\gamma_0\), \(\xi_0\), and initial conditions to test how pesticide management and ecological interactions affect pollinator recovery.

The codes can be used to reproduce the main figures, explore extinction–persistence transitions, compute basin stability, and study how policy regulation and IPM adoption reshape pollinator recovery in tripartite ecological networks.
