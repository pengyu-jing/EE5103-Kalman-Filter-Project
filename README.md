# EE5103 Project: Kalman Filter (Jing Pengyu, A0326293B)

MATLAB code and report for the Kalman filter project of EE5103 Computer Control Systems (NUS): Part 1 (random-walk target), Part 2 (constant-velocity target) and Part 3 (indoor positioning of a target in circular motion).

- `EE5103_Project_Report.pdf`: the report (IEEE conference template).
- `matlab/`: MATLAB code (base MATLAB only, no toolbox needed).
  - `run_all.m`: runs Parts 1–3, prints every number of the report, saves the figures to `matlab/figures/` and the Part 3 tables to `matlab/results/part3_tables.txt`. At the end it prints a self-check against the report (`PASS`).
  - `part1.m`, `part2.m`, `part3.m`: the three parts (Graphs 1–9, bias/variance, Figure 5.1 with the estimates, tables of Kf(k), K(k), P(k|k), P(k+1|k), Monte Carlo studies).
  - `kf.m`: the Kalman filter, eqs. (1.52)–(1.57) of the notes. `kf_gains.m`: the same filter with precomputed gains, used for the Monte Carlo runs.
  - `kf_steady.m`: steady-state Riccati solution. `kf_freq.m`: H(e^jw) and the H2 norm.
  - `cv_model.m`, `osc_model.m`: the models. `identify_circle.m`: circle and angular-rate identification. `save_fig.m`: figure export.

## Running the code

Open MATLAB in the `matlab` folder and run `run_all`.

The figures and tables in the report were produced with an equivalent Python implementation of the same algorithm, with fixed random seeds. Part 3 uses the measured data, so the MATLAB output is identical to the report. Parts 1 and 2 use random noise. MATLAB's random-number generator differs from Python's, so their simulated bias and variance differ slightly from Table I, within the run-to-run spread. The analytic values (steady-state gains and covariances) are identical.
