## Project Overview
The code in this project underlies the numerical experiments from the article<br>
**"Optimal Instrumental Variable Selection for Closed-loop Data-Driven Predictive Control"** by <br>
Rogier Dinkla<sup>1</sup>, Tom Oomen<sup>1,2</sup>, Sebastiaan P. Mulders<sup>1</sup>, and Jan-Willem van Wingerden<sup>1</sup>.<br>

*Affiliations*:<br>
<sup>1</sup> Delft Center for Systems and Control, Faculty of Mechanical Engineering, Delft University of Technology, The Netherlands<br>
<sup>2</sup> Control Systems Technology Group, Department of Mechanical Engineering, Eindhoven University of Technology, The Netherlands

The above work derives an optimal instrumental variable (IV) for closed-loop Data-Driven Predictive Control (DDPC) that minimizes the asymptotic variance of employed system estimates that are asymptotically unbiased. Furthermore, approximations of optimal IVs are provided that either do or do not rely on controller knowledge.

The primary intent of the code is to facilitate reproduction of results presented in the paper, which stem from a batch of Monte Carlo simulations that perform a parametric sweep over a range of innovation noise variances (Re). The code can furthermore be used to generate results from a wider variety of systems and initial controllers, and perform batches of Monte Carlo simulations that sweep over the past window length (p) or number of data matrix columns (N). Since the analysis of results that were not presented in the paper is not the primary purpose of this repository, manual adjustments may be needed to visualize other data.

## Main Components

### 1. **Monte Carlo Simulation Engine** (`src/MC_sim/`)
- `main_MC.m` — Orchestrates a full Monte Carlo run: generates data, computes IVs, builds predictors, runs closed-loop simulation
- `CaseDefinitions.m` — Defines 15 controllers (12 IV-based DDPC variants + 3 benchmarks)
- `IV_4_DDPC.m` — class to efficiently store the similarly structured IVs
- `get_Z.m` — Creates IV matrices for all IV-based DDPC controller cases.
- `get_Lf_Cz.m` — Generates (estimates of) the matrix $L_f$, and uses this to generate an SPC-type of controller by means of `Lf_2_SPC.m`.
- `get_solver.m` — function to find how the optimal input from unconstrained SPC depends on $L_f$, past IO data, and future reference signals.
- `approx_IV_no_controller_info.m`, `approx_IV_controller_info.m` — Respectively Algorithms 1 and 2 from the paper. Used to approximate optimal IVs.

### 2. **Simulation Entry Points** (`src/`)
- `main.m` — Single Monte Carlo simulation with configurable parameters (Re, N, p, f, seed, sys, ...).
- `main_dRe.m` — Batch experiments sweeping innovation noise variance (Re). Single MC runs performed by `run_Re.m`.
- `main_dp.m` — Batch experiments sweeping past horizon length (p). Single MC runs performed by `run_p.m`.
- `main_dN.m` — Batch experiments sweeping number of data matrix columns (N). Single MC runs performed by `run_N.m`.

### 3. **Plant Models & Controllers** (`src/systems/`)
- `get_sys_info.m` — Route to the selected plant model & initial controller based on `sys` parameter
- `init_sims.m` — Initialize simulation: select plant model and initial controller using `get_sys_info.m` and construct initial closed-loop system.
- `tune_Cz0.m` — Tune initial controller `Cz0` via H-infinity mixed-sensitivity design
- Details of the system configurations are available in directories that are named after the used plant model:
  | sys | Plant/Directory | Initial Controller (`Cz0`) | Notes |
  |-----|-------|------------|-------|
  | 1 | Landau1995 [<a href="#ref29">29</a>,<a href="#ref30a">30</a>] | `Cz0_Landau1995_D0.mat` | **Results presented in paper.**<br> 4th-order system with input-output delay,<br> 5th order `Cz0` without direct feedthrough |
  | 2 | Landau1995 [<a href="#ref29">29</a>,<a href="#ref30a">30</a>] | `Cz0_Landau1995_D0_n50.mat` | 50th order `Cz0` without direct feedthrough |
  | 3 | Landau1995 [<a href="#ref29">29</a>,<a href="#ref30a">30</a>] | `Cz0_Landau1995.mat` | 5th order `Cz0` with direct feedthrough ($D_c$ ≠ 0) |
  | 4 | Bemporad2002 [<a href="#ref32">32</a>] | `Cz0_Bemporad2002.mat` | 2nd-order unstable system, 4th order `Cz0`|
  | 5 | Favoreel1999 [<a href="#ref33">33</a>] | `Cz0_Favoreel1999.mat` | Marginally stable 5th-order system, 7th order `Cz0` |
  | 6 | Wang2023 [<a href="#ref8">8</a>] | `Cz0_Wang2023.mat` | Unstable 3rd-order plant, 5th order `Cz0` |
  | 7 | Wang2023 [<a href="#ref8">8</a>] | `Cz0_Wang2023_provided.mat` | Unstable 3rd-order plant, `Cz0` of 2nd order provided by authors and `.mat` file generated from `model_Wang2023.m` |

### 4. **Data Processing** (`src/processing/`)
- `main_processing.m` — Aggregates raw Monte Carlo batch results into summary statistics, saves these to `processed_data.mat` files. Assumes 'Re' type data unless otherwise specified in script.
- `process_dX.m` — Batch aggregation logic (computes m1, mLf, m4 metrics across seeds)
- `util_fun/` — Helper functions for processing

### 5. **Visualization** (`src/plot_figs/`)
- `plot_figs_paper.m` — Generates the figures from the paper.
- `Fig_IV_approx.m`, `Fig_Lf_estimates.m`, `Fig_prediction_quality.m`, `Fig_sim_example.m` — Functions to help plot individual figures.

### 6. **Utility Functions** (`src/util/`)
- `blk_toeplitz_mean.m` — Performs averaging over block-diagonals corresponding to a block-Toeplitz structure.
- `make_Hankel.m` — Hankel matrix construction from data.
- `make_TrajMat.m` — Trajectory matrix construction from data.
- `get_subdir1.m` — Helps navigation to the right data directory. Called by `main_processing.m`.
- `make_blk_tril.m` — Block-lower-triangularize provided matrix.
- `make_blk_tril_toeplitz.m` — Block Toeplitz/triangular matrix operations.
- `make_ext_ctrb.m`, `make_ext_obsv.m` — Extended controllability/observability matrices.
- `make_reference.m` — Replicates reference from [<a href="#ref30a">30</a>]
- `plant2ABCDK.m` — Gets A, B, C, D, and K matrices from the specified plant.
- `ss2lag.m` — Computes the lag of a state-space system.
- `findCases2Sim.m` — To obtain cases to simulate.

## Simulated Controllers (15 cases)

| ID | Type | Description |
|----|------|-------------|
| `iv1` | Baseline | Open-loop IV |
| `iv2` | Optimal (exact) | Denoised IV |
| `iv3` | LCF-IV [<a href="#ref8">8</a>] | IV based on the left coprime factorization  |
| `iv4a`, `iv4b`, `iv4c`, `iv4d` | Approx. (no `Cz0` info) | Approximated optimal IV without controller knowledge (2SLS, 2SLS + causal + time invariance, row-by-row, row-by-row + time invariance) |
| `iv5a`, `iv5b`, `iv5c`, `iv5d` | Approx. (with `Cz0` info) | Approximated optimal IV with controller knowledge (2SLS, 2SLS + causal + time invariance, row-by-row, row-by-row + time invariance) |
| `iv6`  | Reference-based | Future reference as IV |
| `CLSPC` | Benchmark | Closed-loop SPC |
| `actLf` | Oracle | SPC with true Lf |
| `TrPred` | Transient | Transient predictor |

## License
This code is released under the **MIT License** (see [LICENSE.md](LICENSE.md)).

## References
<a id="ref8"></a>[8] Y. Wang, Y. Qiu, M. Sader, D. Huang, and C. Shang, "Data-Driven Predictive Control Using Closed-Loop Data: An Instrumental Variable Approach," *IEEE Control Systems Letters*, vol. 7, pp. 3639–3644, 2023, doi: [10.1109/LCSYS.2023.3340444](https://doi.org/10.1109/LCSYS.2023.3340444).<br>
<a id="ref29"></a>[29] I. D. Landau, D. Rey, A. Karimi, A. Voda, and A. Franco, "A Flexible Transmission System as a Benchmark for Robust Digital Control," *European Journal of Control*, vol. 1, no. 2, pp. 77–96, Jan. 1995, doi: [10.1016/S0947-3580(95)70011-5](https://doi.org/10.1016/S0947-3580(95)70011-5).<br>
<a id="ref30a"></a>[30] A. Chiuso, M. Fabris, V. Breschi, and S. Formentin, "Harnessing uncertainty for a separation principle in direct data-driven predictive control," *Automatica*, vol. 173, p. 112070, Mar. 2025, doi: [10.1016/j.automatica.2024.112070](https://doi.org/10.1016/j.automatica.2024.112070).<br>
<a id="ref32"></a>[32] A. Bemporad, M. Morari, V. Dua, and E. N. Pistikopoulos, "The explicit linear quadratic regulator for constrained systems," *Automatica*, vol. 38, no. 1, pp. 3–20, Jan. 2002, doi: [10.1016/S0005-1098(01)00174-1](https://doi.org/10.1016/S0005-1098(01)00174-1).<br>
<a id="ref31"></a>[33] W. Favoreel, B. De Moor, M. Gevers, and P. Van Overschee, "Closed-Loop Model-Free Subspace-Based LQG-Design," in *Proceedings of the Mediterranean Conference on Control and Automation*, Jan. 1999.