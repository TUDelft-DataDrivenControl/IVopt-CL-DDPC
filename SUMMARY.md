# SUMMARY: Optimal Instrumental Variables for Closed-Loop Data-Driven Predictive Control

## Project Overview

This repository contains the complete implementation of **optimal instrumental variable (IV) selection for Closed-loop Data-Driven Predictive Control (DDPC)**. The code enables researchers to reproduce the theoretical results and numerical experiments from the paper *"Optimal Instrumental Variable Selection for Closed-Loop Data-Driven Predictive Control"*.

The work develops both exact and approximate formulations of optimal instrumental variables specifically designed for data-driven predictive control in closed-loop scenarios, with validated implementations across multiple benchmark dynamical systems.

## Main Components

### 1. **Core Algorithms**
- `@Generalized_DeePC/` — Foundation class implementing data-enabled predictive control (DeePC) framework
- `@DeePC/` — Standard DeePC implementation without instrumental variables
- `@IVopt_DDPC/` — **Main contribution**: DeePC with optimal instrumental variable selection and adaptive IV refinement
- `@CL_DeePC/` — Closed-loop DeePC variant
- **Key files**: `make_solver.m` and `optimizer_solve.m` handle CasADi-based optimization formulation and solution

### 2. **Simulation Framework**
- **Entry points**:
  - `main.m` — Single Monte Carlo simulation with configurable parameters
  - `main_dN.m` — Batch experiments varying data matrix size (N)
  - `main_dRe.m` — Batch experiments varying noise level (innovation noise variance)
  - `main_dp.m` — Batch experiments varying prediction/control horizons (p, f)
- **Initialization**: `init_sims.m` sets up plant models, controllers, and simulation parameters
- **Data management**: Organized into raw/processed data directories with Monte Carlo replication tracking

### 3. **Plant Models**
- `model_Landau1995.m` — **Primary benchmark**: 5th-order system with input delay (used in main paper results)
- `model_Bemporad2002.m` — 2nd-order Unstable system with input-output characteristics
- `model_Favoreel1999.m` — 5th-order system for subspace predictive control comparison
- `model_Wang2023.m` — Additional benchmark system
- `model_DoubleTank.m` — Physical two-tank system configuration

### 4. **Analysis & Visualization**
- `analysis/process_dX.m` — Post-processing of experimental results (batches of Monte Carlo runs)
- `analysis/analyze_dX.m` — Statistical analysis and aggregation
- `plot_figs_paper.m` — **Primary script**: Reproduces all paper figures from processed data
- `analysis/plotting_fun/` — Supporting plotting utilities
- `analysis/processing_fun/` — Statistical processing functions

### 5. **Utility Functions**
- **Data handling**: `get_Z.m`, `get_uy_iv.m`, `get_subdir1.m`
- **Mathematical operations**: `make_Hankel.m`, `make_Page.m`, `make_ext_ctrb.m`, `make_ext_obsv.m`, `make_blk_tril_toeplitz.m`
- **Control theory** : `plant2ABCDK.m`, `Lf_2_SPC.m`, `ss2lag.m`, `constrained_LS.m`
- **Numerical analysis**: `calc_var_mean_est.m`, `diag_stats.m`, `make_boxplot.m`

## Data Flow

```
1. SETUP & INITIALIZATION
   └─ init_sims.m (configure system, controller, parameters)
   
2. MONTE CARLO EXPERIMENTS
   └─ main.m / main_dN.m / main_dRe.m / main_dp.m
      ├─ Generate/load closed-loop data
      ├─ Run IVopt-DDPC and DeePC algorithms
      ├─ Evaluate prediction quality and control performance
      └─ Save results to data/sys<#>/<ref>/<category>/

3. POST-PROCESSING
   └─ analysis/process_dX.m
      ├─ Aggregate Monte Carlo batches
      ├─ Compute statistics (mean, std, percentiles)
      └─ Save to data/processed/

4. VISUALIZATION & ANALYSIS
   └─ plot_figs_paper.m
      ├─ Generates publication-ready figures
      └─ Saves to results/figures/
```

## Key Technical Contributions

1. **Optimal Instrumental Variables**: Derives closed-form expressions for optimal IVs in DDPC that minimize prediction variance
2. **Adaptive Refinement**: Implements iterative IV refinement to progressively improve estimation quality
3. **Computational Efficiency**: Provides approximations of exact optimal IVs for reduced computational burden
4. **Closed-Loop Validation**: Demonstrates performance on standard benchmark systems with varying noise and data availability

## Experimental Parameters

The code explores three main dimensions:

| Dimension | Script | Range | Purpose |
|-----------|--------|-------|---------|
| **Data size** | `main_dN.m` | N = 50–10,000 | Assess impact of training data availability |
| **Noise level** | `main_dRe.m` | Re = 10⁻⁴–10⁻¹ | Evaluate robustness to measurement noise |
| **Horizons** | `main_dp.m` | p, f = 5–50 | Test sensitivity to prediction/control window |

Each dimension runs multiple Monte Carlo replications (default: ~50–100 per parameter set) to quantify variability.

## External Dependencies

### **MATLAB Toolboxes** (R2024b or compatible)
- Control System Toolbox
- Robust Control Toolbox
- Statistics and Machine Learning Toolbox
- System Identification Toolbox

### **External Software**
- **CasADi v3.6.7** — Symbolic differentiation and optimization (users must download; see README)
- **IPOPT** — Nonlinear programming solver (bundled with CasADi)
- **Crameri colormaps** (optional) — Perceptually uniform scientific color schemes for visualization (v1.09 from Zenodo)

### **Included External Code**
- **crameri_colours/** — Color palette utilities (if already present)

## Output Artifacts

- **Raw results**: `data/sys<#>/<reference_type>/<category>/<timestamp>/` — Individual Monte Carlo runs
- **Processed data**: `data/processed/` — Aggregated statistics and summary tables
- **Figures**: `results/figures/` — All paper figures as `.pdf`, `.png` (generated by `plot_figs_paper.m`)
- **Pre-computed data**: `*.mat` files in `src/` contain cached initial controllers (e.g., `Cz0_Landau1995.mat`)

## System Configuration

The README focuses on **System Configuration 6** (defined in `init_sims.m`), which uses:
- **Plant**: Landau1995 model with D=0 (proper system)
- **Control objective**: Minimize output tracking error and input effort
- **Noise**: Additive innovation noise with tunable variance (Re parameter)

## Citation & License

This code is released under the **MIT License** (see LICENSE.md). For citations, please refer to CITATION.md.

---

**Authors**: Rogier Dinkla, Tom Oomen, Sebastiaan P. Mulders, J.W. van Wingerden
**Institution**: Delft University of Technology
**Last Updated**: 2026