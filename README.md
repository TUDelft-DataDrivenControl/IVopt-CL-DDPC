# Optimal Instrumental Variables for Closed-Loop Data-Driven Predictive Control

This repository provides a complete, reproducible implementation of optimal instrumental variable (IV) selection for data-driven predictive control (DDPC). It contains MATLAB code for simulating, benchmarking, and analyzing the proposed IV-enhanced control algorithms across multiple dynamical systems.

**Paper**: *Optimal Instrumental Variable Selection for Closed-Loop Data-Driven Predictive Control*  
**Authors**: Rogier Dinkla, Tom Oomen, Sebastiaan P. Mulders, J.W. van Wingerden

---

## Quick Start

### 1. Prerequisites
- **MATLAB R2024b** (or compatible recent version) with the following toolboxes:
  - Control System Toolbox
  - Robust Control Toolbox  
  - Statistics and Machine Learning Toolbox
  - System Identification Toolbox

- **CasADi v3.6.7** with IPOPT solver (see installation below)
- **Crameri colormaps** (optional, for paper figure visualization)

### 2. Installation

#### Step 1: Set Up Repository Path
The codebase requires a single setup step to define the repository root directory:

```matlab
% In MATLAB, run:
pdir = 'c:\path\to\IVopt-DDPC\01 - coding';  % Use your actual repo path
save pdir.mat pdir
```

This creates a `pdir.mat` file (already in `src/` if you cloned the repo) that all scripts use to find dependencies.

#### Step 2: Download CasADi

The code requires **CasADi v3.6.7** compiled with IPOPT. Since CasADi cannot be freely redistributed as binaries, users must:

1. Download CasADi from [casadi.org](https://casadi.org/) (version 3.6.7)
2. Extract to: `bin/external/casadi-v3.6.7/`
3. Ensure the compiled MEX file `casadiMEX.mexw64` (Windows) or equivalent is present

**Note**: The repository already includes `.la` library files and supporting files; you only need to add the compiled MEX binaries.

#### Step 3 (Optional): Install Crameri Colormaps

For publication-quality figures with perceptually uniform colormaps:

```matlab
% Download from: https://zenodo.org/records/8409685
% Extract to: bin/external/crameri_colours/
```

If not installed, plots will use default MATLAB colors without error.

---

## Running the Code

### Basic Single Simulation

To run a single Monte Carlo experiment with default settings:

```matlab
cd src/
main()  % or: main(Re=0.0481, p=20, f=20, N=1e4, seed=1)
```

**Typical runtime**: 1-2.5 seconds per run (including I/O)  
**Output**: Data saved to `data/raw/sys6/prbs/<timestamp>/`

### Batch Experiments (Paper Results)

To reproduce the paper's main experimental results:

```matlab
% Experiment 1: Varying data size (N)
main_dN()

% Experiment 2: Varying noise level (Re)
main_dRe()

% Experiment 3: Varying prediction/control horizons (p, f)
main_dp()
```

These scripts run multiple parameter configurations with multiple Monte Carlo replications each. **Expected duration**: 1–4 hours (depending on hardware and `parcluster` settings in the scripts).

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Re` | 4.81e-2 | Innovation noise variance |
| `p` | 20 | Past data horizon length |
| `f` | 20 | Future prediction horizon length |
| `N` | 1e4 | Number of data samples for learning |
| `Ncl` | 1500 | Closed-loop simulation length |
| `seed` | 1 | Random seed for reproducibility |
| `sys` | 6 | System configuration (focus on config 6 = Landau1995 with D=0) |

---

## Project Organization

```
.
├── SUMMARY.md                           <- High-level overview (start here)
├── README.md                            <- This file
├── CITATION.md                          <- How to cite this code
├── LICENSE.md                           <- MIT License
├── src/
│   ├── add_paths.m                      <- Path initialization (auto-loaded by scripts)
│   ├── init_sims.m                      <- System setup & parameter initialization
│   ├── main.m                           <- Single Monte Carlo experiment
│   ├── main_dN.m / main_dRe.m / main_dp.m  <- Batch experiments (varying N, Re, p/f)
│   ├── plot_figs_paper.m                <- Generate all paper figures from processed data
│   ├── @IVopt_DDPC/                     <- Main algorithm: DeePC with optimal IVs
│   ├── @DeePC/                          <- Standard DeePC (baseline)
│   ├── @Generalized_DeePC/              <- Base class for all DeePC variants
│   ├── @CL_DeePC/                       <- Closed-loop DeePC variant
│   ├── fun/                             <- Utility functions (Hankel, Toeplitz, statistics, etc.)
│   ├── analysis/
│   │   ├── process_dX.m                 <- Aggregate Monte Carlo results
│   │   ├── analyze_dX.m                 <- Statistical analysis
│   │   ├── plotting_fun/                <- Figure generation helpers
│   │   └── processing_fun/              <- Data aggregation helpers
│   └── *.mat                            <- Pre-computed initial controllers (Cz0_*.mat)
├── bin/
│   ├── model_Landau1995.m               <- Primary benchmark plant
│   ├── model_Bemporad2002.m             <- Alternative benchmark systems
│   ├── model_Favoreel1999.m
│   ├── model_Wang2023.m
│   ├── model_DoubleTank.m
│   └── external/
│       └── casadi-v3.6.7/               <- CasADi binaries (user must download)
├── data/ (git-ignored)
│   ├── raw/                             <- Individual Monte Carlo run outputs
│   │   └── sys6/prbs/<timestamp>/       <- Organized by system, reference, timestamp
│   └── processed/                       <- Aggregated results per parameter
├── results/ (git-ignored)
│   └── figures/                         <- Paper figures (PDF, PNG)
└── pdir.mat                             <- Repository root path (created on setup)
```

---

## Reproducing Paper Results

### Step-by-Step Guide

1. **Setup** (one-time only):
   ```matlab
   % Set pdir and download CasADi as described in "Installation" section above
   pdir = '...';
   save pdir.mat pdir
   ```

2. **Run batch experiments** (this takes 1-4 hours):
   ```matlab
   cd src/
   main_dN()   % Collect results for varying N
   main_dRe()  % Collect results for varying Re
   main_dp()   % Collect results for varying p, f
   ```

3. **Process results**:
   ```matlab
   cd src/analysis/
   process_dX()  % Aggregates Monte Carlo data for all experiments
   ```

4. **Generate figures** (matches paper exactly):
   ```matlab
   cd src/
   plot_figs_paper()  % Creates all paper figures in results/figures/
   ```

### Expected Outputs

- **Raw data**: `data/raw/sys6/<timestamp>/` — Contains subdirectories for each parameter setting, each with `≈50–100` Monte Carlo runs
- **Processed data**: `data/processed/` — Summary statistics (mean, std, percentiles) for each experiment
- **Figures**: `results/figures/` — Publication-ready `.pdf` and `.png` files for all paper plots

---

## Key Experimental Details

### System Configuration (sys=6, used throughout)
- **Plant**: Landau 1995 5th-order model with input delay
- **Features**: Proper system (D=0), realistic input/output constraints
- **Nominal controller**: Designed via mixed-sensitivity (H∞) synthesis
- **Noise**: Additive Gaussian innovation noise (variance = Re)

### Algorithms Compared
1. **IVopt-DDPC** — Proposed method with optimal instrumental variables
2. **Standard DeePC** — Baseline without IV weighting
3. **Initial Controller** — Existing plant model-based controller

### Performance Metrics
- **Prediction quality**: Variance of one-step-ahead output prediction errors
- **Closed-loop performance**: Reference tracking error + input energy
- **Data efficiency**: Performance vs. number of training samples (N)
- **Noise robustness**: Performance vs. noise level (Re)

---

## Advanced Usage

### Changing System or Parameters

Modify the system and initial conditions in `init_sims.m`:

```matlab
opts.sys  = 6;          % Choose plant model (1-8, see init_sims.m for details)
opts.Re   = 1e-2;       % Noise variance
opts.p    = 20;         % Past horizon
opts.f    = 20;         % Future horizon
opts.N    = 1e4;        % Number of data samples
opts.seed = 42;         % Reproducibility
opts.plot = true;       % Visualize closed-loop trajectories
```

### Parallel Execution on Clusters

The scripts `DB_dN.sh`, `DB_dRe.sh` are SLURM job submission scripts for high-performance clusters (e.g., Delft BlueGene). Edit and submit directly:

```bash
sbatch DB_dN.sh
```

### Accessing Pre-computed Controllers

Initial closed-loop controllers are stored as:
- `Cz0_Landau1995_D0.mat` — For system configuration 6 (primary)
- `Cz0_*.mat` — For other systems

These are generated automatically by `init_sims.m` if missing (uses H∞ mixed-sensitivity design).

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `CasADi not found` | Ensure `bin/external/casadi-v3.6.7/casadiMEX.mexw64` exists |
| `pdir.mat not found` | Run `pdir = '...'; save pdir.mat pdir` in MATLAB in the `src/` directory |
| `Missing toolboxes` | Verify all 4 MATLAB toolboxes are installed: `ver` command lists installed toolboxes |
| `Out of memory` | Reduce `N` or `Ncl` parameters; run experiments in smaller batches |
| `Permission denied (bin/external/)` | CasADi files may be read-only; ensure write access to workspace |

---

## Dependencies & Licenses

### MATLAB Toolboxes (Proprietary)
- Control System Toolbox (R2024b)
- Robust Control Toolbox (R2024b)
- Statistics and Machine Learning Toolbox (R2024b)
- System Identification Toolbox (R2024b)

### External Open-Source
- **CasADi** — GNU Lesser General Public License v3 (https://github.com/casadi/casadi)
- **IPOPT** — Eclipse Public License v1.0 (bundled with CasADi)
- **Crameri Colormaps** — https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps

### This Repository
- **MIT License** — See [LICENSE.md](LICENSE.md)

---

## Contact & Support

For questions about implementation or reproducibility:
- **Repository Issues**: Use the GitHub issue tracker
- **Paper Questions**: Contact the authors via Delft University of Technology

---

## License

This project is licensed under the terms of the [MIT License](LICENSE.md)
