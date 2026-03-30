# Optimal Instrumental Variables for Closed-Loop Data-Driven Predictive Control

This repository provides a complete, reproducible implementation of optimal instrumental variable (IV) selection for data-driven predictive control (DDPC). It contains MATLAB code for simulating, benchmarking, and analyzing the proposed IV-enhanced control algorithms across multiple dynamical systems.

**Paper**: *Optimal Instrumental Variable Selection for Closed-loop Data-Driven Predictive Control*  
**Authors**: Rogier Dinkla<sup>1</sup>, Tom Oomen<sup>1,2</sup>, Sebastiaan P. Mulders<sup>1</sup>, and Jan-Willem van Wingerden<sup>1</sup><br>
*Affiliations*:<br>
<sup>1</sup> Delft Center for Systems and Control, Faculty of Mechanical Engineering, Delft University of Technology, The Netherlands<br>
<sup>2</sup> Control Systems Technology Group, Department of Mechanical Engineering, Eindhoven University of Technology, The Netherlands<br>

**License**: MIT (see [LICENSE.md](LICENSE.md))

---
## Prerequisites
### Required Software:
- **MATLAB R2024b** (or compatible) with toolboxes:

  | Toolbox | Version |
  |---------|---------|
  | Control System Toolbox | v24.2 |
  | Robust Control Toolbox | v24.2 |
  | Statistics and Machine Learning Toolbox | v24.2 |
  | System Identification Toolbox | v24.2 |
  | Parallel Computing Toolbox | v24.2 |
   - **Check installation**: Run `ver` in MATLAB to verify toolbox availability
   - **License**: Commercial (MATLAB license required)
- **CasADi v3.6.7**
  - License: [LGPL v3.0](https://github.com/casadi/casadi/blob/3.6.7/LICENSE.txt) (open source)
  - Purpose: Used for symbolic computations. Employed here to obtain a SPC-type of controllers.
- **Crameri colormaps v7.0+**
  - License: MIT (open source, see [MATLAB file exchange](https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps/files/README.md) -> [Zenodo](https://doi.org/10.5281/zenodo.1243862))
  - Purpose: Used for plotting with perceptually uniform colormaps.

## Installation

### Automated installation using `INIT.m`:
```matlab
cd src/
INIT  % Downloads CasADi & Crameri, configures paths, optionally sets up SLURM
```
- Edit SLURM cluster configuration settings in `INIT.m` before running (optional — only needed if using a SLURM cluster)

### Manual installation:

##### Step 1: Add the src directory to the path
From the project directory in MATLAB, run:
```matlab
addpath(genpath('src'));
```

##### Step 2: SLURM Configuration (optional)
Specify the following parameters in MATLAB and save them to `src/SlurmSettings.mat`:
```matlab
ProfileName = ''; % Enter profile name used for the SLURM cluster (e.g., 'MyClusterProfile')
account     = ''; % Enter account e.g. 'institution-department'
partition   = ''; % Enter partitions to use, comma-separated if multiple (e.g., 'partition1,partition2')
save(fullfile('src','SlurmSettings.mat'),'account','partition','ProfileName');
```
**N.B.**: The above parameters may be left empty if no SLURM cluster is used, but the file must be saved as shown.

##### Step 3: Download CasADi
The code requires **CasADi v3.6.7**:
1. Download CasADi from [casadi.org](https://web.casadi.org/get/) (version 3.6.7)
2. Extract the contents to: `bin/casadi-v3.6.7/`
3. Add the directory to the path:
   ```matlab
   addpath(fullfile('bin','casadi-v3.6.7'));
   ```
4. Verify by running:
   ```matlab
   x = casadi.SX.sym('x');
   fprintf('CasADi loaded successfully.\n');
   ```

##### Step 4: Install Crameri Colormaps
1. Download [Crameri colormaps](https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)
2. Extract contents to: `bin/crameri_colours/`
3. Add to path:
   ```matlab
   addpath(fullfile('bin','crameri_colours'));
   ```
4. Verify by running `crameri` in MATLAB (opens an illustration of available colormaps)

---

## Running the Code

### Basic Single Simulation

To run a single Monte Carlo experiment with default settings:

```matlab
cd src/
main()  % or: main(Re=0.01, p=20, f=20, N=1000, seed=1)
```

**Typical runtime**: 1–2.5 seconds per run (including I/O)  
**Output**: Data saved to `data/sys1/ref0_prbs/Re_1e-02_N_1e03_p_20_f_20_Ncl_1500_Qk_100_Rk_1_dRk_1/seed_1.mat`

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Re` | 0.01 | Innovation noise variance |
| `N` | 1000 | Number of columns in Hankel data matrix |
| `p` | 20 | Past horizon (learning window length) |
| `f` | 20 | Future horizon (prediction/control window) |
| `Ncl` | 1500 | Closed-loop simulation length |
| `seed` | 1 | Random seed for reproducibility |
| `sys` | 1 | System configuration (1 = Landau1995, primary) |
| `ref0` | 'prbs' | Reference signal type: 'make' or 'prbs' |
| `Qk`, `Rk`, `dRk` | 100, 1, 1 | Controller cost weights |

### Batch Experiments (Paper Results)

To reproduce the paper's main experimental results:

```matlab
cd src/
main_dRe()  % Sweep innovation noise variance (Re)
main_dp()   % Sweep past horizon (p)
main_dN()   % Sweep number of Hankel columns (N)
```

These scripts run multiple parameter configurations with 100 Monte Carlo seeds each.

---

## Reproducing Paper Results

### Step-by-Step Guide

1. **Setup** (one-time only):
   ```matlab
   cd src/
   INIT  % Or follow manual installation steps above
   ```

2. **Run batch experiments**:
   ```matlab
   cd src/
   main_dRe()  % Primary: varying noise level
   ```

3. **Process results**:
   ```matlab
   cd src/processing/
   main_processing  % Aggregates Monte Carlo data, saves processed_data.mat
   ```

4. **Generate figures** (matches paper exactly):
   ```matlab
   cd src/plot_figs/
   plot_figs_paper  % Creates publication-ready figures in results/
   ```

### Expected Outputs

- **Raw data**: `data/sys1/ref0_prbs/dRe/<batch_name>/` — Subdirectories per parameter value, each with 100 seed files
- **Processed data**: `processed_data.mat` — Summary statistics (mean, std, median, percentiles) saved in the batch directory
- **Figures**: `results/` — Publication-ready PDF files (`dinkl2.pdf` through `dinkl5.pdf`)

---

## Project Organization

```
.
├── SUMMARY.md                              ← High-level overview of components
├── README.md                               ← This file: installation & reproduction guide
├── CITATION.md                             ← How to cite this code
├── LICENSE.md                              ← MIT License
│
├── src/
│   ├── INIT.m                              ← One-time setup (CasADi, Crameri, paths, SLURM)
│   ├── main.m                              ← Single Monte Carlo experiment
│   ├── main_dRe.m / main_dp.m / main_dN.m ← Batch experiments (varying Re, p, N)
│   ├── run_Re.m / run_p.m / run_N.m       ← Batch runner helpers
│   ├── run_X_ParCluster.m                  ← SLURM parallel cluster submission
│   ├── SlurmSettings.mat                   ← SLURM cluster configuration
│   │
│   ├── MC_sim/                             ← Monte Carlo simulation engine
│   │   ├── main_MC.m                       ← Orchestrates full MC run
│   │   ├── CaseDefinitions.m              ← Defines 17 IV cases
│   │   ├── IV_4_DDPC.m                    ← IV matrix storage class
│   │   ├── get_Z.m                        ← Creates IV matrices
│   │   ├── get_Lf_Cz.m                   ← Generates SPC predictors & controllers
│   │   ├── get_solver.m                   ← CasADi optimization solver
│   │   ├── get_actual_matrices.m          ← Oracle/true matrices
│   │   ├── Lf_2_SPC.m                    ← Convert Lf to SPC controller
│   │   ├── approx_IV_no_controller_info.m ← Approximate IVs (no controller info)
│   │   └── approx_IV_controller_info.m    ← Approximate IVs (with controller info)
│   │
│   ├── processing/                         ← Stage B: Data aggregation
│   │   ├── main_processing.m              ← Entry point
│   │   ├── process_dX.m                   ← Batch aggregation logic
│   │   └── util_fun/                      ← Helper functions
│   │
│   ├── plot_figs/                          ← Stage C: Visualization
│   │   ├── plot_figs_paper.m              ← Generates all paper figures
│   │   ├── Fig_IV_approx.m               ← IV approximation comparison figure
│   │   ├── Fig_Lf_estimates.m            ← Lf matrix estimate figure
│   │   ├── Fig_prediction_quality.m      ← Prediction quality figure
│   │   └── Fig_sim_example.m             ← Simulation example figure
│   │
│   ├── systems/                            ← Plant models and controllers
│   │   ├── init_sims.m                    ← Initialize simulation
│   │   ├── get_sys_info.m                 ← Route to plant model
│   │   ├── tune_Cz0.m                    ← Tune initial controller
│   │   ├── Landau1995/                    ← Primary benchmark (3 controllers)
│   │   ├── Bemporad2002/                  ← 2nd-order unstable system
│   │   ├── Favoreel1999/                  ← 5th-order system
│   │   └── Wang2023/                      ← Additional benchmark
│   │
│   └── util/                               ← General utilities
│       ├── make_Hankel.m                  ← Hankel matrix construction
│       ├── make_blk_tril_toeplitz.m       ← Block Toeplitz matrices
│       ├── make_blk_tril.m               ← Block triangular matrices
│       ├── make_ext_ctrb.m / make_ext_obsv.m ← Extended matrices
│       ├── plant2ABCDK.m                  ← State-space extraction
│       └── ss2lag.m                       ← Controller lag
│
├── bin/                                    ← External dependencies (user must download)
│   ├── casadi-v3.6.7/                     ← CasADi optimization library
│   └── crameri_colours/                   ← Scientific colormaps
│
├── data/ (git-ignored)                     ← Simulation results
│   └── sys1/ref0_prbs/                    ← Organized by system, reference type
│       ├── dRe/<batch>/                   ← Batch: varying noise level
│       ├── dp/<batch>/                    ← Batch: varying past horizon
│       └── dN/<batch>/                    ← Batch: varying data size
│
└── results/ (git-ignored)                  ← Generated paper figures (PDFs)
```

---

## Key Experimental Details

### System Configuration (sys=1, primary paper results)
- **Plant**: Landau 1995 5th-order model with input delay
- **Initial controller**: 5-state controller designed via mixed-sensitivity (H∞) synthesis
- **Noise**: Additive Gaussian innovation noise (variance = Re)

### System Configurations

| sys | Plant | Controller | Notes |
|-----|-------|------------|-------|
| 1 | Landau1995 | 5-state initial controller | **Primary paper results** |
| 2 | Landau1995 | 50-state initial controller | High-order comparison |
| 3 | Landau1995 | 5-state, Dc ≠ 0 | Direct feedthrough |
| 4 | Bemporad2002 | — | Unstable system |
| 5 | Favoreel1999 | — | Alternative benchmark |
| 6 | Wang2023 | — | Recent benchmark |
| 7 | Wang2023 | — | Uses provided controller |

### Controllers Compared (17 cases)

| ID | Type | Description |
|----|------|-------------|
| `iv1` | Baseline | Open-loop IV |
| `iv2a`, `iv2b`, `iv2c` | Optimal (exact) | Optimal IV with 0, 1, 2 refinement iterations |
| `iv3a`, `iv3c` | LCF-IV | Linear conjugate form IV |
| `iv4a`, `iv4b`, `iv4c` | Approx (no info) | Approximated optimal IV without controller knowledge |
| `iv5a`, `iv5b`, `iv5c` | Approx (with info) | Approximated optimal IV with controller knowledge |
| `iv6a`, `iv6c` | Reference-based | Future reference as IV |
| `CLSPC` | Benchmark | Standard closed-loop SPC |
| `actLf` | Oracle | True transfer matrix (upper bound) |
| `TrPred` | Transient | Transient predictor |

IVs with "b" suffix incorporate (approximations of) future denoised outputs.  
IVs with "c" suffix are variants that employ 2SLS (two-stage least squares).

---

## SLURM Cluster Execution

For large batch experiments on HPC clusters:

1. Configure SLURM settings in `INIT.m` (or create `src/SlurmSettings.mat` manually)
2. Running `main_dRe`, `main_dp`, or `main_dN` automatically:
   - Detects SLURM settings from `SlurmSettings.mat`
   - Submits jobs via `run_X_ParCluster.m`
   - Manages job runtime and returns results to standard data directory

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| CasADi not found | Ensure `bin/casadi-v3.6.7/casadiMEX.mexw64` exists and is on the MATLAB path |
| Missing toolboxes | Run `ver` in MATLAB to check which toolboxes are installed |
| Out of memory | Reduce `N` or `Ncl` parameters; run experiments in smaller batches |
| SLURM jobs fail | Check `SlurmSettings.mat` has correct `ProfileName`, `account`, `partition` |

---

## Software Dependencies & Licenses

### **MATLAB (Required)**
- **Version**: R2024b (or compatible)
- **License**: Commercial (proprietary)
- **Toolboxes**: Control System, Robust Control, Statistics & ML, System Identification, Parallel Computing (all v24.2)

### **CasADi (Required)**
- **Version**: v3.6.7
- **License**: LGPL v3.0 (open source, see [CasADi license](https://github.com/casadi/casadi/blob/master/LICENSE.txt))
- **Includes**: IPOPT nonlinear solver (Eclipse Public License v1.0)
- **Purpose**: Symbolic differentiation and nonlinear optimization
- **Installation**: Download from [casadi.org](https://web.casadi.org/get/)

### **Crameri Colormaps (Required for paper figures)**
- **License**: MIT (open source)
- **Source**: [MathWorks File Exchange](https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)
- **Purpose**: Publication-quality scientific colormaps

### **This Code Repository**
- **License**: MIT License (see [LICENSE.md](LICENSE.md))
- **Copyright**: (c) 2024 Data-Driven Control, TU Delft
- **Permissions**: Freely usable for research and commercial purposes with attribution

---

## Hardware Requirements

### **Minimum Specification**
- **Processor**: Any modern multi-core CPU (≥2 cores)
- **RAM**: 4 GB (minimum); 8 GB recommended
- **Disk**: 10 GB free space (for raw and processed data)
- **GPU**: Not required; computation uses CPU

### **Recommended Specification for Batch Experiments**
- **Processor**: ≥4 cores (for parallel Monte Carlo runs via MATLAB `parpool`)
- **RAM**: ≥16 GB (for full batch sweeps)
- **Disk**: ≥50 GB (for large-scale parameter sweeps)

### **Specialized Equipment**
- None required. All simulations are purely computational.

---

## Contact & Support

For questions about implementation or reproducibility:
- **Repository Issues**: Use the GitHub issue tracker
- **Paper Questions**: Contact the authors via Delft University of Technology

---

## License

This project is licensed under the terms of the [MIT License](LICENSE.md)
