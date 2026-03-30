# Optimal Instrumental Variable Selection for Closed-loop Data-Driven Predictive Control

This repository provides the code code that was employed to arrive at the results presented in *Optimal Instrumental Variable Selection for Closed-loop Data-Driven Predictive Control* by Rogier Dinkla<sup>1</sup>, Tom Oomen<sup>1,2</sup>, Sebastiaan P. Mulders<sup>1</sup>, and Jan-Willem van Wingerden<sup>1</sup>.<br>

*Affiliations*:<br>
<sup>1</sup> Delft Center for Systems and Control, Faculty of Mechanical Engineering, Delft University of Technology, The Netherlands<br>
<sup>2</sup> Control Systems Technology Group, Department of Mechanical Engineering, Eindhoven University of Technology, The Netherlands<br>

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
  - Purpose: Provides perceptually uniform colormaps.

---
## Installation

### Automated installation using `INIT.m`:
- Edit SLURM cluster configuration settings in `INIT.m` before running (optional — only needed if using a SLURM cluster)
- Run the `INIT.m` file from the project directory. The output should resemble (using Windows here)
```matlab
MATLAB release : 2024b
Platform       : windows
Download URL   : https://github.com/casadi/casadi/releases/download/3.6.7/casadi-3.6.7-windows64-matlab2018b.zip
Downloading CasADi ...
Extracting to  : <project directory>\bin\casadi-v3.6.7
Done.
CasADi added to MATLAB path.
CasADi 3.6.7 loaded successfully.
Downloading Crameri colour binaries ...
Extracting to  : <project directory>\bin\crameri_colours
Done.
Crameri colour binaries added to MATLAB path.
Crameri colour binaries loaded successfully.
```

### Manual installation:

#### Step 1: Add the src directory to the path
From the project directory in MATLAB, run:
```matlab
addpath(genpath('src'));
```

#### Step 2: SLURM Configuration (optional)
Specify the following parameters in MATLAB and save them to `src/SlurmSettings.mat`:
```matlab
ProfileName = ''; % Enter profile name used for the SLURM cluster (e.g., 'MyClusterProfile')
account     = ''; % Enter account e.g. 'institution-department'
partition   = ''; % Enter partitions to use, comma-separated if multiple (e.g., 'partition1,partition2')
save(fullfile('src','SlurmSettings.mat'),'account','partition','ProfileName');
```
**N.B.**: The above parameters may be left empty if no SLURM cluster is used, but the file must be saved as shown.

#### Step 3: Download CasADi
The code project uses **CasADi v3.6.7**:
1. Download CasADi from [casadi.org](https://web.casadi.org/get/) (version 3.6.7)
2. Extract the contents to: `bin/casadi-v3.6.7/`
3. Add the directory to the path:
   ```matlab
   addpath(fullfile('bin','casadi-v3.6.7'));
   ```
4. Verify CasADi was installed properly by running:
   ```matlab
   x = casadi.SX.sym('x');
   fprintf('CasADi loaded successfully.\n');
   ```

#### Step 4: Install Crameri Colormaps
1. Download [Crameri colormaps](https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)
2. Extract the contents to: `bin/crameri_colours/`
3. Add to path:
   ```matlab
   addpath(fullfile('bin','crameri_colours'));
   ```
4. Verify installation by running the below command (should open an illustration of available colormaps).
   ```matlab
   crameri
   ```
---
## Running the code
### Reproducing results from the article
**Option A: Run Simulations -> Process Data -> Visualize results** <br>
- step 1) Run `main_dRe()`
  - This generates a batch of Monte Carlo simulations for varying Re (innovation noise variance).
  - The parent directory for all of this data should be of the form <br>
    `data/sys1/ref0_prbs/dRe/Re_1e-05_1e-01_15_p_20_N_1e03_f_20_YYYYMMDD_HHmm/`
- step 2) Run `main_processing`
  - This produces a `processed_data.mat` file in the parent directory mentioned directly above.
- step 3) Edit the timestamp of the repository referenced by the `subdir1` variable in  `plot_figs_paper.m`.
- step 4) Run `plot_figs_paper`
  - This will produce all of the plots from the paper.

**Option B: Download (Processed) Data -> Visualize results**
- step 1) Download the already processed data from Zenodo.
- step 2) Run `plot_figs_paper`

### Generating other results
Important parameters that the user can change for different data sets are

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Re` | 0.01 | Innovation noise variance |
| `N` | 1000 | Number of columns in Hankel data matrix |
| `p` | 20 | Past horizon (learning window length) |
| `f` | 20 | Future horizon (prediction/control window) |
| `Ncl` | 1500 | Closed-loop simulation length |
| `seed` | 1 | Random seed for reproducibility |
| `sys` | 1 | System configuration (1 = Landau1995, primary) |
| `ref0` | 'prbs' | Initial reference signal type: 'make*' or 'prbs' |
| `Qk`, `Rk`, `dRk` | 100, 1, 1 | Controller cost weights |

\*'make' specifies an *initial* reference that is based on the reference used in [<a href="#ref25">25</a>]. The output reference trajectory that is employed by the different controllers is always based on such a reference and not a 'prbs' variant.

A single Monte Carlo simulation can be run using the `main` function. Alternatively, batches of Monte Carlo simulations that perform a parameter sweep are performed by the functions
- `main_dRe` - Parameter sweep over Re (innovation noise variance).
- `main_dp` - Parameter sweep over p (past window length).
- `main_dN` - Parameter sweep over N (number of Hankel matrix columns).

Each of the `main_d<X>` functions by default performs 100 Monte Carlo simulations per value of Re, p, or N.<br> For execution on a SLURM cluster:
1. Configure SLURM settings in `INIT.m` (or create `src/SlurmSettings.mat` manually, see installation section above)
2. Running `main_dRe`, `main_dp`, or `main_dN` will automatically:
   - Detect SLURM settings from `SlurmSettings.mat`
   - Submit and manage job runtime via `run_X_ParCluster.m`
   - Save results to an appropriately named data directory
---

## Software Dependencies & Licenses

### **MATLAB (Required)**
- **Version**: R2024b (or compatible)
- **License**: Commercial (proprietary)
- **Toolboxes**: Control System, Robust Control, Statistics & ML, System Identification, Parallel Computing (all v24.2)

### **CasADi (Required)**
- **Version**: v3.6.7
- **License**: LGPL v3.0 (open source, see [CasADi license](https://github.com/casadi/casadi/blob/master/LICENSE.txt))
- **Purpose**: Symbolic differentiation. Used in `src/MC_sim/get_solver.m`

### **Crameri Colormaps (Required for paper figures)**
- **License**: MIT (open source)
- **Source**: [MathWorks File Exchange](https://nl.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)
- **Purpose**: Perceptually uniform scientific colormaps (linear colour gradients among other features).

### **This Code Repository**
- **License**: MIT License (see [LICENSE.md](LICENSE.md))
- **Copyright**: (c) 2026 Data-Driven Control, TU Delft

---

## References
<a id="ref6"></a>[6] Y. Wang, Y. Qiu, M. Sader, D. Huang, and C. Shang, "Data-Driven Predictive Control Using Closed-Loop Data: An Instrumental Variable Approach," *IEEE Control Systems Letters*, vol. 7, pp. 3639–3644, 2023, doi: [10.1109/LCSYS.2023.3340444](https://doi.org/10.1109/LCSYS.2023.3340444).<br>
<a id="ref25"></a>[25] A. Chiuso, M. Fabris, V. Breschi, and S. Formentin, "Harnessing uncertainty for a separation principle in direct data-driven predictive control," *Automatica*, vol. 173, p. 112070, Mar. 2025, doi: [10.1016/j.automatica.2024.112070](https://doi.org/10.1016/j.automatica.2024.112070).