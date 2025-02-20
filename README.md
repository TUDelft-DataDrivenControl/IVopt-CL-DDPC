# Data-driven predictive control with optimal IVs
This repository contains the necessary components of a data-driven predictive control implementation that employs optimal IVs

## Project organization
- PG = project-generated
- HW = human-writable
- RO = read only
```
.
├── .gitignore
├── CITATION.md
├── LICENSE.md
├── README.md
├── bin                         <- Compiled and external code, ignored by git (PG)
│   ├── external/Casadi         <- When running simulations on a cluster, place compiled CasADi code with IPOPT solver here. (RO)
|   └── model_Favoreel1999.m    
├── data                        <- All project data, ignored by git
│   ├── processed               <- The final, canonical data sets for modeling. (PG)
│   ├── raw                     <- The original, immutable data dump. (RO)
│   └── temp                    <- Intermediate data that has been transformed. (PG)
├── results         
│   └── figures                 <- Figures for the manuscript or reports (PG)
└── src                         <- Source code for this project (HW)

```

## License

This project is licensed under the terms of the [MIT License](/LICENSE.md)
