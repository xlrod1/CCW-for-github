# Clone–Censor–Weight (CCW) Method for Target Trial Emulation

This repository contains R code used in the **target trial emulation study** investigating the causal effect of **induction timing after prelabor rupture of membranes (PROM)** on delivery outcomes, implemented using the **Clone–Censor–Weight (CCW)** method.  
It accompanies the analyses presented in *Axelrod et al., “Trial Emulation to Compare Induction Timing Strategies After Prelabor Rupture of Membranes”* (manuscript in preparation).

---

## 📂 Repository Structure
code/
└── CCW_method/
├── Target_Trial_1/
│ ├── balance_over_time/
│ ├── boot_SW_CIF_TV.R
│ ├── f_compute_survival_for_clone_TV.R
│ └── main_code.R
│
├── Target_Trial_2/
│ ├── boot_SW_CIF_TV.R
│ ├── cloning_functions.R
│ ├── f_compute_survival_for_clone_TV.R
│ └── main_code.R
│
└── Non_causal_analysis/
├── melamed.R
├── melamed_CIF.R
└── melamed_CIF_diff.R

---

## 🔍 Main Components

### Target_Trial_1/
Contains the full implementation of **Target Trial 1**, including:
- `main_code.R`: complete workflow for cloning, censoring, weighting, and survival estimation.  
- `f_compute_survival_for_clone_TV.R`: function for weighting and estimating the survival curves and CIFs curves time‐varying covariates.  
- `boot_SW_CIF_TV.R`: bootstrap procedure for estimating uncertainty in weighted CIFs.  
- `balance_over_time/`: scripts for assessing covariate balance over time (unweighted vs. weighted).

---

### Target_Trial_2/
Implementation of a second trial scenario using similar CCW steps but different treatment strategies.
- `main_code.R`  
- `f_compute_survival_for_clone_TV.R`  
- `cloning_functions.R`  
- `boot_SW_CIF_TV.R`

---

### Non_causal_analysis/
Classical (non-causal) analyses for comparison:
- `melamed.R`: replication of the classical Melamed et al. (2023) analysis.  
- `melamed_CIF.R`: computes CIFs using standard survival methods.  
- `melamed_CIF_diff.R`: computes contrasts in CIFs between groups.

---


📊 Outputs

The scripts produce:

Weighted Kaplan–Meier survival curves and CIFs for delivery outcomes.

Standardized mean difference (SMD) plots over time (for assessing positivity).

Bootstrap confidence intervals for CIF contrasts across treatment regimes.
