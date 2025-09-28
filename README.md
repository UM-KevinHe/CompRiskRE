# Relative Entropy-Based Competing Risk Models

This software package is for the paper "Cause-Specific Multinomial Relative Entropy-Based Discrete Relative Risk Models for Integrated Prediction of Competing Risk".

## Introduction

Competing risks are frequently encountered in clinical studies, but prediction models face challenges such as low event rates, limited sample sizes, and high dimensionality, leading to unstable parameter estimation and poor performance. Existing survival data integration methods are primarily designed for single-event outcomes, assume homogeneity across data sources, and are built on classical regression frameworks with linear covariate effects. These limitations prevent direct incorporation of heterogeneous or partially overlapping information from historical models, especially when new data provide event-type–specific outcomes but historical models report only aggregated failures.

To overcome these limitations, we propose a multinomial relative entropy–based discrete competing risk integration procedure. This method enables robust integration of published prediction models with new competing risk data, adaptively adjusts the weight of historical information, accommodates partial covariate overlap, handles incomplete historical information, requires only summary-level inputs to protect patient privacy, and is computationally efficient and readily adaptable to various machine learning frameworks.

## Citation

TBD

## Installation

```
#Install the package, need to install the devtools packages:
devtools::install_github("UM-KevinHe/CompRiskRE")
```

## Using CompRiskRE

```
#Use CompRiskRE estimate
eta <- generate_eta(method = "exponential", n = 30, max_eta = 30)
CompRiskRE_FT(
  beta_cor = 0.9,
  eta = eta,
  N_ext = 5000,
  N_loc = 2000,
  N_val = 200,
  N_test = 5000,
  nCause = 2,
  mu = 1,
  sigma = 0.05,
  seed = 2024
)
```
- `beta_cor`: numeric correlation level between external and local models.
- `eta`: numeric vector of Relative Entropy regularization parameters.  
- `N_ext`: sample size for external data (default: 5000).  
- `N_loc`: sample size for local data (default: 2000).  
- `N_val`: sample size for validation data (default: 200).  
- `N_test`: sample size for test data (default: 5000).  
- `nCause`: number of competing causes (default: 2).  
- `mu`: mean of covariates (default: 1).  
- `sigma`: standard deviation of covariates (default: 0.05).  
- `seed`: random seed for reproducibility.  

