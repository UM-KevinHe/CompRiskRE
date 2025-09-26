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
install.packages("CompRiskRE")
devtools::install_github("UM-KevinHe/CompRiskRE")
```

## Using CompRiskRE

```
#Use CompRiskRE estimate

```

  
## Simulation example

