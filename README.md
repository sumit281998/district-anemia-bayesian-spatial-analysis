# Bayesian Spatial Analysis of Childhood Anaemia in India

This repository contains the R code used for the district-level
Bayesian spatial analysis of childhood anaemia using NFHS-5
(2019–21) data.

## Main analysis

The analysis includes:

- District-level descriptive analysis
- Global Moran's I
- Local Indicators of Spatial Association (LISA)
- Beta regression
- Bayesian spatial beta regression
- BYM2 spatial random effects
- DIC and WAIC model comparison
- Fixed-effect posterior estimates
- District-level posterior spatial effects
- Spatial visualisation

## Statistical model

The primary spatial model uses a BYM2 specification based on
queen-contiguity between districts.

The mosquito-net exposure is specified as a three-level categorical
variable:

1. No net
2. Untreated mosquito net
3. Insecticide-treated mosquito net (ITN)

The reference category is no net.

## Software

- R
- INLA
- sf
- spdep
- dplyr
- ggplot2
- tmap

See the R script for the complete analysis workflow and
session information.

## Data availability

The underlying NFHS-5 data are not redistributed in this repository.
Users should obtain the data through the appropriate official
data-access procedures and construct the analysis dataset according
to the variable definitions documented in the manuscript.

## Reproducibility

The analysis script is designed to run after the required analysis
dataset and district boundary data have been supplied.
