# Cassiopea diel pulsation analysis

This repository contains the dataset and R code associated with the manuscript:

**Dawn peaks and individual variation shape functional pulsation rhythms in the photosymbiotic jellyfish Cassiopea***

## Repository content

- `cassiopea_pulsation_dataset.xlsx`  
  Clean dataset used for all analyses.

- `analysis_cassiopea_pulsation.R`  
  R script reproducing the analyses and figures reported in the manuscript.  
  The script also includes additional exploratory/robustness analyses at the end, under the section **ADDITIONAL CODE**.

## Data description

The dataset includes pulsation frequency measurements from 30 Cassiopea spp. individuals monitored over 72 h under semi-natural conditions.

Main variables include:

- individual identity
- experimental day
- time point: Dawn, Noon, Dusk, Midnight
- pulsation frequency
- umbrella diameter
- illuminance
- temperature
- conductivity
- salinity
- pH
- dissolved oxygen

## Analyses

The R script includes:

- data import and cleaning
- descriptive plots of pulsation frequency over 72 h
- collinearity assessment among environmental variables
- linear mixed-effects models
- post hoc comparisons using estimated marginal means
- cyclic GAMM analysis of 24 h temporal structure
- Lomb–Scargle periodogram
- variance decomposition and repeatability analysis
- additional exploratory analyses

## Requirements

Analyses were performed in R version 4.5.1.

Required R packages:

readxl
dplyr
tidyr
ggplot2
lme4
lmerTest
emmeans
car
rptR
janitor
lubridate
lomb
mgcv
MASS
