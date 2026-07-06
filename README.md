# Lupine IPM Uncertainty Analysis

This repository contains scripts for analyzing uncertainty in integral projection models of *Lupinus tidestromii* from Point Reyes National Seashore, collected from 2005 to 2018. Models are adapted from [Compagnoni et al. 2021](https://doi.org/10.1002/ecs2.3454).

### Data Availability

Due to restrictions associated with sharing the underlying raw demographic dataset, the complete raw data used in this study are not publicly distributed as part of this repository. However, this repository contains the full analytical workflow used to process the data, fit vital rate models, construct integral projection models, bootstrap resample the parameter values, conduct the full uncertainty analysis, and generate manuscript figures.

To maximize transparency and reproducibility, several upstream scripts that rely on the raw dataset are included for reference. These scripts document the preprocessing, model fitting, and model selection workflow used in the case study, including generation of exploratory figures used to evaluate model fit. Although these scripts cannot be executed without access to the raw data, they are provided to allow readers and reviewers to fully inspect the analytical pipeline.

All downstream analyses beginning with the mean parameter sets are fully reproducible using the files provided in this repository together with the datasets hosted on Zenodo, which also include an example dataframe illustrating the structure of the raw demographic data.

Zenodo DOI: [DOI placeholder]

A workflow overview is provided below, with fully-reproducible scripts indicated.


### Software Environment

Analyses were conducted in R version 4.3.0 "Already Tomorrow"

Package dependencies include:

| Package | Version |
|---|---|
| bbmle | 1.0.25 |
| binom | 1.1.1.1 |
| dplyr | 1.1.2 |
| extraDistr | 1.10.0 |
| ggplot2 | 3.4.4 |
| ggthemes | 4.2.4 |
| goodpractice | 1.0.5 |
| Hmisc | 5.0.1 |
| ipmr | 0.0.7 |
| lintr | 3.1.2 |
| lme4 | 1.1.33 |
| mgcv | 1.8.42 |
| patchwork | 1.3.0 |
| readxl | 1.4.2 |
| rgdal | 1.6.6 |
| sensemakr | 0.1.6 |
| testthat | 3.1.7 |
| tidyr | 1.3.0 |
| tidyverse | 2.0.0 |


### Workflow
- `01_VR_formatting.R`: Construct functions to subset dataframes from the entire dataset to fit vital rate models on, and visualize the raw data
- `02_VR_models.R`: Fit vital rate models, plot models over the raw data, perform model selection, and export model parameters
- `03_IPMs.R`: Construct IPMs from mean model parameter values (fully reproducible)
- `04_Sampling.R`: Bootstrap resample demographic data and refit vital rate models to generate a dataframe of sampled parameter values
- `05_Uncertainty.R`: Perform the uncertainty analysis (fully reproducible)
- `06_Simulation.R`: Simulate germination data and analyze contribution to uncertainty (fully reproducible)
- `07_Figures.R`: Prepare figures for publication (fully reproducible)
