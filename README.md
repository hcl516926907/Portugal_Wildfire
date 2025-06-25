# XGBoost meets INLA: a two-stage spatio-temporal forecasting of wildfires in Portugal

This repository contains the code to reproduce the results presented in the paper:

**"XGBoost meets INLA: a two-stage spatio-temporal forecasting of wildfires in Portugal"**

We propose a two-stage probabilistic forecasting framework for predicting monthly wildfire activity at the municipality (council) level across Portugal. The goal is to forecast both the number of fire events and the total burnt area, incorporating environmental drivers and spatio-temporal dependencies.


## Overview of the Methodology

1. **Stage 1 – Deterministic Forecasting with XGBoost**  
   An XGBoost model is trained to produce **point forecasts** of fire counts and burnt area using meteorological and land surface covariates. The model is trained at the council-month level.

2. **Stage 2 – Probabilistic Modelling with INLA**  
   The outputs from XGBoost are used as covariates in a **latent Gaussian model** fitted using **INLA (Integrated Nested Laplace Approximation)**. This stage also incorporates **spatial adjacency information** between councils and **temporal encoding** to generate **full predictive distributions**, allowing for uncertainty quantification in the forecasts.


## Repository Structure

### `Data_Prep/`

Scripts for data preprocessing, cleaning, and covariate integration.

- **`Covariates_Preparation.R`**  
  - Processes environmental covariates (e.g., temperature, precipitation) from **ERA5 reanalysis data**.  
  - Merges these covariates with wildfire occurrence records and aggregates all data to the **council-month** level.

- **`Data_Quality_Check.R`**  
  - Assesses the quality and completeness of wildfire records from **ICNF (Instituto da Conservação da Natureza e das Florestas)**.  
  - Includes checks for duplicates, missing entries, and spatial inconsistencies.

- **`Source_Functions.R`**  
  - Supporting functions used by `Data_Quality_Check.R`.  
  - Adapted from [UrbanFiresData by rb1970](https://github.com/rb1970/UrbanFiresData).

---

### Files

- **`INLA_Council_h1.R`**  
  - Fits the second-stage **hurdle model** using INLA, modelling both the fire count and burnt area of wildfires.  
  - Incorporates spatial and temporal random effects, as well as XGBoost-predicted covariates.

- **`Model_Diagnostics.R`**  
  - Performs model comparison and validation.  
  - Generates diagnostic plots and scoring statistics to evaluate the predictive performance of different model variants.


- **`Posterior_Prediction.R`**  
  - Performs posterior sampling of latent effects and hyperparameters from the fitted INLA model.  
  - Generates probabilistic wildfire forecasts at the council-month level.

- **`XGBoost_Bayesian_Optimization_Tuning.R`**  
  - Tune hyperparameters for the XGBoost model using **Bayesian optimisation**.  
  - Selects optimal parameter values based on cross-validation performance.  
  - Produces one-month-ahead point forecasts of fire count and burnt area.
