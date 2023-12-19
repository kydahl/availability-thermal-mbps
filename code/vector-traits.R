## Title: Mosquito Thermal Trait Data processing ###############################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Conduct all parts of the analysis
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Load empirical trait data
##           2) Fit trait thermal performance curves to trait data
##           3) Translate traits into model parameters
##           4) Build data set incorporating all axes of variation
##           5) Calculate model outputs
##           6) Illustrate model outputs (might be done separately)
##           7) Conduct sensitivity analysis
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023


# 0) Set-up, load in necessary packages and data-sets ---------------------
library(tidyverse)

# 1) Load empirical trait data --------------------------------------------

# # Run this to produce workable dataset from raw data for the first time
plot_bool <- FALSE # decide whether you'd like to generate a diagnostic plot
source("code/data-cleaning.R")

# # Alternatively, run this to load pre-processed dataset
# # load dataset for fitting trait TPCs
# data.in.TPC <- read_rds("data/clean/data_for_TPC_fitting.rds")

# 2) Fit trait thermal performance curves to trait data -------------------

# Run this to generate samples of trait TPC parameters from informed posterior distributions for the first time
# Set parameters for MCMC
n.chains <- 2 #5
n.adapt <- 1000 #5000
n.samps <- 1000 #5000
# Do you want to look at diagnostic plots?
plot_bool <- FALSE
source("code/get-thermal-trait-priors.R")
write_rds(data.in.transform, "results/TPC_param_samples.rds")

# Get summary statistics of each trait for each system
table <- read_rds("results/TPC_param_samples.rds") %>% 
  pivot_longer(cols = c("c","T0", "Tm"), names_to = "hyperparameter") %>% 
  group_by(system_name = system_ID, trait, hyperparameter, `function` = func) %>% 
  filter(!is.na(value)) %>% # removes T0 for Culex lifespan, which is linear and hence only has 2 parameters
  mutate(mean = mean(value),
            median = median(value),
            lowHCI_89 = quantile(value, 0.055),
            highHCI_89 = quantile(value, 0.945)) %>% 
  dplyr::select(-c(value, sample_num, func, system_ID)) %>% 
  distinct()

write_csv(table, "results/TPC_fit_summary.csv")

# # Alternatively, run this to load pre-processed data set
# data.in.transform <- read_rds("results/TPC_param_samples.rds")

# 3) Translate traits into model parameters -------------------------------

# Define temperature range of study
Temps <- seq(10, 40, length.out = 601) # full: length.out = 601, thin: length.out = 301

# Thin sample size
thin_size <- 1000 # full = 5000, thin = 1000

plot_bool = FALSE
source("code/trait-transform.R")

write_rds(data.in.params, "results/VecTPC_vals.rds", compress = "gz")

# Save trait mean values
mean.Vec <- data.in.params %>%
  dplyr::select(-c(mosquito_species, pathogen, muL, etaL)) %>%
  pivot_longer(cols = lf:V0, names_to = "variable", values_to = "value") %>%
  group_by(system_ID, Temperature, variable) %>%
  summarise(mean = mean(value), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(names_from = "variable", values_from = "mean") %>%
  mutate(etaL = rhoL / (deltaL + eps)) %>%
  mutate(muL = etaL - rhoL) %>%
  mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                     0,
                     KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL))
  ))

write_rds(mean.Vec, "results/meanVec_vals.rds", compress = "gz")

# Save trait median values (this one is used more frequently)
median.Vec <- data.in.params %>%
  dplyr::select(-c(mosquito_species, pathogen, muL, etaL)) %>%
  pivot_longer(cols = lf:V0, names_to = "variable", values_to = "value") %>%
  group_by(system_ID, Temperature, variable) %>%
  summarise(median = median(value), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(names_from = "variable", values_from = "median") %>%
  mutate(etaL = rhoL / (deltaL + eps)) %>%
  mutate(muL = etaL - rhoL) %>%
  mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                     0,
                     KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL))
  ))
write_rds(median.Vec, "results/medianVec_vals.rds", compress = "gz")

# remove work sets (these may be taking up a lot of memory)
rm("combined_df", "Infection_df", "noInfection_df", "TPC_df", "data.in.transform")
