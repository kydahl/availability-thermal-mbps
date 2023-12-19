## Title: Sensitivity and uncertainty analysis of Topt, CTmin, CTmax, and CTwidth
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Determine which traits most underlie the variation in key outputs
##
## Contents: 0) Set up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Set up data frames
##           3) Topt local, global sensitivity, uncertainty
##           4) CTmin local, global sensitivity, uncertainty
##           5) CTmax local, global sensitivity, uncertainty
##           6) CTwidth local, global sensitivity, uncertainty
##           *) Diagnostics and visualizations
##
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023
# _______________________________________________________________________________

# 0) Set-up, load in necessary packages and data-sets ---------------------

# Load Libraries
library(tidyverse)
library(reshape2)
library(multidplyr)
library(foreach)
library(doParallel)
library(progress)
library(latex2exp) # nec
library(viridis) # nec
library(cowplot) # nec
library(MetBrewer) #nec
library(svglite) #nec
library(HDInterval)

# Set up parallel
if (!exists("cl")) {
  cl <- parallel::makeCluster(parallel::detectCores() - 2)
  registerDoSNOW(cl)
}

eps <- .Machine$double.eps

# Initialize a data frame to save our analyses
init.df <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                  sigmaH = c(), KH = c(), variable = c(),
                  lowHCI = c(), highHCI = c(), mean = c(), median = c())

# Load needed datasets

# Host trait dataset
data.Host <- read_rds("results/Host_vals.rds")
# Vector trait dataset
data.Vec <- read_rds("results/VecTPC_vals.rds")

# R0 data
data.R0 <- read_rds("results/R0_vals.rds")

# 1) Define accessory functions -------------------------------------------

# Function: Slice data to optimize usage of parallel processing and memory
slice_func<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

# Label functions: 
appender_sigmaH <- function(string) {
  unname(TeX(paste("$\\sigma_H = $", string)))}
appender_KH <- function(string) {
  unname(TeX(paste("$K_H = $", string)))}

# 2) Set up mosquito trait data frames ----------------------------------------
## Get summary statistics of vector traits across samples
# Mean values
mean.Vec <- read_rds("results/meanVec_vals.rds")

# Median values (this one is used more frequently)
median.Vec <- read_rds("results/medianVec_vals.rds")

# 3) R0 uncertainty -------------------------------------------------------

# Calculate the width of the 95% HPD for the full posterior of R0 across temperatures

sigmaH_vec <- c(10, 100, Inf)
KH_vec <- 10^seq(-2,4)

# Set up host trait data set
data.R0HPD <- dplyr::filter(data.Host,
                            sigmaH %in% sigmaH_vec,
                            KH %in% KH_vec)

# Calculate R0 HPD width with full variation
full.R0.HPD <- expand_grid(data.R0HPD, 
                           data.Vec) %>% 
  data.table::data.table() %>%
  mutate(RV = ifelse(is.infinite(sigmaH),
                     sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                     sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
  mutate(bV = ifelse(is.infinite(sigmaH),
                     sigmaV, # Ross-Macdonald model
                     sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
  mutate(RH = ifelse(V0 == 0,
                     0,
                     bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(RV*RH)) %>%
  select(system_ID, sample_num, sigmaH, KH, Temperature, R0) %>% 
  group_by(system_ID, sigmaH, KH, Temperature) %>%
  summarise(
    HPD_low = hdi(R0, credMass = 0.95)[1],
    HPD_high = hdi(R0, credMass = 0.95)[2],
    HPD_width = max(eps,HPD_high-HPD_low),
    .groups = "keep"
  ) %>% 
  group_by(system_ID, sigmaH, KH) %>%
  mutate(norm_HPD_width = HPD_width / max(HPD_width))

# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0, lf)) %>% 
  colnames() %>% c("muV")

# Initialize data frame
R0.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                 sigmaH = c(), KH = c(), 
                 HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c())

# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    mutate(muV = 1/lf) %>% 
    # Remove non-focal mosquito traits
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    # Add back in the median values of non-focal mosquito traits
    full_join(mean.Vec %>% 
                mutate(muV = 1/lf) %>% 
                select(-all_of(var_name))%>% 
                distinct()) %>% 
    # Re-calculate lifespan and adult mosquito density
    mutate(lf = 1/muV) %>%
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>% 
    mutate(
      HPD_width_var = hdi(!! sym(var_name), credMass = 0.95)[2]-hdi(!! sym(var_name), credMass = 0.95)[1]
    )
  
  # b) Get posterior samples of R0 (as a function of temperature)
  temp_df <- expand_grid(data.R0HPD, data.HPD.Vec) %>%
    mutate(lf = 1/muV) %>% 
    data.table::data.table() %>%
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
    mutate(RH = ifelse(V0 == 0,
                       0,
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    select(system_ID, sample_num, sigmaH, KH, Temperature, R0, HPD_width_var) %>% 
    # c) Calculate the 95% HPD at each temperature
    group_by(system_ID, sigmaH, KH, Temperature) %>% 
    mutate(
      HPD_low = hdi(R0, credMass = 0.95)[1],
      HPD_high = hdi(R0, credMass = 0.95)[2],
      HPD_width = HPD_high-HPD_low,
      .groups = "keep"
    ) %>% 
    select(system_ID, sigmaH, KH, Temperature, HPD_width, HPD_width_var) %>% 
    right_join(full.R0.HPD %>% 
                 select(-c(HPD_low, HPD_high)) %>% 
                 rename(full_HPD_width = HPD_width)) %>% 
    mutate(focal_var = var_name) %>% 
    group_by(system_ID, sigmaH, KH, focal_var, Temperature) %>% 
    # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
    mutate(rel_HPD_width = ifelse(full_HPD_width %in% c(0,eps), 0, HPD_width / full_HPD_width)) %>%
    mutate(rel_HPD_width_var = HPD_width / HPD_width_var) %>% 
    distinct()
  
  R0.HPD <- rbind(R0.HPD, temp_df)
}

# Save R0 relative highest posterior density data
write_rds(R0.HPD, "results/R0_HPD_unc.rds")

# 4) ddT R0 uncertainty -------------------------------------------------------

# Calculate the width of the 95% HPD for the full posterior of ddT R0 across temperatures

# Set up host trait data frame
sigmaH_vec <- c(100, Inf)
KH_vec <- 10^seq(-2,4)

data.ddTR0HPD <- dplyr::filter(data.Host,
                               sigmaH %in% sigmaH_vec,
                               KH %in% KH_vec)
# Initialize full HPD dataframe
full.ddTR0.HPD <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                         sigmaH = c(), KH = c(), 
                         HPD_low = c(), HPD_high = c(), HPD_width = c())
# Calculate full ddTR0 HPD dataframe
full.ddTR0.HPD <- expand_grid(data.ddTR0HPD, 
                              data.Vec) %>% 
  data.table::data.table() %>%
  mutate(RV = ifelse(is.infinite(sigmaH),
                     sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                     sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
  mutate(bV = ifelse(is.infinite(sigmaH),
                     sigmaV, # Ross-Macdonald model
                     sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
  mutate(RH = ifelse(V0 == 0,
                     0,
                     bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(RV*RH)) %>%
  ungroup() %>% 
  arrange(system_ID, sample_num, sigmaH, KH, Temperature) %>% 
  # Calculate temperature derivative of R0
  mutate(ddTR0 = (R0 - lag(R0)) / (Temperature - lag(Temperature))) %>% 
  select(system_ID, sample_num, sigmaH, KH, Temperature, ddTR0) %>% 
  group_by(system_ID, sigmaH, KH, Temperature) %>% 
  summarise(
    HPD_low = hdi(ddTR0, credMass = 0.95)[1],
    HPD_high = hdi(ddTR0, credMass = 0.95)[2],
    HPD_width = max(eps,HPD_high-HPD_low),
    .groups = "keep"
  )

# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0, lf)) %>% 
  colnames() %>% c("muV")

# Initialize data frame
ddTR0.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                    sigmaH = c(), KH = c(), 
                    HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c())

# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    mutate(muV = 1/lf) %>% 
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    full_join(mean.Vec %>% 
                mutate(muV = 1/lf) %>% 
                select(-all_of(var_name))%>% 
                distinct()) %>% 
    mutate(lf = 1/muV) %>%
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>%
    select(-lf)
  # b) Get posterior samples of R0 (as a function of temperature)
  temp_df <- expand_grid(data.ddTR0HPD, data.HPD.Vec) %>%
    mutate(lf = 1/muV) %>% 
    data.table::data.table() %>%
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
    mutate(RH = ifelse(V0 == 0,
                       0,
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    ungroup() %>% 
    arrange(system_ID, sample_num, sigmaH, KH, Temperature) %>% 
    # Calculate temperature derivative of R0
    mutate(ddTR0 = (R0 - lag(R0)) / (Temperature - lag(Temperature))) %>% 
    select(system_ID, sample_num, sigmaH, KH, Temperature, ddTR0) %>% 
    # c) Calculate the 95% HPD at each temperature
    group_by(system_ID, sigmaH, KH, Temperature) %>% 
    summarise(
      HPD_low = hdi(ddTR0, credMass = 0.95)[1],
      HPD_high = hdi(ddTR0, credMass = 0.95)[2],
      HPD_width = HPD_high-HPD_low,
      .groups = "keep"
    ) %>% 
    select(system_ID, sigmaH, KH, Temperature, HPD_width) %>% 
    right_join(full.ddTR0.HPD %>% 
                 select(-c(HPD_low, HPD_high)) %>% 
                 rename(full_HPD_width = HPD_width)) %>% 
    mutate(focal_var = var_name) %>% 
    group_by(system_ID, sigmaH, KH, focal_var, Temperature) %>% 
    # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
    mutate(rel_HPD_width = ifelse(full_HPD_width %in% c(0,eps), 0, HPD_width / full_HPD_width))
  
  ddTR0.HPD <- rbind(ddTR0.HPD, temp_df)
}

# Save R0 relative highest posterior density data
write_rds(ddTR0.HPD, "results/ddTR0_HPD_unc.rds")

# 5) Topt uncertainty -------------------------------------------------------

# Calculate the width of the 95% HPD for the full posterior of Topt across vertebrate host abundance

# Set up host trait dataframe
KH_vec <- 10^seq(-1.04, 3.02, by = .02)
KH_length = length(KH_vec)
sigmaH_vec <- 100
data.ToptHPD <- dplyr::filter(data.Host, sigmaH %in% sigmaH_vec) %>% 
  select(-KH) %>% 
  cross_join(as_tibble(list(KH = KH_vec))) %>% 
  distinct()

# Initialize full Topt HPD dataframe
full.Topt.HPD <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                        sigmaH = c(), KH = c(), 
                        HPD_low = c(), HPD_high = c(), HPD_width = c())

# Set up progress bar
pb <- txtProgressBar(max = KH_length, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Calculate full Topt HPD widths
full.Topt.HPD <- foreach(index_KH = 1:KH_length,
                         .packages = c("tidyverse", "HDInterval"),
                         .options.snow = opts,
                         .combine = 'rbind') %dopar%  {
                           expand_grid(dplyr::filter(data.ToptHPD, 
                                                     KH %in% KH_vec[index_KH]), 
                                       data.Vec) %>%  
                             ungroup() %>% 
                             mutate(RV = ifelse(is.infinite(sigmaH),
                                                sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                                                sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
                             mutate(bV = ifelse(is.infinite(sigmaH),
                                                sigmaV, # Ross-Macdonald model
                                                sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
                             mutate(RH = ifelse(V0 == 0,
                                                0,
                                                bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
                             # Basic reproduction number
                             mutate(R0 = sqrt(RV*RH))  %>%
                             # filter to maximum value of R0
                             dplyr::filter(R0>0) %>% 
                             group_by(system_ID, sample_num, sigmaH, KH) %>%
                             dplyr::filter(R0 == max(R0)) %>%
                             distinct() %>% 
                             # Get temperature at which R0 is maximized
                             rename(Topt = Temperature) %>% 
                             select(-R0) %>%
                             ungroup() %>% 
                             # select(system_ID, sample_num, sigmaH, Topt) %>%
                             group_by(system_ID, sigmaH) %>%
                             mutate(HPD_width = max(eps, hdi(Topt, credMass = 0.95)[2] - hdi(Topt, credMass = 0.95)[1])) %>%
                             mutate(KH = KH_vec[index_KH])
                         } %>% 
  select(system_ID, sigmaH, KH, HPD_width) %>% 
  distinct()
# Close progress bar
close(pb)

# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0,)) %>% 
  colnames()

# Initialize data frame
Topt.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                   sigmaH = c(), KH = c(), 
                   HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c(), rel_HPD_width_var = c())

# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  print(paste0("Focal variable: ", var_name))
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    full_join(mean.Vec %>% 
                select(-all_of(var_name))%>% 
                distinct(),
              by = join_by(system_ID, Temperature, KL)) %>% 
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL)))) %>% 
    ungroup()
  
  pb <- txtProgressBar(max = KH_length, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # b) Get posterior samples of R0 (as a function of temperature)
  temp_df <- foreach(index_KH = 1:KH_length,
                     .packages = c("tidyverse", "HDInterval"),
                     .options.snow = opts,
                     .combine = 'rbind') %dopar%  {
                       expand_grid(dplyr::filter(data.ToptHPD, 
                                                 KH %in% KH_vec[index_KH]), 
                                   data.HPD.Vec) %>%  
                         ungroup() %>% 
                         mutate(RV = ifelse(is.infinite(sigmaH),
                                            sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                                            sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
                         mutate(bV = ifelse(is.infinite(sigmaH),
                                            sigmaV, # Ross-Macdonald model
                                            sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
                         mutate(RH = ifelse(V0 == 0,
                                            0,
                                            bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
                         # Basic reproduction number
                         mutate(R0 = sqrt(RV*RH))  %>%
                         # dplyr::filter to maximum value of R0
                         dplyr::filter(R0>0) %>% 
                         group_by(system_ID, sample_num, sigmaH, KH) %>%
                         dplyr::filter(R0 == max(R0)) %>%
                         distinct() %>% 
                         # Get temperature at which R0 is maximized
                         rename(Topt = Temperature) %>% 
                         select(-R0) %>%
                         ungroup() %>% 
                         group_by(system_ID, sigmaH) %>%
                         mutate(HPD_width = max(eps, hdi(Topt, credMass = 0.95)[2] - hdi(Topt, credMass = 0.95)[1])) %>%
                         mutate(
                           HPD_width_var = hdi(!! sym(var_name), credMass = 0.95)[2] - hdi(!! sym(var_name), credMass = 0.95)[1]
                         ) %>% 
                         ungroup() %>% 
                         mutate(KH = KH_vec[index_KH])
                     } %>% 
    select(system_ID, sigmaH, KH, HPD_width, HPD_width_var) %>% 
    distinct() %>% 
    right_join(full.Topt.HPD %>% 
                 rename(full_HPD_width = HPD_width),
               by = join_by(system_ID, sigmaH, KH)) %>% 
    distinct() %>% 
    mutate(focal_var = var_name) %>% 
    group_by(system_ID, sigmaH, KH, focal_var) %>% 
    # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
    mutate(rel_HPD_width = ifelse(full_HPD_width < 1.1 * eps, 0, HPD_width / full_HPD_width)) %>% 
    mutate(rel_HPD_width_var = HPD_width / HPD_width_var)
  close(pb)
  
  Topt.HPD <- rbind(Topt.HPD, temp_df)
}

# Save Topt relative highest posterior density data
write_rds(Topt.HPD, "results/Topt_HPD_unc.rds")


# 6) dTopt dKH uncertainty -------------------------------------------------------
# Calculate the width of the 95% HPD for the full posterior of dTopt/dKH across vertebrate host abundance

# Set up host trait data frame
KH_vec <- 10^seq(-1.04, 3.02, by = .02)
KH_length = length(KH_vec)
sigmaH_vec <- 100
data.dToptdKH.HPD <- dplyr::filter(data.Host, sigmaH %in% sigmaH_vec) %>% 
  select(-KH) %>% 
  full_join(as_tibble(list(KH = KH_vec)), by = character()) %>% 
  distinct()

# Initialize full dToptdKH HPD width dataframe
full.dToptdKH.HPD <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                            sigmaH = c(), KH = c(), 
                            HPD_low = c(), HPD_high = c(), HPD_width = c())

# Set up progress bar
iterations <- KH_length-1
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Function: calculate dToptdKH HPD width
dToptdKH.HPD_func <- function(Vector_data, index_KH) {
  expand_grid(dplyr::filter(data.dToptdKH.HPD, 
                            KH %in% KH_vec[(index_KH-1):index_KH]), 
              Vector_data) %>%  
    ungroup() %>% 
    mutate(RV = sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0))) %>%
    mutate(bV = sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps)) %>%
    mutate(RH = ifelse(V0 == 0,
                       0,
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    # Filter to maximum value of R0
    dplyr::filter(R0>0) %>% 
    select(system_ID, sample_num, sigmaH, KH, Temperature, R0) %>% 
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    dplyr::filter(R0 == max(R0)) %>%
    distinct() %>% 
    # Get temperature at which R0 is maximized
    rename(Topt = Temperature) %>% 
    select(-R0) %>%
    ungroup() %>% 
    arrange(system_ID, sample_num, KH) %>%
    # Calculate numerical derivative of Topt wrt KH
    mutate(dToptdKH = (Topt - lag(Topt)) / (KH_vec[index_KH] - KH_vec[index_KH-1])) %>%
    ungroup() %>%
    select(system_ID, sample_num, dToptdKH) %>%
    group_by(system_ID) %>%
    summarise(
      HPD_low = hdi(dToptdKH, credMass = 0.95)[1],
      HPD_high = hdi(dToptdKH, credMass = 0.95)[2],
      HPD_width = max(eps, HPD_high-HPD_low),
      .groups = "keep"
    ) %>%
    mutate(sigmaH = sigmaH_vec,
           KH = KH_vec[index_KH])
}

# Compute full dToptdKH HPD width
full.dToptdKH.HPD <- foreach(index_KH = 2:KH_length,
                             .packages = c("tidyverse", "HDInterval"),
                             .options.snow = opts,
                             .combine = 'rbind') %dopar%  {
                               dToptdKH.HPD_func(data.Vec, index_KH)
                             }

# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0, lf)) %>% 
  colnames() %>% c("muV")

# Initialize data frame
dToptdKH.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                       sigmaH = c(), KH = c(), 
                       HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c())

# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  print(paste0("Focal variable: ", var_name))
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    mutate(muV = 1/lf) %>% 
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    full_join(mean.Vec %>% 
                mutate(muV = 1/lf) %>% 
                select(-all_of(var_name))%>% 
                distinct()) %>% 
    mutate(lf = 1/muV) %>%
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL))))
  
  # Set up progress bar
  iterations <- KH_length-1
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  temp_df <- foreach(index_KH = 2:KH_length,
                     .packages = c("tidyverse", "HDInterval"),
                     .options.snow = opts,
                     .combine = 'rbind') %dopar%  {
                       dToptdKH.HPD_func(data.HPD.Vec, index_KH) %>%
                         right_join(full.dToptdKH.HPD %>% 
                                      select(-c(HPD_low, HPD_high)) %>% 
                                      rename(full_HPD_width = HPD_width) %>% 
                                      filter(KH == KH_vec[index_KH]),
                                    by = c("system_ID", "sigmaH", "KH")) %>% 
                         mutate(focal_var = var_name) %>% 
                         group_by(system_ID, sigmaH, KH, focal_var) %>% 
                         # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
                         mutate(rel_HPD_width = ifelse(full_HPD_width < 1.1 * eps, 0, HPD_width / full_HPD_width))
                     }
  dToptdKH.HPD <- rbind(temp_df, dToptdKH.HPD)
  close(pb)
}

# Save Topt relative highest posterior density data
write_rds(dToptdKH.HPD, "results/dToptdKH_HPD_unc.rds")


# 7) CTmin/max/width uncertainty ------------------------------------------
# Calculate the width of the 95% HPD for the full posterior of Topt across vertebrate host abundance

# Set up host trait data frame
sigmaH_vec <- c(10, 100, Inf)
KH_vec <- 10^seq(-2, 4, length.out = 601)
data.CTHPD <- dplyr::filter(data.Host, sigmaH %in% sigmaH_vec) %>% 
  select(-KH) %>% 
  full_join(as_tibble(list(KH = KH_vec)), by = character()) %>% 
  distinct()

# Initialize full CT HPD width dataframe
full.CT.HPD <- tibble(system_ID = c(), Temperature = c(), Model = c(),
                      sigmaH = c(), KH = c(), variable = c(),
                      HPD_low = c(), HPD_high = c(), HPD_width = c())

# Calculate full CT HPD width
for (index_KH in unique(data.CTHPD$KH)) {
  full.CT.HPD <- expand_grid(dplyr::filter(data.CTHPD, KH == index_KH), 
                             data.Vec) %>%
    data.table::data.table() %>%
    mutate(RV = ifelse(is.infinite(sigmaH),
                       sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                       sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
    mutate(bV = ifelse(is.infinite(sigmaH),
                       sigmaV, # Ross-Macdonald model
                       sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
    mutate(RH = ifelse(V0 == 0,
                       0,
                       bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
    # Basic reproduction number
    mutate(R0 = sqrt(RV*RH)) %>%
    # Filter to maximum value of R0
    group_by(system_ID, sample_num, sigmaH, KH) %>%
    dplyr::filter(R0 > 1) %>%
    # Get lowest temperature at which R0 exceeds one
    mutate(CTmin = min(Temperature)) %>%
    # Get highest temperature at which R0 exceeds one
    mutate(CTmax = max(Temperature)) %>%
    # Get width of critical thermal interval
    mutate(CTwidth = CTmax - CTmin) %>% 
    # Add back values removed from filtering R0>1 above
    # and assign the right NA values
    full_join(expand_grid(
      dplyr::filter(data.CTHPD, KH == index_KH), data.Vec) %>% 
                select(c(Model, system_ID, sample_num, sigmaH, KH)) %>% 
                distinct(),
      by = c("sigmaH", "Model", "KH", "system_ID", "sample_num")
      ) %>% 
    replace_na(list(CTmin = Inf,
                    CTmax = -Inf,
                    CTwidth = 0)) %>% 
    pivot_longer(cols = c(CTmin, CTmax, CTwidth), names_to = "variable", values_to = "value") %>% 
    select(system_ID, sample_num, sigmaH, KH, variable, value) %>%
    group_by(system_ID, sigmaH, KH, variable) %>% 
    summarise(
      HPD_low = hdi(value, credMass = 0.95)[1],
      HPD_high = hdi(value, credMass = 0.95)[2],
      HPD_width = max(eps, HPD_high-HPD_low),
      .groups = "keep"
    ) %>% 
    rbind(full.CT.HPD)
}

# Get focal parameter names
temp_vars <- data.Vec %>% 
  select(-c(system_ID, mosquito_species, pathogen, Temperature, sample_num,
            etaL, muL, KL, V0, lf)) %>% 
  colnames() %>% c("muV")

# Initialize data frame
CT.HPD <- tibble(system_ID = c(), Temperature = c(), focal_var = c(),
                 sigmaH = c(), KH = c(), variable = c(),
                 HPD_width = c(), full_HPD_width = c(), rel_HPD_width = c())

# For each focal parameter, calculate relative HPD width
for (var_name in temp_vars) {
  print(paste0("Focal variable: ", var_name))
  # a) Set all but focal parameter to its posterior mean
  data.HPD.Vec <- data.Vec %>% 
    mutate(muV = 1/lf) %>% 
    select(system_ID, Temperature, sample_num, all_of(var_name), KL) %>% 
    full_join(mean.Vec %>% 
                mutate(muV = 1/lf) %>% 
                select(-all_of(var_name))%>% 
                distinct(),
              by = c("system_ID", "Temperature", "KL")) %>% 
    mutate(lf = 1/muV) %>%
    mutate(V0 = ifelse(sigmaV_f * deltaL < (1 / lf),
                       0,
                       KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL))))
  # Set up progress bar
  pb <- progress_bar$new(
    format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
    total = length(unique(data.CTHPD$KH)),
    width = 120)    
  for (index_KH in unique(data.CTHPD$KH)) {
    pb$tick()
    # b) Get posterior samples of R0 (as a function of temperature)
    temp_df <- expand_grid(dplyr::filter(data.CTHPD, KH == index_KH), 
                           data.HPD.Vec) %>%
      mutate(RV = ifelse(is.infinite(sigmaH),
                         sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                         sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
      mutate(bV = ifelse(is.infinite(sigmaH),
                         sigmaV, # Ross-Macdonald model
                         sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
      mutate(RH = ifelse(V0 == 0,
                         0,
                         bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
      # Basic reproduction number
      mutate(R0 = sqrt(RV*RH)) %>%
      # Filter to maximum value of R0
      group_by(system_ID, sample_num, sigmaH, KH) %>%
      dplyr::filter(R0 > 1) %>%
      # Get lowest temperature at which R0 exceeds one
      mutate(CTmin = min(Temperature)) %>%
      # Get highest temperature at which R0 exceeds one
      mutate(CTmax = max(Temperature)) %>%
      # Get width of critical thermal interval
      mutate(CTwidth = CTmax - CTmin) %>% 
      select(system_ID, sample_num, sigmaH, KH, CTmin, CTmax, CTwidth) %>% 
      distinct() %>% 
      replace_na(list(CTmin = Inf,
                      CTmax = -Inf,
                      CTwidth = 0)) %>% 
      pivot_longer(cols = c(CTmin, CTmax, CTwidth), names_to = "variable", values_to = "value") %>% 
      group_by(system_ID, sigmaH, KH, variable) %>% 
      summarise(
        HPD_low = hdi(value, credMass = 0.95)[1],
        HPD_high = hdi(value, credMass = 0.95)[2],
        HPD_width = max(eps, HPD_high-HPD_low),
        .groups = "keep"
      ) %>% 
      select(system_ID, sigmaH, KH, variable, HPD_width) %>% 
      right_join(full.CT.HPD %>% 
                   select(-c(HPD_low, HPD_high)) %>% 
                   rename(full_HPD_width = HPD_width) %>% 
                   filter(KH == index_KH),
                 by = c("system_ID", "sigmaH", "KH", "variable")) %>% 
      mutate(focal_var = var_name) %>% 
      group_by(system_ID, sigmaH, KH, focal_var) %>% 
      # d) Normalize HPD width by the full HPD width when all parameters are allowed to vary
      mutate(rel_HPD_width = ifelse(full_HPD_width < 1.1 * eps, 0, HPD_width / full_HPD_width))
    
    CT.HPD <- rbind(CT.HPD, temp_df)
  }
  pb$terminate()
}

# Save CT relative highest posterior density data
write_rds(CT.HPD, "results/CT_HPD_unc.rds", compress = 'gz')

# 8) Densities of Topt, CT min, max, and width ----------------------------

###* Density of Topt ----
# Set up host trait dataframe
data.Topt <- data.Host %>%
  dplyr::filter(sigmaH %in% c(100, Inf)) %>%
  dplyr::filter(KH %in% 10^seq(-2,4))

Topt.density.df <- expand_grid(data.Topt, data.Vec) %>% 
  data.table::data.table() %>%
  mutate(RV = ifelse(is.infinite(sigmaH),
                     sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                     sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
  mutate(bV = ifelse(is.infinite(sigmaH),
                     sigmaV, # Ross-Macdonald model
                     sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
  mutate(RH = ifelse(V0 == 0,
                     0,
                     bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(RV*RH)) %>%
  # filter to maximum value of R0
  dplyr::filter(R0>0) %>%
  select(system_ID, sample_num, Model, sigmaH, KH, Temperature, R0) %>% 
  group_by(system_ID, sample_num, sigmaH, KH) %>%
  dplyr::filter(R0 == max(R0)) %>%
  distinct() %>%
  # Get temperature at which R0 is maximized
  rename(Topt = Temperature) %>%
  dplyr::select(system_ID, sample_num, Model, sigmaH, KH, Topt)

# Add back in rows removed from filtering
Topt.density.df <- full_join(Topt.density.df,
                             expand_grid(data.Topt, data.Vec) %>%
                               select(c(Model, system_ID, sample_num, sigmaH,KH)) %>%
                               distinct())

###* Densities of CTmin,max,width ----

# Initialize data frame
CT.density.df <- init.df
# Set up host trait data frame
data.CT <- data.Host %>%
  dplyr::filter(sigmaH %in% c(100, Inf)) %>%
  dplyr::filter(KH %in% 10^seq(-2,4))

# Slice host trait data
cluster_size <- parallel::detectCores() - 1
sigmaH_slices <- slice_func(unique(data.CT$sigmaH), 2)
KH_slices <- slice_func(unique(data.CT$KH), cluster_size)

CT.density.df <- expand_grid(data.CT, data.Vec) %>% 
  data.table::data.table() %>%
  mutate(RV = ifelse(is.infinite(sigmaH),
                     sigmaV * betaH / (1 / (lf + eps)), # Ross-Macdonald
                     sigmaH * sigmaV * betaH * KH / ((1 / (lf + eps)) * (sigmaH * KH + sigmaV * V0)))) %>%
  mutate(bV = ifelse(is.infinite(sigmaH),
                     sigmaV, # Ross-Macdonald model
                     sigmaV * sigmaH * KH / (sigmaH * KH + sigmaV * V0 + eps))) %>%
  mutate(RH = ifelse(V0 == 0,
                     0,
                     bV * betaV * V0 * exp(-1 / (lf * etaV)) / (KH * (gammaH + muH) + eps))) %>%
  # Basic reproduction number
  mutate(R0 = sqrt(RV*RH)) %>%
  # filter to temperatures where R0 exceeds one
  group_by(system_ID, sample_num, sigmaH, KH) %>%
  dplyr::filter(R0 > 1) %>%
  # Get lowest temperature at which R0 exceeds one
  mutate(CTmin = min(Temperature)) %>%
  # Get highest temperature at which R0 exceeds one
  mutate(CTmax = max(Temperature)) %>%
  # Get width of critical thermal interval
  mutate(CTwidth = CTmax - CTmin) %>%
  dplyr::select(system_ID, sample_num, Model, sigmaH, KH, CTmin, CTmax, CTwidth) %>% 
  distinct()

# Combine CT and Topt data frames and add informative labels
all.density.df <- right_join(Topt.density.df, CT.density.df) %>% 
  select(-CTwidth) %>% 
  pivot_longer(cols = c(CTmin, Topt, CTmax), names_to = "variable", values_to = "value") %>% 
  mutate(system_label = case_when(
    system_ID == "Aedes aegypti / DENV" ~ "(a) *Ae. aegypti* and DENV",
    system_ID == "Aedes aegypti / ZIKV" ~ "(b) *Ae. aegypti* and ZIKV",
    system_ID == "Aedes albopictus / DENV" ~ "(c\u200b) *Ae. albopictus* and DENV",
    system_ID == "Anopheles gambiae / Plasmodium falciparum" ~ "(d) *An. gambiae* and *P. falciparum*",
    system_ID == "Culex quinquefasciatus / WNV" ~ "(e) *Cx. quinquefasciatus* and WNV"
  )) %>%
  mutate(variable = case_when(
    variable == "CTmin" ~ "Critical thermal minimum",
    variable == "Topt" ~ "Thermal optimum",
    variable == "CTmax" ~ "Critical thermal maximum"
  )) %>% 
  mutate(KH_factor = factor(KH, levels = c("0.01", "0.1", "1", "10", "100", "1000", "10000")))

all.density.df$variable <- factor(all.density.df$variable, 
                                  levels = c("Critical thermal minimum", "Thermal optimum","Critical thermal maximum"))
all.density.df$sigmaH <- factor(all.density.df$sigmaH, 
                                levels = c(100, Inf))

write_rds(all.density.df, "results/CTminmaxTopt_densities.rds")
