## Title: Calculate Topt and critical thermal extrema #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Transform trait values, which we have fit TPCs to from data, in
##          order to use in model analyses
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Calculate R0 TPCs
##           3) Calculate Topt
##           4) Calculate Topt alternate: restrict to R0 > 1
##           5) Calculate CTmin, CTmax, and CTwidth
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

eps <- .Machine$double.eps

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# Initialize a data frame to save our analyses
init.df <- tibble(system_ID = character(), Temperature = double(), 
                  Model = character(), sigmaH = double(), KH = double(), 
                  variable = character(), lowHCI = double(), highHCI = double(),
                  mean = double(), median = double())


# Load in vector trait data
data.Vec <- read_rds("results/VecTPC_vals.rds")
# # Thin out the vector trait data samples
# thin_size <- 2000
# samples <- unique(data.Vec$sample_num)
# num_samples <- length(samples)
# sample_inds <- sample(samples, min(num_samples, thin_size), replace = FALSE)
# 
# data.Vec <- data.Vec %>%
#   filter(sample_num %in% sample_inds)

# Load in host trait data
data.Host <- read_rds("results/Host_vals.rds")

# 1) Define accessory functions -------------------------------------------

# Function: Slice data to optimize usage of parallel processing and memory
slice_func<-function(x,n) {
  N<-length(x);
  lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

# Function: evaluate R0 (as a function of temperature), its mean, median, and 
#           highest prob. intervals across samples
R0_TPC_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>%
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
    mutate(R0 = sqrt(RV*RH))%>%
    dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, R0) %>%
    # Normalize R0 across temperature
    group_by(system_ID, Model, sigmaH, KH) %>%
    ungroup() %>%
    pivot_longer(cols = c(R0), names_to = "variable", values_to = "value") %>%
    group_by(system_ID, Temperature, Model, sigmaH, KH, variable) %>%
    # partition(cluster) %>%
    summarise(
      lowHCI = quantile(value, 0.055),
      highHCI = quantile(value, 0.945),
      mean = mean(value),
      median = median(value)
    ) #%>%
  # collect()
}

# Function: evaluate Topt and its mean, median, and highest prob. intervals across
#           samples
Topt_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>% 
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
    group_by(system_ID, sigmaH, KH, sample_num) %>%
    filter(R0 == max(R0)) %>%
    # Get temperature at which R0 is maximized
    rename(Topt = Temperature) %>%
    dplyr::select(system_ID, Model, sigmaH, KH, Topt, sample_num) %>%
    ungroup() %>%
    group_by(system_ID, Model, sigmaH, KH) %>%
    summarise(
      lowHCI = quantile(Topt, 0.055),
      highHCI = quantile(Topt, 0.945),
      mean = mean(Topt),
      median = median(Topt),
      .groups = "keep"
    )
}

# Function: evaluate CTmin, CTmax, CTwidth and their means, medians, and highest
#           prob. intervals across samples
CT_heat_func <- function(in_df, system_name) {
  out_df <- in_df %>%
    expand_grid(filter(data.Vec, system_ID == system_name)) %>%
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
    filter(R0 > 1) 
  
  # If R0 < 1 across all temperatures, report the following:
  #  CTmin = Inf, CTmax = -Inf, CTwidth = 0
  if (dim(out_df)[1] == 0) {
    out_df <- expand_grid(select(in_df, sigmaH, KH, Model),
                          system_ID = system_name,
                          variable = c("CTmax", "CTmin", "CTwidth")) %>%
      distinct() %>%
      mutate(mean = case_when(
        variable == "CTmax" ~ -Inf,
        variable == "CTmin" ~ Inf,
        variable == "CTwidth" ~ 0
      )) %>%
      mutate(median = mean) %>%
      mutate(highHCI = mean) %>%
      mutate(lowHCI = mean)
  } else {
    out_df <- out_df %>%
      # Get lowest temperature at which R0 exceeds one
      mutate(CTmin = min(Temperature)) %>%
      # Get highest temperature at which R0 exceeds one
      mutate(CTmax = max(Temperature)) %>%
      # Get width of critical thermal interval
      mutate(CTwidth = ifelse(is.finite(CTmax), 
                              CTmax - CTmin,
                              0))%>%
      ungroup() %>%
      dplyr::select(system_ID, Model, sigmaH, KH, CTmin, CTmax, CTwidth) %>%
      pivot_longer(cols = c(CTwidth, CTmin, CTmax), names_to = "variable", values_to = "value") %>% 
      group_by(system_ID, Model, sigmaH, KH, variable) %>%
      summarise(
        lowHCI = quantile(value, 0.055),
        highHCI = quantile(value, 0.945),
        mean = mean(value),
        median = median(value),
        .groups = "keep"
      )
  }
  return(out_df)
}

# Function: Evaluate local derivatives of R0 with respect to mosquito traits
R0_deriv_func <- function(in_df) {
  
  sigmaH_val <- unique(in_df$sigmaH)
  
  in_df <- 
    
    # Infinite sigmaH case
    if (is.infinite(sigmaH_val)) {
      out_df <- in_df %>% 
        mutate(lf = ifelse(is.infinite(muV), 0, 1/muV)) %>% 
        # Derivative of V0 wrt focal variable
        mutate(V0_deriv = case_when(
          focal_var == "sigmaV" ~ rhoL * KL / (sigmaV_f * deltaL * (sigmaV)^2 + eps),
          focal_var == "sigmaV_f" ~ (rhoL * KL)/ (deltaL * (sigmaV_f^2) + eps),
          focal_var == "deltaL" ~ (rhoL * KL)/ (sigmaV_f * (deltaL^2) + eps),
          focal_var == "rhoL" ~ KL * (lf - (1/(sigmaV_f + eps))),
          focal_var == "muV" ~ -1 * rhoL * KL * lf^2,
          TRUE ~ 0
        )) %>% 
        # Derivative of R0 wrt focal variable
        mutate(R0_deriv = case_when(
          focal_var == "sigmaV" ~ 0.5 * (sigmaV * V0_deriv + 2 * V0) / (sigmaV * V0 + eps),
          focal_var == "betaV" ~ 0.5 / (betaV + eps),
          focal_var == "sigmaV_f" ~ 0.5 * V0_deriv / (V0 + eps),
          focal_var == "deltaL" ~ 0.5 * V0_deriv / (V0 + eps),
          focal_var == "rhoL" ~ 0.5 * V0_deriv / (V0 + eps),
          focal_var == "etaV" ~ 0.5 * muV / (etaV^2 + eps),
          focal_var == "muV" ~ ifelse(etaV > 0 & lf > 0,
                                      0.5 * ((V0_deriv / (V0 + eps)) - (1/(etaV + eps)) - lf),
                                      0),
          TRUE ~ 0
        )) %>% 
        mutate(Temp_deriv = R0_deriv * dpardT)
      
      # Finite sigmaH case
    } else {
      out_df <- in_df %>% 
        mutate(lf = ifelse(is.infinite(muV), 0, 1/muV)) %>% 
        # Derivative of V0 wrt focal variable
        mutate(V0_deriv = case_when(
          focal_var == "sigmaV" ~ rhoL * KL / (sigmaV_f * deltaL * (sigmaV)^2 + eps),
          focal_var == "sigmaV_f" ~ (rhoL * KL)/ (deltaL * (sigmaV_f^2) + eps),
          focal_var == "deltaL" ~ (rhoL * KL)/ (sigmaV_f * (deltaL^2) + eps),
          focal_var == "rhoL" ~ KL * (lf - (1/(sigmaV_f + eps))),
          focal_var == "muV" ~ -1 * rhoL * KL * lf^2,
          TRUE ~ 0
        ))  %>% 
        # Intermediate quantity
        mutate(K = V0 * (sigmaH * KH + sigmaV * V0 + eps)^(-2)) %>% 
        # Intermediate derivative
        mutate(K_deriv = (sigmaH * KH + sigmaV * V0 + eps)^(-2) * case_when(
          focal_var == "sigmaV" ~ V0_deriv - 2 * V0 * (sigmaV * V0_deriv + V0) / (sigmaH * KH + sigmaV * V0),
          focal_var %in% c("betaV", "etaV") ~ 0,
          TRUE ~ V0_deriv * (1 - 2 * sigmaV * V0 / (sigmaH * KH + sigmaV * V0))
        )) %>% 
        # Derivative of R0 wrt focal variable
        mutate(R0_deriv = case_when(
          focal_var == "sigmaV" ~ (1 / (sigmaV + eps) + 0.5 * K_deriv / (K+eps)),
          focal_var == "betaV" ~ 0.5 / (betaV + eps),
          focal_var == "sigmaV_f" ~ 0.5 * K_deriv / (K + eps),
          focal_var == "deltaL" ~ 0.5 * K_deriv / (K + eps),
          focal_var == "rhoL" ~ 0.5 * K_deriv / (K + eps),
          focal_var == "etaV" ~ 0.5 * muV / (etaV^2 + eps),
          focal_var == "muV" ~ ifelse(etaV > 0 & lf > 0,
                                      0.5 * ((K_deriv / (K + eps)) - (1/(etaV + eps)) - lf),
                                      0),
          TRUE ~ 0
        )) %>% 
        mutate(Temp_deriv = R0_deriv * dpardT)
    }
}

# 2) Calculate R0 TPCs -----------------------------------------

# # Set up parallel
# if (!exists("cluster")) {
#   cluster_size <- min(21, parallel::detectCores() - 2)
#   cluster <- new_cluster(cluster_size)
#   cluster_library(cluster, c("dplyr", "tidyr"))
# }

# Start new cluster for doParallel
cluster_size <- parallel::detectCores()-1
my.cluster <- parallel::makeCluster(cluster_size, type = "PSOCK"
)
# Register cluster for doParallel
doSNOW::registerDoSNOW(cl = my.cluster)

# Set up host trait data frame
data.R0 <- as.data.frame(data.Host) %>%
  filter(sigmaH %in% c(10^seq(-1,2), Inf,
                       unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 51)])) %>%
  filter(KH %in% c(10^seq(-2,5) ,
                   unique(KH)[seq(1, length(unique(KH)), length.out = 51)]))

# Initialize R0 data frame
R0.df <- init.df

gc()

# # Collect R0 TPC data across systems and host trait values
# for (i in 1:length(unique(data.Vec$system_ID))) {
#   R0.df <- init.df
#   system_name = unique(data.Vec$system_ID)[i]
#   print(paste0("R0 vals: ", system_name))
#   KHslice_num <- 1
#   # if (exists("pb")) {pb$terminate();rm(pb)}
#   pb <- progress_bar$new(
#     format = ":spin :system progress = :percent [:bar] :elapsedfull | eta: :eta",
#     total = length(KH_slices) * length(sigmaH_slices),
#     width = 120)  
#   for(index_KH in KH_slices) {
#     sigmaHslice_num <- 1
#     for (index_sigmaH in sigmaH_slices) {
#       temp_df <- data.R0 %>%
#         filter(sigmaH %in% index_sigmaH,
#                KH %in% index_KH) %>%
#         R0_TPC_func(., system_name)
#       
#       R0.df <- rbind(temp_df, R0.df)
#       rm(temp_df)
#       gc()
#       
#       pb$tick()
#       sigmaHslice_num <- sigmaHslice_num  + 1
#     }
#     KHslice_num <- KHslice_num +1
#   }
#   
#   # Save data separately to reduce RAM usage
#   R0_TPCs.df <- R0.df %>% 
#     ungroup() %>% 
#     filter(sigmaH %in% c(100, Inf))
#   save_name <- paste0("results/R0_TPC_vals_", i,".rds")
#   write_rds(R0_TPCs.df, save_name, compress = "gz")
#   
#   R0_heat.df <- R0.df %>% 
#     ungroup() %>% 
#     filter(KH %in% c(0.1, 1, 10, 100))
#   save_name <- paste0("results/R0_vals_", i,".rds")
#   write_rds(R0_heat.df, save_name, compress = "gz")
#   
#   rm(R0.df)
#   gc()
# }

# Collect R0 TPC data across systems and host trait values
for (i in 1:length(unique(data.Vec$system_ID))) {
  
  R0.df <- init.df
  system_name = unique(data.Vec$system_ID)[i]
  print(paste0("R0 vals: ", system_name))
  # if (exists("pb")) {pb$terminate();rm(pb)}
  
  # Slice host trait data
  sigmaH_slices <- slice_func(unique(data.R0$sigmaH), 1)
  KH_slices <- slice_func(unique(data.R0$KH), floor(cluster_size/3))
  
  # Set up iteration grid
  iter_grid <- expand_grid(system_ID = system_name,
                           tibble(KH = KH_slices)) %>% 
    expand_grid(sigmaH = sigmaH_slices)
  
  # Set up progress bar
  iterations <- dim(iter_grid)[1]
  pb <- progress_bar$new(
    format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
    total = iterations,
    width = 120)                                                                                                        
  progress <- function(n){
    pb$tick()
  }
  opts <- list(progress = progress)
  
  # Calculate Topt values
  R0.df <- foreach(
    # system_name = iter_grid$system_ID,
    index_KH = iter_grid$KH, #seq(1,length(KH_slices)),
    index_sigmaH = iter_grid$sigmaH, #seq(1,length(sigmaH_slices)),
    .packages = "tidyverse",
    .combine = rbind,
    .options.snow = opts) %dopar% {
      # sigmaH_in = sigmaH_slices[[index_sigmaH]]
      # KH_in = KH_slices[[index_KH]]
      sigmaH_in = index_sigmaH
      KH_in = index_KH
      # tibble(sigmaH = sigmaH_in, KH = KH_in)
      
      data.R0 %>%
        filter(sigmaH %in% sigmaH_in,
               KH %in% KH_in) %>%
        R0_TPC_func(., system_name)
    }
  
  pb$terminate()
  # Save data separately to reduce RAM usage
  R0_TPCs.df <- R0.df %>% 
    ungroup() %>% 
    filter(sigmaH %in% c(100, Inf))
  save_name <- paste0("results/R0_TPC_vals_", i,".rds")
  write_rds(R0_TPCs.df, save_name, compress = "gz")
  
  R0_heat.df <- R0.df %>% 
    ungroup() %>% 
    filter(KH %in% c(0.1, 1, 10, 100))
  save_name <- paste0("results/R0_vals_", i,".rds")
  write_rds(R0_heat.df, save_name, compress = "gz")
  
  rm(R0.df)
  gc()
}

# Save data
# proper_dim <- dim(data.R0)[1] * length(unique(data.Vec$system_ID)) * length(unique(data.Vec$Temperature))
# 
# (dim(R0.df)[1] == proper_dim)

R0_TPCs.df <- init.df
R0_heat.df <- init.df

# Combine data into a single file
for (i in 1:length(unique(data.Vec$system_ID))) {
  save_name <- paste0("results/R0_TPC_vals_", i,".rds")
  temp_df <- read_rds(save_name)
  R0_TPCs.df <- rbind(R0_TPCs.df, temp_df)
  
  save_name <- paste0("results/R0_vals_", i,".rds")
  temp_df <- read_rds(save_name)
  R0_heat.df <- rbind(R0_heat.df, temp_df)
}

write_rds(R0_TPCs.df, "results/R0_TPC_vals.rds", compress = "gz")
write_rds(R0_heat.df, "results/R0_vals.rds", compress = "gz")
rm(R0_TPCs.df,R0_heat.df)
# Reset cluster
rm(cluster)
gc()
# 3) Calculate R0 derivatives ---------------------------------------------

# # Reset cluster
# rm(cluster)
# # Set up parallel
# cluster_size <- parallel::detectCores() - 2
# cluster <- new_cluster(cluster_size)
# cluster_library(cluster, c("dplyr", "tidyr"))

# Set up host trait data frame (for future visualization)
data.dR0dpar <- data.Host %>%
  filter(sigmaH %in% c(10, Inf)) %>%
  filter(KH %in% 10^seq(-2,3))#-2,5))

# Slice host trait data
sigmaH_slices <- slice_func(unique(data.dR0dpar$sigmaH), 1)
KH_slices <- slice_func(unique(data.dR0dpar$KH), 1)

# Initialize R0 data frame
dR0dvar.df <- init.df %>% rename(focal_var = variable)

gc()

# Parameters of interest
temp_vars <- c("sigmaV", "betaV", "sigmaV_f", "deltaL", "rhoL", "etaV", "muV")

# Collect R0 TPC data across systems and host trait values
print_index = 1
for (system_name in unique(data.Vec$system_ID)) {
  print(paste0("(",print_index, "/", length(unique(data.Vec$system_ID)),") dR0/dvar: ", system_name))
  
  pb <- progress_bar$new(
    format = ":spin :system progress = :percent [:bar] :elapsedfull | eta: :eta",
    total = length(unique(data.dR0dpar$KH)) * length(sigmaH_slices),
    width = 120)  
  for(index_KH in KH_slices) {
    for (index_sigmaH in sigmaH_slices) {
      pb$tick()
      
      dpardT_df <- data.Vec %>%  
        filter(system_ID == system_name) %>%
        ungroup() %>%
        select(-c(KL, V0, etaL, muL)) %>% 
        mutate(muV = ifelse(lf != 0, 1/ (lf + eps), Inf), .keep = "unused") %>% 
        # mutate(muV = ifelse(muV > 1/eps, Inf, muV)) %>% 
        pivot_longer(cols = all_of(temp_vars), 
                     names_to = "focal_var") %>% 
        group_by(system_ID, mosquito_species, pathogen, sample_num, focal_var) %>% 
        arrange(system_ID, sample_num, focal_var, Temperature) %>% 
        mutate(lag_val = lag(value)) %>% 
        # If trait is "infinite" (say because lifespan is 0), then there is no
        # change, so derivative is zero
        mutate(numerator = ifelse(is.infinite(value) & is.infinite(lag(value)),
                                  0,
                                  value - lag(value))
        ) %>% 
        mutate(dpardT = (numerator / (Temperature - lag(Temperature)))) %>% 
        mutate(dpardT = ifelse(dpardT > 1/(10 * eps), 0, dpardT))
      
      
      temp_df <- data.dR0dpar %>%
        filter(sigmaH %in% index_sigmaH,
               KH %in% index_KH) %>%
        expand_grid(filter(data.Vec, system_ID == system_name)) %>%
        cross_join(data.frame(focal_var = temp_vars)) %>% 
        left_join(dpardT_df, by = join_by(system_ID, mosquito_species, pathogen, Temperature, sample_num, focal_var)) %>% 
        filter(Temperature > 10) %>% 
        mutate(muV = ifelse(lf != 0, 1/ (lf + eps), Inf), .keep = "unused") %>% 
        # mutate(muV = ifelse(muV > 1/eps, NA, muV)) %>% 
        # Calculate derivatives
        R0_deriv_func(.) %>% 
        filter(value != 0)
      
      temp_df2 <- temp_df %>%
        dplyr::select(system_ID, sample_num, Temperature, Model, sigmaH, KH, focal_var, value,  Temp_deriv) %>%
        group_by(system_ID, Temperature, Model, sigmaH, KH, focal_var) %>%
        # Filter out na values: these correspond to zero lifespan cases
        filter(is.finite(value)) %>%
        # partition(cluster) %>%
        summarise(
          lowHCI = quantile(Temp_deriv, 0.055, na.rm = TRUE),
          highHCI = quantile(Temp_deriv, 0.945, na.rm = TRUE),
          mean = mean(Temp_deriv),
          median = median(Temp_deriv),
          .groups = "drop"
        ) #%>%
      # collect()
      
      # Get total derivative
      totalderiv_df <- temp_df2 %>% 
        # select(system_ID, Model, Temperature, sigmaH, KH, sample_num, focal_var, Temp_deriv) %>% 
        group_by(system_ID, Model, Temperature, sigmaH, KH) %>% 
        summarise(median = sum(median),
                  .groups = "drop") %>% 
        mutate(focal_var = "total",
               lowHCI = NA, highHCI = NA, mean = NA)
      
      dR0dvar.df <- rbind(temp_df2, totalderiv_df) %>% 
        rbind(dR0dvar.df)
    }
  }
  gc()
  print_index <- print_index + 1
}
pb$terminate()

# Save data
proper_dim <- dim(data.dR0dpar)[1] * length(unique(data.Vec$system_ID)) * length(unique(data.Vec$Temperature)) * (length(temp_vars) + 1)

(dim(dR0dvar.df)[1] == proper_dim)

write_rds(dR0dvar.df, "results/dR0dk_vals.rds", compress = "gz")

rm(dR0dvar.df)

# 4) Calculate Topt -------------------------------------------------------

# Set up new cluster
gc()
# Start new cluster for doParallel
cluster_size <- parallel::detectCores()-1
my.cluster <- parallel::makeCluster(cluster_size, type = "PSOCK"
)
# Register cluster for doParallel
doSNOW::registerDoSNOW(cl = my.cluster)

# Set up host trait data frame
data.Topt <- as.data.frame(data.Host) %>%
  filter(sigmaH %in% c(10^seq(-1,2), Inf,
                       unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 51)])) %>%
  filter(KH %in% c(10^seq(-2,5) ,
                   unique(KH)[seq(1, length(unique(KH)), length.out = 51)]))

# Slice host trait data
sigmaH_slices <- slice_func(unique(data.Topt$sigmaH), 1)
KH_slices <- slice_func(unique(data.Topt$KH), floor(cluster_size/2))

# Set up iteration grid
iter_grid <- expand_grid(system_ID = unique(data.Vec$system_ID),
                         tibble(KH = KH_slices)) %>% 
  expand_grid(sigmaH = sigmaH_slices)

# Set up progress bar
iterations <- dim(iter_grid)[1]
pb <- progress_bar$new(
  format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
  total = iterations,
  width = 120)                                                                                                        
progress <- function(n){
  pb$tick(tokens = list(system = iter_grid$system_ID[n]))
}
opts <- list(progress = progress)

# Calculate Topt values
Topt.df <- foreach(
  system_name = iter_grid$system_ID,
  index_KH = iter_grid$KH,
  index_sigmaH = iter_grid$sigmaH,
  .packages = "tidyverse",
  .combine = rbind,
  .options.snow = opts) %dopar% {
    data.Topt %>%
      filter(sigmaH %in% index_sigmaH,
             KH %in% index_KH) %>%
      Topt_heat_func(., system_name)
  }

pb$terminate()


# Save Topt data
proper_dim <- (dim(iter_grid)[1])
dim(Topt.df)[1] == proper_dim

# if (exists("Topt.df") & dim(Topt.df)[1] == proper_dim) {
  write_rds(Topt.df, "results/Topt_vals.rds", compress = "gz")
# } else {
  # warning("No file written. Topt.df either empty or not complete.")
# }

# End cluster
stopCluster(my.cluster) 

rm(Topt.df)

# 5) Calculate CTmin, CTmax, and CTwidth ----------------------------------
# Set up new cluster
gc()
# Start new cluster for doParallel
cluster_size <- parallel::detectCores()-1

my.cluster <- parallel::makeCluster(
  cluster_size, 
  type = "PSOCK"
)
# Register cluster for doParallel
doSNOW::registerDoSNOW(cl = my.cluster)

# Set up host trait data frame
data.CT <- as.data.frame(data.Host) %>%
  filter(sigmaH %in% c(10^seq(-1,2), Inf,
                       unique(sigmaH)[seq(1, length(unique(sigmaH)), length.out = 51)])) %>%
  filter(KH %in% c(10^seq(-2,5) ,
                   unique(KH)[seq(1, length(unique(KH)), length.out = 51)]))

# Slice host trait data
sigmaH_slices <- slice_func(unique(data.CT$sigmaH), 1)
KH_slices <- slice_func(unique(data.CT$KH), floor(cluster_size/2))

# Set up iteration grid
iter_grid <- expand_grid(system_ID = unique(data.Vec$system_ID),
                         tibble(KH = KH_slices)
                         ) %>% 
  expand_grid(sigmaH = sigmaH_slices)

# Set up progress bar
iterations <- dim(iter_grid)[1]
pb <- progress_bar$new(
  format = ":spin :system progress = :percent [:bar] :elapsed | eta: :eta",
  total = iterations,
  width = 120)
progress <- function(n){
  pb$tick(tokens = list(system = iter_grid$system_ID[n]))
}
opts <- list(progress = progress)

# Calculate critical thermal minimum, maximum, and width
CT.df <- foreach(
  system_name = iter_grid$system_ID,
  index_KH = iter_grid$KH,
  index_sigmaH = iter_grid$sigmaH,
  .packages = "tidyverse",
  .combine = rbind,
  .options.snow = opts) %dopar% {
    data.CT %>%
      filter(sigmaH %in% index_sigmaH,
             KH %in% index_KH) %>%
      CT_heat_func(., system_name)
  }

# CT.df <- tibble(sigmaH = c(), KH = c(), Model = c(), system_ID = c(),
#                 variable = c(), mean = c(), median = c(), highHCI = c(), lowHCI = c())
# 
# progress_index = 1
# 
# for (system_name in iter_grid$system_ID) {
#   for (index_KH in iter_grid$KH) {
#     for (index_sigmaH in iter_grid$sigmaH) {
#       print(100 * progress_index/(dim(iter_grid)[1]))
#       temp_df <- data.Host %>%
#         filter(sigmaH %in% index_sigmaH,
#                KH %in% index_KH) %>%
#         CT_heat_func(., system_name)
#       
#       rbind(CT.df, temp_df)
#       progress_index = progress_index + 1
#     }
#   }
# }

pb$terminate()

# Save CT data
(proper_dim <- 3 * dim(iter_grid)[1])

(proper_dim == dim(CT.df)[1])

# if (exists("CT.df") & dim(CT.df)[1] == proper_dim) {
  write_rds(CT.df, "results/CT_vals.rds", compress = "gz")
# } else {
  # warning("No file written. CT.df either empty or not complete.")
# }

# End cluster
stopCluster(my.cluster) 

rm(CT.df)
gc()