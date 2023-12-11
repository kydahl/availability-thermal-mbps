## Title: Prior distributions of mosquito thermal traits #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Compute prior distributions for thermal trait parameters and provide
##          functions for sampling from these distributions
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Calculate prior distributions of thermal trait parameters from data
##           3) Save thermal trait parameter distributions
##           4) Functions to sample from thermal trait parameter distributions
##
##
## Inputs:  data - data/clean/data_for_TPC_fitting.csv
##                 data/clean/trait_TPC_forms.csv
##                 data/clean/gamma_fits.csv
##
##          code - code/Mordecai_2017/mcmc_utils_all.R
##                 code/Mordecai_2017/temp_functions_all.R
##
## Outputs: data - data/clean/TPC_param_samples.csv
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized February 2023
## Modified from original code provided in the following articles:
# * Mordecai, E. A., J. M. Caldwell, M. K. Grossman, C. A. Lippi, L. R. Johnson, M. Neira, J. R. Rohr, S. J. Ryan, V. Savage, M. S. Shocket, R. Sippy, A. M. Stewart Ibarra, M. B. Thomas, and O. Villena. 2019. Thermal biology of mosquito-borne disease. Ecology letters 22:1690–1708.
# * Mordecai, E. A., J. M. Cohen, M. V. Evans, P. Gudapati, L. R. Johnson, C. A. Lippi, K. Miazgowicz, C. C. Murdock, J. R. Rohr, S. J. Ryan, V. Savage, M. S. Shocket, A. Stewart Ibarra, M. B. Thomas, and D. P. Weikel. 2017. Detecting the impact of temperature on transmission of Zika, dengue, and chikungunya using mechanistic models. PLoS neglected tropical diseases 11:e0005568.
# * Mordecai, E. A., K. P. Paaijmans, L. R. Johnson, C. Balzer, T. Ben-Horin, E. de Moor, A. McNally, S. Pawar, S. J. Ryan, T. C. Smith, and K. D. Lafferty. 2013. Optimal temperature for malaria transmission is dramatically lower than previously predicted. Ecology letters 16:22–30.
# * Shocket, M. S., S. J. Ryan, and E. A. Mordecai. 2018. Temperature explains broad patterns of Ross River virus transmission. eLife 7.
# * Shocket, M. S., A. B. Verwillow, M. G. Numazu, H. Slamani, J. M. Cohen, F. El Moustaid, J. Rohr, L. R. Johnson, and E. A. Mordecai. 2020. Transmission of West Nile and five other temperate mosquito-borne viruses peaks at temperatures between 23°C and 26°C. eLife 9.
# * Tesla, B., L. R. Demakovsky, E. A. Mordecai, S. J. Ryan, M. H. Bonds, C. N. Ngonghala, M. A. Brindley, and C. C. Murdock. 2018. Temperature drives Zika virus transmission: evidence from empirical and mathematical models. Proceedings. Biological sciences / The Royal Society 285:20180795.
# _______________________________________________________________________________

# 0) Set-up, load in necessary packages and data-sets ---------------------

# Load Libraries
library(tidyverse)
library(IDPmisc)
library(rjags) # Make sure you have installed JAGS-4.x.y.exe (for any x >=0, y>=0) from http://www.sourceforge.net/projects/mcmc-jags/files
library(MASS)
library(here)
library(retry)


# Load utility functions (from Mordecai et al., 2017)
# This file contains tools for analysis and visualization.
source("code/Mordecai_2017/mcmc_utils_all.R")

# This file contains the thermal response functions and their derivatives.
source("code/Mordecai_2017/temp_functions_all.R")

### * Load in data ----
# Trait data
data.in <- data.in.TPC

# 1) Define accessory functions -------------------------------------------

# Function: get samples of thermal trait parameters
thermtrait.prior.sample <- function(data_in, trait_in, mosquito_in, pathogen_in,
                                    n.chains = 5, n.adapt = 1000, n.samps = 1000,
                                    old_informative = FALSE) {
  # restrict the dataset to the trait we care about
  working_data <- dplyr::filter(data_in, trait.name == trait_in)
  
  # check if the trait is a probability
  prob_bool <- trait_in %in% c("b", "c", "bc", "e2a", "pLA", "pRH", "pO", "EV")
  
  # Get the proper TPC function
  TPC_function <- read_csv("data/clean/trait_TPC_forms.csv", show_col_types = FALSE) %>%
    # when the TPC function can't be found from literature (this occurs for 10 systems)
    # use Briere (less restrictive than quadratic)
    mutate(TPC.function = ifelse(is.na(TPC.function), "Briere", TPC.function)) %>%
    # For cases with multiple possible TPC functions, choose the one used in later studies
    mutate(TPC.function = str_extract(TPC.function, "[^/]+$")) %>%
    dplyr::filter(
      trait.name == trait_in,
      mosquito_species == mosquito_in,
      pathogen == pathogen_in
    ) %>%
    dplyr::select(TPC.function) %>%
    as.character()
  
  # set initial values for jags model
  inits_list <- if (TPC_function == "Briere") {
    if (trait_in == "EFD") { # use revised initial values from Mordecai 2017
      list(T0 = 15, Tm = 34, c = 0.01)
    }  else  if (trait_in == "bc") { # use revised initial values from Mordecai 2017
      list(T0 = 20, Tm = 35, c = 1e-4)
    } else if (trait_in == "EV") {
      list(T0 = 5, Tm = 31, c = 0.00007)
    } else {
      list(T0 = 5, Tm = 38, c = 0.00007)
    }
  } else if (TPC_function == "Quadratic") {
    list(c = 0.005, Tm = 33, T0 = 5)
  } else if (TPC_function == "AltQuadratic") {
    list(c = 0.0001, Tm = 0.5, T0 = 0.04)
  } else if (TPC_function == "Linear") {
    list(z = 0.2, c = 0.005) # corresponds to Tm = 40
  }
  
  # designate the TPC hyperparameters to be fit
  variable_names <- if (TPC_function %in% c("Briere", "Quadratic", "AltQuadratic")) {
    c("c", "Tm", "T0", "sigma")
  } else if (TPC_function == "Linear") {
    c("c", "z", "sigma")
  }
  
  # If you choose not to use previously published TPC hyperparameters from Mordecai et al., 2017:
  # 1) use default priors
  # 2) update these using any related species
  # 3) sequentially update with more recent datasets
  if (old_informative == FALSE) {
    # Gather data from other species of the same genus
    other_species <- working_data %>%
      dplyr::filter(stringr::word(mosquito_species, 1, 1) == stringr::word(mosquito_in, 1, 1)) %>%
      dplyr::filter(mosquito_species != mosquito_in) %>%
      distinct(T, trait, .keep_all = TRUE)
    
    # If data was recorded for infections, use infections of the same pathogen
    # in other species to inform priors
    if (pathogen_in != "none") {
      other_species <- working_data %>%
        dplyr::filter(mosquito_species != mosquito_in) %>%
        dplyr::filter(stringr::word(pathogen, 1, 1) == stringr::word(pathogen_in, 1, 1)) %>%
        rbind(other_species)
    }
    
    # if data from congeners is unavailable, use data from species outside of the focal genera (that is, not Aedes, Culex, or Anopheles)
    if (dim(other_species)[1] == 0) {
      other_species <- working_data %>%
        dplyr::filter(mosquito_species == "Other spp.")
    }
    
    # get hyperparameters of TPC parameter priors using this data, if available
    if (dim(other_species)[1] > 0) {
      other_data <- list(
        "Y" = other_species$trait, "T" = other_species$T,
        "N" = length(other_species$T)
      )
      
      jags_other <- if (prob_bool) {
        # for probabilities (truncates at 1)
        case_when(
          TPC_function == "Briere" ~ "code/jags-models/jags-briere-trunc.bug",
          TPC_function == "Quadratic" ~ "code/jags-models/jags-quad-neg-trunc.bug",
          TPC_function == "Linear" ~ "code/jags-models/jags-linear-trunc.bug"
        )
      } else {
        case_when(
          TPC_function == "Briere" & trait_in != "EFD" ~ "code/jags-models/jags-briere.bug",
          TPC_function == "Briere" & trait_in == "EFD" ~ "code/jags-models/jags-briere-EFD.bug",
          TPC_function == "Quadratic" ~ "code/jags-models/jags-quad-neg.bug",
          TPC_function == "AltQuadratic" ~ "code/jags-models/jags-quad-alt.bug",
          TPC_function == "Linear" ~ "code/jags-models/jags-linear.bug"
        )
      }
      
      print(paste0("Getting hyperparameters for informed prior distribution from ", dim(distinct(other_species, mosquito_species, pathogen))[1], " other systems..."))
      # get hyperparameters of informed TPC parameter distributions
      # by resampling, this function repeats the fitting process until it converges
      prev_hypers <- retry(
        get.prior_hyperparams(
          other_data, TPC_function, variable_names,
          jags_other, inits_list,
          n.chains, n.adapt, 500, prob_bool
        ),
        until = function(val, cnd) {
          !is.null(val)
        }
      )
    } else { # if no trait data is available for any other species, start with uninformed priors
      prev_hypers <- c()
    }
    
    # Use previously published TPC hyperparameters from Mordecai et al., 2017
  } else if (old_informative == TRUE) {
    # get previous hyperparameters
    prev_hypers <- read_csv("data/clean/gamma_fits.csv") %>%
      dplyr::filter(trait %in% trait_in) %>%
      unique() %>%
      # adjust by the appropriate multiplier
      mutate(value = multiplier * value) %>%
      dplyr::select(-multiplier) %>%
      # all the rest of this is to get the data back in the form that "jags" wants
      dplyr::select(-trait) %>%
      pivot_wider(names_from = Var2, values_from = value) %>%
      as.data.frame() %>%
      `rownames<-`(.[, 1]) %>%
      dplyr::select(-Var1)
    
    if (dim(prev_hypers)[1] == 0) {
      stop("No saved TPC hyperparameter data. Switch old_informative to false")
    }
    
    jags_choice <- if (prob_bool) {
      # probabilities (truncates at 1)
      case_when(
        TPC_function == "Briere" ~ "code/jags-models/jags-briere-trunc-informative.bug",
        TPC_function == "Quadratic" ~ "code/jags-models/jags-quad-neg-trunc-informative.bug",
        TPC_function == "Linear" ~ "code/jags-models/jags-linear-trunc-informative.bug"
      )
    } else {
      case_when(
        TPC_function == "Briere" & trait_in != "EFD" ~ "code/jags-models/jags-briere-informative.bug",
        TPC_function == "Briere" & trait_in == "EFD" ~ "code/jags-models/jags-briere-EFD-informative.bug",
        TPC_function == "Quadratic" ~ "code/jags-models/jags-quad-neg-informative.bug",
        TPC_function == "AltQuadratic" ~ "code/jags-models/jags-quad-alt-informative.bug",
        TPC_function == "Linear" ~ "code/jags-models/jags-linear-informative.bug"
      )
    }
    
    jags_data <- list(
      "Y" = working_data$trait, "T" = working_data$T,
      "N" = length(working_data$T), "hypers" = prev_hypers
    )
    
    samps <- run.jags(
      jags_data, TPC_function, variable_names,
      jags_choice, inits_list,
      n.chains, n.adapt, n.samps, prob_bool
    )
  }
  # Begin using data on the focal species to create and generate samples from
  # the TPC parameter posterior distributions
  
  # get all of the studies with data reported for the particular trait and mosquito species
  proper_studies <- working_data %>%
    dplyr::filter(mosquito_species == mosquito_in, pathogen %in% c(pathogen_in, "none")) %>%
    distinct(trait,T)
  
  
  # Select the appropriate bugs model
  jags_choice <- if (prob_bool) {
    # probabilities (truncates at 0 and 1)
    case_when(
      TPC_function == "Briere" & is.null(prev_hypers) ~ "code/jags-models/jags-briere-trunc.bug",
      TPC_function == "Briere" & !is.null(prev_hypers) ~ "code/jags-models/jags-briere-trunc-informative.bug",
      TPC_function == "Quadratic" & is.null(prev_hypers) ~ "code/jags-models/jags-quad-neg-trunc.bug",
      TPC_function == "Quadratic" & !is.null(prev_hypers) ~ "code/jags-models/jags-quad-neg-trunc-informative.bug",
      TPC_function == "Linear" & is.null(prev_hypers) ~ "code/jags-models/jags-linear-trunc.bug",
      TPC_function == "Linear" & !is.null(prev_hypers) ~ "code/jags-models/jags-linear-trunc-informative.bug"
    )
  } else {
    case_when(
      TPC_function == "Briere" & is.null(prev_hypers) & trait_in != "EFD" ~ "code/jags-models/jags-briere.bug",
      TPC_function == "Briere" & is.null(prev_hypers) & trait_in == "EFD" ~ "code/jags-models/jags-briere-EFD.bug",
      TPC_function == "Briere" & !is.null(prev_hypers) & trait_in != "EFD" ~ "code/jags-models/jags-briere-informative.bug",
      TPC_function == "Briere" & !is.null(prev_hypers) & trait_in == "EFD" ~ "code/jags-models/jags-briere-EFD-informative.bug",
      TPC_function == "Quadratic" & is.null(prev_hypers) ~ "code/jags-models/jags-quad-neg.bug",
      TPC_function == "Quadratic" & !is.null(prev_hypers) ~ "code/jags-models/jags-quad-neg-informative.bug",
      TPC_function == "AltQuadratic" & is.null(prev_hypers) ~ "code/jags-models/jags-quad-alt.bug",
      TPC_function == "AltQuadratic" & !is.null(prev_hypers) ~ "code/jags-models/jags-quad-alt-informative.bug",
      TPC_function == "Linear" & is.null(prev_hypers) ~ "code/jags-models/jags-linear.bug",
      TPC_function == "Linear" & !is.null(prev_hypers) ~ "code/jags-models/jags-linear-informative.bug"
    )
  }
  
  
  if (is.null(prev_hypers)) {
    jags_data <- list(
      "Y" = proper_studies$trait, "T" = proper_studies$T,
      "N" = length(proper_studies$T)
    )
  } else {
    # Following past work, reduce variance of hypers by a given factor
    # These mostly come from Shocket 2020 and Mordecai 2017
    hyper_relax_factor <- case_when(
      trait_in == "PDR" ~ 1,
      trait_in == "e2a" ~ 0.1,
      trait_in == "MDR" ~ 0.1,
      trait_in == "EFD" ~ 0.5,
      trait_in == "lf" ~ 0.5,
      trait_in == "c" ~ 0.5,
      trait_in == "b" ~ 0.5,
      trait_in == "a" ~ 0.5,
      trait_in == "EFOC" ~ 1,
      trait_in == "bc" ~ 0.1,
      trait_in == "EPR" ~ 1,
      trait_in == "pO" ~ 0.5,
      trait_in == "EV" ~ 0.1,
      trait_in == "pLA" ~ 0.1
    )
    jags_data <- list(
      "Y" = proper_studies$trait, "T" = proper_studies$T,
      "N" = length(proper_studies$T), "hypers" = prev_hypers * hyper_relax_factor
    )
  }
  
  # generate the samples from the jags.model
  samps <- run.jags(
    jags_data, TPC_function, variable_names,
    jags_choice, inits_list,
    n.chains, n.adapt, n.samps, prob_bool
  )
  
  # for the linear model (which doesn't have T0 as a parameter), add NAs)
  if (is.null(samps$T0)) {
    samps$T0 <- NA
  }
  
  samps <- mutate(samps, func = as.character(TPC_function))
  out.samps <- dplyr::select(samps, -c(sigma, tau))
  
  return(out.samps)
}


# Function: get hyperparameters of informed TPC parameter prior distributions
get.prior_hyperparams <- function(in_data, TPC_function, variable_names,
                                  jags_choice, inits_list,
                                  n.chains, n.adapt,
                                  scale_factor = 100, prob_bool) {
  # Initialize the hyperparameter list
  hypers <- NULL
  
  # get samples from trait TPC parameter distributions
  temp_samples <- run.jags(
    in_data, TPC_function, variable_names,
    jags_choice, inits_list,
    n.chains, n.adapt, n.samps, prob_bool
  )
  
  # transform TPC parameters for the Linear TPC
  if (TPC_function == "Linear") {
    temp_samples <- mutate(temp_samples, z = c * Tm) %>%
      dplyr::select(-c(Tm, T0))
  }
  
  # rescale TPC parameters to ease fitting
  temp_samples <- temp_samples / scale_factor
  temp_samples <- dplyr::select(temp_samples, -sample_num)
  
  # calculate hyperparameters of gamma distribution fit to TPC parameter distributions
  hypers <- apply(temp_samples, 2, function(df) {
    fitdistr(df, "gamma",
             method = "Nelder-Mead",
             hessian = FALSE
    )$estimate
  })
  
  # return to original scale
  hypers["rate", ] <- hypers["rate", ] / scale_factor
  
  return(hypers)
}

# Function: get samples of trait TPC hyperparameters from running jags
run.jags <- function(jags_data, TPC_function, variable_names,
                     jags_choice, inits_list,
                     n.chains, n.adapt, n.samps, prob_bool = FALSE) {
  # initialize jags model
  print("Initializing JAGS model...")
  jags <- jags.model(jags_choice,
                     data = jags_data,
                     n.chains = n.chains, inits = inits_list,
                     n.adapt = n.adapt,
                     quiet = FALSE # switch to FALSE to show messages and progress bars
  )
  # get n.samps new samples from the trait TPC parameter posterior distributions
  print("Sampling...")
  coda.samps <- coda.samples(jags, variable_names, n.samps)
  
  # put coda.samples into the format used for further analyses.
  if (TPC_function == "Briere") {
    samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
  } else if (TPC_function %in% c("Quadratic", "AltQuadratic")) {
    samps <- make.quad.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE) %>%
      mutate(c = qd, .keep = "unused")
  } else {
    samps <- make.linear.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE) %>%
      mutate(Tm = n.inter / slope) %>%
      mutate(c = slope) %>%
      mutate(T0 = NA) %>%
      dplyr::select(c, T0, Tm, sigma)
  }
  
  # Remove samples where the minimum temperature exceeds the maximum
  if (TPC_function != "AltQuadratic") {
    samps <- dplyr::filter(samps, is.na(T0) | Tm > T0)
    }
  
  # extra processing if the trait is a probability
  if (prob_bool) {
    quad_max_func <- function(c, T0, Tm) {
      c * ((Tm - T0) / 2)^2
    }
    briere_max_func <- function(c, T0, Tm) {
      T_star <- (1 / 10) * ((4 * Tm + 3 * T0) + sqrt((4 * Tm + 3 * T0)^2 - 40 * T0 * Tm))
      maxB <- c * T_star * (T_star - T0) * sqrt(Tm - T_star)
    }
    linear_max_func <- function(c, T0, Tm) {
      c * Tm
    }
    # remove samples which lead to maximum values exceeding one
    samps <- samps %>%
      mutate(test = case_when(
        TPC_function == "Quadratic" ~ quad_max_func(c, T0, Tm),
        TPC_function == "Briere" ~ briere_max_func(c, T0, Tm),
      )) %>%
      dplyr::filter(between(test, 0, 1)) %>%
      dplyr::select(-test)
  }
  
  # Resample until we obtain enough real samples
  while (dim(samps)[1] < n.chains * n.samps) {
    print("Resampling to ensure T0 < Tm ...")
    coda.samps <- coda.samples(jags, variable_names, n.samps)
    if (TPC_function == "Briere") {
      temp_samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
    } else if (TPC_function == "Quadratic") {
      temp_samps <- make.quad.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE) %>%
        mutate(c = qd, .keep = "unused")
    } else if (TPC_function == "Linear") {
      temp_samps <- make.linear.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps), sig = TRUE) %>%
        mutate(Tm = n.inter / slope) %>%
        mutate(c = slope) %>%
        dplyr::select(c, Tm, sigma)
    }
    temp_samps <- dplyr::filter(temp_samps, is.na(T0) | Tm > T0) %>%
      dplyr::filter(!is.na(Tm))
    
    # Resample if the max value for a probability trait exceeds 1
    if (prob_bool) {
      print("...and maximum does not exceed 1 for a probability")
      temp_samps <- temp_samps %>%
        mutate(test = case_when(
          TPC_function == "Quadratic" ~ quad_max_func(c, T0, Tm),
          TPC_function == "Briere" ~ briere_max_func(c, T0, Tm),
        )) %>%
        dplyr::filter(between(test, 0, 1)) %>%
        dplyr::select(-test) %>%
        dplyr::filter(!is.na(Tm))
    }
    
    print(paste0(
      "***Added an additional ",
      ifelse(dim(samps)[1] + dim(temp_samps)[1] > n.chains * n.samps,
             n.chains * n.samps - dim(samps)[1],
             dim(temp_samps)[1]
      ),
      " samples from re-sampling***"
    ))
    
    samps <- rbind(samps, temp_samps)
    print(paste0(min(100, 100 * round(dim(samps)[1] / (n.chains * n.samps), 2)), "% complete with resampling"))
  }
  samps <- samps[1:(n.chains * n.samps), ]
  
  samps$tau <- 1 / samps$sigma
  if (is.null(samps$c)) {
    samps <- mutate(samps, c = qd, .keep = "unused")
  }
  samps$sample_num <- 1:dim(samps)[1]
  return(samps)
}

# 2) Calculate prior distributions of thermal trait parameters from data ----

# Identify all distinct combinations of traits and transmission systems
data_in <- data.in.TPC %>% 
  # We don't need to keep track of where we got the data for anymore
  dplyr::select(-c(lead_author, year)) %>% 
  distinct()

distinct_combos <- data_in %>%
  filter(system_ID %in% c(
    "Aedes aegypti / DENV", "Aedes aegypti / none",
    "Aedes aegypti / ZIKV", "Aedes aegypti / none",
    "Aedes albopictus / DENV", "Aedes albopictus / none",
    "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
    "Culex spp. / WNV",
    "Anopheles gambiae / Plasmodium falciparum",
    "Anopheles gambiae / none"
  )) %>%
  distinct(trait.name, system_ID) %>%
  # Remove unused Culex spp. / WNV data (we only need b, c, and bc)
  filter(system_ID != "Culex spp. / WNV" | trait.name != "MDR") %>%
  filter(system_ID != "Culex spp. / WNV" | trait.name != "PDR")

samples <- tibble(
  trait = as.character(),
  system_ID = as.character(),
  T0 = as.double(),
  Tm = as.double(),
  c = as.double(), # note that we're using c as a generic parameter for Briere or Quadratic
  sample_num = as.double(),
  func = as.character()
)

# Go through all trait/system combinations to generate TPC parameter posterior samples
for (system_index in 1:dim(distinct_combos)[1]) {
  # Pull system information
  system_sample <- distinct_combos$system_ID[system_index]
  trait_in <- distinct_combos$trait.name[system_index]
  mosquito_in <- filter(data_in, system_ID == system_sample) %>%
    dplyr::select(mosquito_species) %>%
    unique() %>%
    as.character()
  pathogen_in <- filter(data_in, system_ID == system_sample, trait.name == trait_in) %>%
    dplyr::select(pathogen) %>%
    unique() %>%
    as.character()
  
  # Give a progress report
  print(paste0(
    "System # ", system_index, " of ", dim(distinct_combos)[1], ": ", mosquito_in, " / ", pathogen_in, " / ", trait_in,
    " -------------------------------------------------------------------------------"
  ))
  
  
  # generate TPC parameter posterior samples
  temp_sample <- thermtrait.prior.sample(data_in, trait_in, mosquito_in, pathogen_in,
                                         n.chains, n.adapt, n.samps,
                                         old_informative = FALSE
  )
  
  # add in name of system and trait to list of samples
  temp_sample <- temp_sample %>%
    mutate(
      trait = trait_in,
      system_ID = system_sample
    )
  # add samples to running list
  samples <- rbind(samples, temp_sample)
}
# )

data.in.transform <- samples

# check that we have all the data we need
distinct_samples <- data.in.transform %>%
  distinct(system_ID, trait)

print(distinct_samples)


# *) Diagnostics & visualizations -----------------------------------------

# Do you just want to look at focal species?
focal_bool <- TRUE

if (plot_bool) {
  plot_samples <- data.in.transform %>%
    filter(sample_num %in% 1:1000)
  
  if (focal_bool) {
    plot_samples <- filter(plot_samples, system_ID %in% c(
      "Aedes aegypti / DENV", "Aedes aegypti / none",
      "Aedes aegypti / ZIKV", "Aedes aegypti / none",
      "Aedes albopictus / DENV", "Aedes albopictus / none",
      "Culex quinquefasciatus / WNV", "Culex quinquefasciatus / none",
      "Culex univittatus / WNV",
      "Anopheles spp. / Plasmodium spp.",
      "Anopheles spp. / none"
    ))
  }
  library(reshape2)
  library(cowplot)
  # Define TPC functions
  # Briere function
  Briere <- function(q, Tmin, Tmax) {
    function(t) {
      pmax(q * t * (t - Tmin) * (Tmax - t)^(1 / 2), 0, na.rm = TRUE)
    }
  }
  # Quadratic function
  Quadratic <- function(q, Tmin, Tmax) {
    function(t) {
      pmax(-q * (t - Tmin) * (t - Tmax), 0, na.rm = TRUE)
    }
  }
  # Linear
  Linear <- function(q, Tmax) {
    function(t) {
      pmax(q * (Tmax - t), 0, na.RM = FALSE)
    }
  }
  # Function to thin intervals for samples
  res_reduce <- function(df, new_length) {
    old_length <- length(unique(df))
    if (old_length < new_length) {
      warning("new_length is larger than old length. Vector will be unchanged")
      new_length <- old_length
    }
    ret <- unique(df)[seq(1, length(unique(df)), length.out = new_length)]
  }
  
  # Figure: Histograms of parameter distributions of thermal traits ----
  # distributions should be clumped (except for c)
  # logc should be clumped
  # temp_diff should be positive
  plot_df <- plot_samples %>%
    dplyr::select(-func) %>%
    mutate(temp_diff = Tm - T0) %>%
    mutate(logc = log(c)) %>%
    melt(id = c("system_ID", "trait", "sample_num"))
  
  label_order =c("Gonotrophic cycle rate",
                 "Eggs per female per day",
                 "Eggs per female per oviposition cycle",
                 "\\% females ovipositing",
                 "Eggs per raft",
                 "\\% of egg rafts that hatch",
                 "\\# of emerging larvae per raft",
                 "Egg viability",
                 "Pr(egg -> adult survival)",
                 "Pr(larva -> adult survival)",
                 "Mosquito development rate",
                 "Lifespan",
                 "Vector competence",
                 "Pr(infectious | infected)",
                 "Pr(infected | exposure)",
                 "Parasite development rate")
  
  plot_df$trait_label <-  case_match(
    plot_df$trait,
    "a" ~ "Gonotrophic cycle rate",
    "EFD" ~ "Eggs per female per day",
    "EFOC" ~ "Eggs per female per oviposition cycle",
    "pO" ~ "\\% females ovipositing",
    "EPR" ~ "Eggs per raft",
    "pRH" ~ "\\% of egg rafts that hatch",
    "nLR" ~ "\\# of emerging larvae per raft",
    "EV" ~ "Egg viability",
    "e2a" ~ "Pr(egg -> adult survival)",
    "pLA" ~ "Pr(larva -> adult survival)",
    "MDR" ~ "Mosquito development rate",
    "lf" ~ "Lifespan",
    "bc" ~ "Vector competence",
    "c" ~ "Pr(infectious | infected)",
    "b" ~ "Pr(infected | exposure)",
    "PDR" ~ "Parasite development rate"
  )
  
  plot_df <- plot_df %>% mutate(trait_label = factor(trait_label, levels = label_order))
  
  parm_hists <- plot_df %>%
    filter(variable %in%  c("T0", "Tm", "logc")) %>%
    mutate(var_label = case_when(
      variable == "T0" ~  "CT_min",
      variable == "Tm" ~  "CT_max",
      variable == "logc" ~  "log(q)"
    )) %>% 
    ggplot(aes(value, color = system_ID, fill = system_ID)) +
    geom_histogram(aes(), bins = 100) +
    facet_grid(trait ~ var_label, scales = "free") +
    theme_minimal_grid(12)
  
  # # Save figure
  # ggsave("figures/imputed_traits/parm_hists.svg",
  #        device = "svg",
  #        width = 16, height = 9, units = "in"
  # )
  
  ###* Figure: TPC curves with 89% high confidence intervals ----
  
  # Temperature vector used for visualizations
  Temps <- seq(0, 50, length.out = 200)
  
  # thin_size <- 100
  
  # For each mosquito species, trait, and sample, get a thermal response curve
  TPC_df <- plot_samples %>%
    # filter(sample_num %in% seq(1, thin_size)) %>%
    full_join(list(Temperature = Temps), by = character(), copy = TRUE) %>%
    mutate(Trait_val = case_when(
      func == "Briere" ~ Briere(c, T0, Tm)(Temperature),
      func == "Quadratic" ~ Quadratic(c, T0, Tm)(Temperature),
      func == "AltQuadratic" ~ AltQuadratic(c, T0, Tm)(Temperature),
      func == "Linear" ~ Linear(c, Tm)(Temperature)
    )) %>%
    dplyr::select(-c("c", "T0", "Tm"))
  
  # get mean TPC from samples
  meanTPC_df <- TPC_df %>%
    group_by(system_ID, trait, Temperature) %>%
    summarise(mean_val = mean(Trait_val), .groups = "keep") %>%
    unique()
  
  # get edges of 89% HCI of samples
  quantsTPC_df <- TPC_df %>%
    group_by(system_ID, trait, Temperature) %>%
    mutate(lowHCI_val = quantile(Trait_val, 0.055)) %>%
    mutate(highHCI_val = quantile(Trait_val, 0.945)) %>%
    dplyr::select(-c("sample_num", "Trait_val", "func")) %>%
    unique()
  
  TPC_plot_all <- TPC_df %>%
    filter(sample_num < 51) %>%
    group_by(sample_num) %>%
    arrange(Temperature, trait, system_ID) %>% 
    # group_by()
    ggplot(aes(x = Temperature, y = Trait_val, color = system_ID, group = sample_num)) +
    # means of TPC curves
    geom_line(alpha = 0.3) +
    ylab("") +
    facet_wrap(~trait, scales = "free", ncol = 2) +
    theme_minimal_grid(12)
  # Save figure
  ggsave("TPC_plot_all.svg", TPC_plot_all,
         device = "svg",
         width = 16, height = 9, units = "in"
  )
  
  TPC_plot <- TPC_df %>%
    group_by(sample_num) %>%
    arrange(Temperature) %>%
    # group_by()
    ggplot() +
    # means of TPC curves
    geom_line(
      data = meanTPC_df,
      aes(x = Temperature, y = mean_val, color = system_ID)
    ) +
    # 89% HCI of TPC curves
    geom_ribbon(
      data = quantsTPC_df,
      aes(x = Temperature, ymin = lowHCI_val, ymax = highHCI_val, fill = system_ID),
      alpha = 0.1
    ) +
    ylab("") +
    facet_wrap(~trait, scales = "free", ncol = 2) +
    theme_minimal_grid(12)
  
  # # Save figure
  # ggsave("figures/imputed_traits/TPC_plot.svg", TPC_plot,
  #        device = "svg",
  #        width = 16, height = 9, units = "in"
  # )
}