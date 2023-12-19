## Title: Vertebrate host traits and parameters #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Record and save host traits for use in analyses
##
## Written and maintained by: Kyle Dahlin, kydahlin@gmail.com
## Initialized March 2023
# _______________________________________________________________________________

library(tidyverse)

## Set resolution for host trait variation
# Host density vector: Number of values to include to consider for vertebrate host density
KH_vec_length <- 30 #300

# Biting tolerance vector: Number of values to consider for biting tolerance
sigmaH_vec_length <- 30 #300

## Host life history & behavioral traits
# Host recruitment rate:
# upper estimate from range for Primate traits: (0.001150685, 0.009624300)
lambdaH_baseline <- .005

# Host mortality rate:
# estimate from range for Primate traits: lifespan (8.6, 60) years
muH_baseline <- 1 / (365 * 20)

# Host maximum biting tolerance (mosquitoes bites per day)
sigmaH_vec <- 10^seq(-0.25, 3.25, length.out = sigmaH_vec_length - 7) %>%
  c(1, 10, 20, 50, 100, 1000, Inf) %>%
  unique() %>% sort()
sigmaH_baseline <- 100

# Host carrying capacity
KH_vec <- 10^seq(-2, 5, length.out = KH_vec_length - 6) %>%
  c(10^seq(-2,5)) %>%
  unique() %>% sort()

## Host-related pathogen parameters
# Probability of becoming infected after bite from infectious mosquito
# operates as a scaling parameter
betaH_baseline <- 1

# Host recovery rate
# plausible estimate for infectious period = (5-14 days)
gammaH_vec <- 1 / c(5, 14)
gammaH_baseline <- 1 / 5

data.Host <- expand_grid(
  # Life-history parameters
  lambdaH = lambdaH_baseline,
  muH = muH_baseline,
  KH = KH_vec,
  sigmaH = sigmaH_vec,
  # Infection-related parameters
  gammaH = gammaH_baseline,
  betaH = betaH_baseline) %>%
  as.data.frame() %>%
  mutate(Model = ifelse(is.infinite(sigmaH), "Ross-Macdonald model", "Chitnis model"))

# Save host traits
write_rds(data.Host, "results/Host_vals.rds")