## Title: Transform TPC traits into model parameters #######################
##
## Project: Global zoonoses - spillover of mosquito-borne pathogens
##
## Purpose: Transform trait values, which we have fit TPCs to from data, in 
##          order to use in model analyses
##
## Contents: 0) Set-up, load in necessary packages and data-sets
##           1) Define accessory functions
##           2) Create TPCs and deal with missing data
##           3) Create parameters data frame
##           4) Save parameters data frame
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

# 1) Define accessory functions -------------------------------------------
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
# Alternate Quadratic function
AltQuadratic <- function(q, Tmin, Tmax) {
  function(t) {
    1/(q * t^2 - Tmin * t + Tmax)
  }
}

# Linear
Linear <- function(q, Tmax) {
  function(t) {
    pmax(q * (Tmax - t), 0, na.RM = FALSE)
  }
}

# Function: designate proper thermal response function
# - output is a function of temperature
get.thermal.response <- function(data_in, Temperature) {
  parms <- dplyr::select(data_in, c, T0, Tm)
  function_type <- dplyr::select(data_in, func)
  
  temp_function <- case_when(
    function_type == "Briere" ~ Briere(parms$c, parms$T0, parms$Tm),
    function_type == "Quadratic" ~ Quadratic(parms$c, parms$T0, parms$Tm),
    function_type == "AltQuadratic" ~ AltQuadratic(parms$c, parms$T0, parms$Tm),
    function_type == "Linear" ~ Linear(parms$c, parms$Tm)
  )
  
  out <- temp_function(Temperature)
}

# 2) Create TPCs and deal with missing data ----

# choose a random set of samples from the TPC parameters
num_samples <- length(unique(data.in.transform$sample_num))
sample_inds <- sample(1:num_samples, thin_size, replace = FALSE)

# Create data frame of TPCs
TPC_df <- data.in.transform %>%
  filter(sample_num %in% sample_inds) %>%
  cross_join(list(Temperature = Temps), copy = TRUE) %>%
  mutate(Trait_val = case_when(
    func == "Briere" ~ Briere(c, T0, Tm)(Temperature),
    func == "Quadratic" ~ Quadratic(c, T0, Tm)(Temperature),
    func == "AltQuadratic" ~ AltQuadratic(c, T0, Tm)(Temperature),
    func == "Linear" ~ Linear(c, Tm)(Temperature)
  )) %>%
  dplyr::select(-c("c", "T0", "Tm")) %>% 
  pivot_wider(id_cols = c("system_ID", "sample_num", "Temperature"), 
              names_from = "trait",
              values_from = "Trait_val")  %>% 
  # Combine traits into intermediate parameters as necessary
  # i.e. putting together reproductive traits to estimate eggs per female per day (EFD)
  mutate(e2a = case_when(
    !(is.na(e2a)) ~ e2a,
    !(is.na(EV * pLA)) ~ EV * pLA
  )) %>% 
  mutate(EFD = case_when(
    !(is.na(EFD)) ~ EFD,
    !(is.na(EFOC * a)) ~ EFOC * a,
    !(is.na(EPR * pO)) ~ EV * EPR * pO * a
  )) %>% 
  mutate(bc = ifelse(is.na(b) & is.na(c), bc, b * c)) %>% 
  # throw out traits that are no longer needed
  dplyr::select(system_ID, sample_num, Temperature,
                a, bc, PDR, e2a, EFD, lf, MDR) %>% 
  # separate out mosquito species and pathogen names
  separate_wider_delim(system_ID, delim = " / ", names = c("mosquito_species","pathogen"))

# Combine parasite relevant traits with mosquito life history traits according to system
noInfection_df <- TPC_df %>% 
  filter(pathogen == "none") %>% 
  dplyr::select(-pathogen) %>% 
  dplyr::select(-c("bc", "PDR"))

Infection_df <- TPC_df %>% 
  filter(pathogen != "none") %>% 
  dplyr::select(-c("a", "e2a", "EFD", "lf", "MDR"))

combined_df <- left_join(Infection_df, noInfection_df)%>% 
  unite(
    col = "system_ID",
    c("mosquito_species", "pathogen"),
    sep = " / ",
    remove = FALSE
  ) %>% 
  relocate(system_ID, mosquito_species, pathogen, Temperature, sample_num)

# Similar to Shocket 2020: use Culex spp. / WNV / bc data for Culex quinquefasciatus / WNV / bc
# This combines data for Culex univittatus, tarsalis, and pipiens
combined_df <- combined_df %>% 
  # temporarily remove Cx. quinquefasciatus / WNV rows
  filter(system_ID != "Culex quinquefasciatus / WNV") %>% 
  # add rows back in after switching in "bc" values from Cx. univittatus / WNV
  rbind(filter(combined_df, system_ID == "Culex quinquefasciatus / WNV") %>% 
          # remove original bc values (all NA)
          dplyr::select(-bc) %>% 
          # join with Cx. univittatus data
          right_join(filter(combined_df, system_ID == "Culex spp. / WNV") %>% 
                       dplyr::select(Temperature:bc))) %>% 
  # remove Cx. univittatus data
  filter(system_ID != "Culex spp. / WNV") 

## Carrying capacity for larval mosquitoes
# NB: In the absence of good estimates for each species or temperature-dependence of this trait, we assume that this parameter is constant. It can be used as a  scaling parameter for overall mosquito abundance
# (it could alternately be used to fix the maximum adult mosquito density across species)
larval_mosquito_carrying_capacity <- 300

# 3) Make parameter dataframe ---------------------------------------------
eps <- .Machine$double.eps

data.in.params <- combined_df %>%
  # Biting rate.
  mutate(sigmaV = a) %>% 
  # Fecundity. Simplified from a * 0.5 * EFD / a
  mutate(sigmaV_f = 0.5 * EFD) %>% # = female eggs per female per day
  # Aquatic-stage mosquito survival probability.
  mutate(deltaL = e2a) %>% 
  # Aquatic-stage mosquito development rate.
  mutate(rhoL = MDR) %>% 
  mutate(etaL = rhoL / (deltaL + eps)) %>%
  # Aquatic-stage mosquito mortality rate. This should be ~infinite if deltaL = 0
  mutate(muL = etaL - rhoL) %>%   
  # Adult mosquito average lifespan
  # Pathogen development rate.
  mutate(etaV = PDR) %>%
  # Mosquito infection probability.
  mutate(betaV = bc) %>% 
  dplyr::select(system_ID:sample_num, lf, sigmaV:betaV) %>%
  mutate(KL = larval_mosquito_carrying_capacity) %>%
  mutate(V0 = ifelse( sigmaV_f * deltaL < (1 / lf),
                      0,
                      KL * rhoL * lf * (1 - 1 / (lf * sigmaV_f * deltaL))))

# *) Diagnostics & visualizations -----------------------------------------

if (plot_bool) {
  
library(cowplot)
library(svglite)  
  
  # Helper function to place legends in empty facets of plot grids
  # Code by: Artem Sokolov, found here: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
  shift_legend <- function(p) {
    pnls <- cowplot::plot_to_gtable(p) %>%
      gtable::gtable_filter("panel") %>%
      with(setNames(grobs, layout$name)) %>%
      purrr::keep(~ identical(.x, zeroGrob()))
    
    if (length(pnls) == 0) {return(p)} #stop("No empty facets in the plot")
    
    lemon::reposition_legend(p, "center", panel = names(pnls))
  }
  
# For each mosquito species, trait, and sample, get a thermal response curve
TPC_df <- data.in.params %>%  
  ungroup() %>% 
  dplyr::select(-c(muL, etaL, mosquito_species, pathogen)) %>%
  melt(id = c("system_ID", "Temperature", "sample_num"),
       variable.name = "trait",
       value.name = "Trait_val") %>% 
  group_by(system_ID, trait, Temperature) %>%
  summarise(mean = mean(Trait_val),
            median = median(Trait_val),
            lowHCI = quantile(Trait_val, 0.05, na.rm = TRUE),
            highHCI = quantile(Trait_val, 0.95, na.rm = TRUE),
            .groups = "keep") %>%
  unique() %>% ungroup()

TPC_plot <- TPC_df %>%
  arrange(Temperature) %>%
  ggplot() +
  # means of TPC curves
  geom_path(aes(x = Temperature, y = mean, color = system_ID)) +
  geom_ribbon(
    aes(x = Temperature, ymin = lowHCI, ymax = highHCI, fill = system_ID),
    alpha = 0.1
  ) +
  ylab("") +
  facet_wrap(~trait, scales = "free", ncol = 2) +
  theme_minimal_grid(12)

TPC_plot <- shift_legend(TPC_plot)

# Save figure
ggsave("figures/imputed_traits/param_TPC_plot.svg",
       plot = TPC_plot,
       device = "svg",
       width = 16, height = 9, units = "in")
}