Mordecai-2017 folder
This folder provides copies of scripts originally provided with the publication:
Mordecai, E. A., Cohen, J. M., Evans, M. V., Gudapati, P., Johnson, L. R., Lippi, C. A., Miazgowicz, K., Murdock, C. C., Rohr, J. R., Ryan, S. J., Savage, V., Shocket, M. S., Stewart Ibarra, A., Thomas, M. B., & Weikel, D. P. (2017). Detecting the impact of temperature on transmission of Zika, dengue, and chikungunya using mechanistic models. PLoS Neglected Tropical Diseases, 11(4), e0005568. https://doi.org/10.1371/journal.pntd.0005568

The code is also made available online at the URL https://figshare.com/s/b79bc7537201e7b5603f


collect_references.R
  * This script collects the all the data sources used to parameterize mosquito traits in this study and organizes them into a table for easy viewing.

# Code for performing analyses
host-traits.R
  * Sets up vertebrate host parameters

vector-traits.R
  * Sets up vector trait parameters, these are functions of temperature and are sampled from posterior distributions constructed in get-thermal-trait-priors.R

  - data-cleaning.R
    * Processes data from data/raw/ folder into a format used to derive mosquito trait TPC posterior distributions in get-thermal-trait-priors.R

  - get-thermal-trait-priors.R
    * Samples from data-informed mosquito trait TPC hyperparameters using models in code/jags-models folder and functions from code/Mordecai_2017

  - trait-transform.R
    * Transform traits into parameters to be input into the model

get-outputs.R
  * Calculate the basic reproduction number (R0), transmission thermal optimum (Topt), critical thermal minimum (CTmin), critical thermal maximum (CTmax), and width of the parasite population thermal niche (CTwidth)

sensitivity.R
  * Sensitivity and uncertainty analysis of R0, Topt, CTmin, CTmax, and CTwidth
