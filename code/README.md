# Code for performing all analyses
## 'Mordecai-2017' folder
This folder provides copies of scripts originally provided with the publication:
Mordecai, E. A., Cohen, J. M., Evans, M. V., Gudapati, P., Johnson, L. R., Lippi, C. A., Miazgowicz, K., Murdock, C. C., Rohr, J. R., Ryan, S. J., Savage, V., Shocket, M. S., Stewart Ibarra, A., Thomas, M. B., & Weikel, D. P. (2017). Detecting the impact of temperature on transmission of Zika, dengue, and chikungunya using mechanistic models. PLoS Neglected Tropical Diseases, 11(4), e0005568. https://doi.org/10.1371/journal.pntd.0005568

The code is also made available online at the URL https://figshare.com/s/b79bc7537201e7b5603f

## 'jags-models' folder
This folder collects all the jags models used in the creation of trait thermal performance curves. These files are modified from those originally provided in Mordecai et al. 

## Original code
To re-create all study analyses, first fun 'vector-traits.R' to obtain the trait thermal performance curves and all vector parameter values, then 'host-traits.R' to obtain all vertebrate host parameter values. Next, 'get-outputs.R' will calculate all the outputs considered in the study and 'sensitivity.R' provides the sensitivity analysis results. Note that some scripts are particularly time and memory intensive. Short descriptions of each R script are provided below.

collect_references.R
  * This script collects the all the data sources used to parameterize mosquito traits in this study and organizes them into a table for easy viewing.

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
