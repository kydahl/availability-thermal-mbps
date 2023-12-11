#####
# ReadMe.txt
This document describes the code and files for the paper ‘Detecting the impact of temperature on transmission of Zika, dengue and chikungunya using mechanistic models.’ In all analyses, the informative prior models are the default. The scripts for key analyses are described below. Remaining R scripts are for data pre-processing and functions called within the other scripts.

#####
# Thermal response fitting from data:
Aegypti_IndividualParameterFitCode.R
- Fits trait thermal responses for Ae. aegypti and DENV from data using uninformative priors.
- Informative_Aegypti_IndividualParameterFitCode.R fits thermal responses using the same data and informative priors fit in Aedes_prior_fits.R

Albopictus_IndividualParameterFitCode.R
- Fits trait thermal responses for Ae. albopictus from data using uninformative priors. 
- Informative_Albopictus_IndividualParameterFitCode.R fits thermal responses using the same data and informative priors fit in Aedes_prior_fits.R

lifespan-comp-script.R
- Fits lifespan thermal responses for both Aegypti and Albopictus using uninformative priors.
- lifespan-comp-script-informative.R does the same using informative priors fit in 

Aedes_prior_fits.R
- Fits the informative priors for the thermal response functions

#####
# Calculating R0 from the thermal responses:
Aegypti_DENV_R0_Analysis.R
- Calculates R0 vs. temperature for all posterior samples for the Aegypti model.
- Informative_Aegypti_DENV_R0_Analysis.R does the same for the informative prior model.

Albopictus_R0_Analysis.R
- Calculates R0 vs. temperature for all posterior samples for the Albopictus model.
- Informative_Albopictus_R0_Analysis.R does the same for the informative prior model.

#####
# Plotting R0 results
R0_plotting.R
- Plots R0 vs. temperature for the Aegypti and Albopictus models, and various other model outputs.

##### 
# Validation data acquisition and analyses
inc_R0_validation_dataprep.R
- Validation data prep and clean-up

validation3.Rmd
- Runs all statistical model validation analyses

#####
# Sensitivity and uncertainty analyses
Aegypti_DENV_R0_Sensitivity.R
- Calculates and plots the derivatives of R0 with respect to each trait for the Aegypti model.

Albopictus_R0_Sensitivity.R
- Calculates and plots the derivatives of R0 with respect to each trait for the Albopictus model.

Aegypti_R0_bc_EIP_Sensitivity.R
- Calculates R0 in the Aegypti models using a number of different assumptions for vector competence and EIP to test its sensitivity.

Albopictus_R0_bc_EIP_Sensitivity.R
- Calculates R0 in the Albopictus models using a number of different assumptions for vector competence and EIP to test its sensitivity.

