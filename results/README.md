# Results from scripts in '/code'
This folder is populated after running the scripts provided in '/code'
These data are compressed as .rds files and may become substantially larger when loaded into the R workspace.

VecTPC_vals.rds
  * Data specifying the mosquito and parasite trait thermal performance curves
  * Created by vector-traits.R

R0_TPC_vals.rds
  * Values of the basic reproduction number R0 as a function of temperature
  * Created by get-outputs.R

R0_vals.rds
  * Values of the basic reproduction number R0 evaluated across a high resolution range of temperature and host availability values
  * Created by get-outputs.R

Topt_vals.rds
  * Values of the thermal optimum evaluated across a high resolution range of host availability values
  * Created by get-outputs.R

CT_vals.rds
  * Values of the critical thermal minimum, maximum, and range width evaluated across a high resolution range of host availability values
  * Created by get-outputs.R

dR0dk_vals.rds
  * Values of the local derivative of R0 with respect to temperature-dependent mosquito parameters
  * Created by get-outputs.R

CTminmaxTopt_densities.rds
  * Distributions of the thermal optimum, critical thermal minimum, critical thermal maximum derived from the posterior distributions of the trait thermal performance curve hyperparameters.
  * Created by sensitivity.R

Topt_HPD_unc.rds
  * Uncertainty in the thermal optimum measured using the widths of the highest posterior density interval
  * Created by sensitivity.R

CT_HPD_unc.rds.rds
  * Uncertainty in the critical thermal minimum and maximum measured using the widths of the highest posterior density interval
  * Created by sensitivity.R