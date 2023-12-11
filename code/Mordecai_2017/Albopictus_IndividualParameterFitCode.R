### This file will process the temperature data and create the individual parameter 
### samples for all the traits that are included in the Aedes albopictus R0 model. 

### This code fits models using UNINFORMATIVE priors

## Remember to set your working directory so the code can access the data and 
## necessary supplementary code.

## Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')

# This file contains tools for analysis and visualization.
source("mcmc_utils_all.R") 

# This file contains the thermal response functions and their derivatives
source("temp_functions_all.R") 

# Creating a small constant to keep denominators from being zero.
ec<-0.000001 

# Specifing the parameters that control the MCMC (these will be used throughout the code). 

n.chains <- 5
n.adapt <- 5000
n.samps <- 5000

## Loading the data; my data will be called albopictusCHIKVmodelTempData.csv 
data.all <- read.csv("albopictusCHIKVmodelTempData_2016-03-26.csv", header=TRUE)

## Now the code will choose all the temperature sensitive traits and fit either a
## Briere or Quadratic model to it using MCMC sampling.

################################
## The first trait the biting rate, a.

# This chooses the trait, GCD, which is used to model the biting rate, a.

data <- data.all[which(data.all$trait.name=="GCD"),]

# Transforming the GCD into a, the biting rate

data$trait = (1/(data$trait))

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.
par(mfrow=c(1,1))
plot(trait ~ T, data = data)

# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere.bug',
                    data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                    n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                    n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
a.samps <- samps

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
priors1<-list()
priors1$names<-c( "T0", "Tm", "c","tau")
priors1$fun<-c( "uniform", "uniform", "gamma","gamma")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(0, 24, NA)
priors1$hyper[,2]<-c(25, 45, NA)
priors1$hyper[,3]<-c(1, 10, NA)
priors1$hyper[,4]<-c(0.0001, 0.0001, NA)
plot.hists(a.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", a.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(10, 40), ylim=c(0,0.4),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

################################################
## The next trait that is fitted is MDR, the mean development time for a mosquito. 

data1 <- data.all[which(data.all$trait.name=="MDR"),]
data2 <- data.all[which(data.all$trait.name=="1/MDR"),]
data2$trait <- 1/data2$trait
data = rbind(data1, data2)

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)
points(trait ~ T, data=subset(data, ref=="Ezeakacha_Dissertation_2015"), col=2, pch=16)
points(trait ~ T, data = subset(data, ref=="Yee_et al_2016 JAE_in_review"), col = 3, pch = 16)

# Given the data the Briere fuction is chosen. Jag-briere.bug contains the specifics of 
# the Briere model with the default priors.

jags <- jags.model('jags-briere.bug',
                    data = list('Y' = data$trait, 'T' = data$T, 'N' = length(data$T)),
                    n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                    n.adapt = n.adapt) 

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c','Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnostic information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps <- make.briere.samps(coda.samps, nchains=n.chains, samp.lims=c(1, n.samps))
samps$tau <- 1/samps$sigma
MDR.samps <-  samps

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
plot.hists(MDR.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", MDR.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(10, 40),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

#################################
## The next trait we will be modeling will be TFD, a measure of total female fecundity,
## which will be used in the R0 model as a place holder for the EFD, another fecundity
## measuere, (EFD = TFD*(1/a), where a is the biting rate). 

# Here we select TFD out of the data, and furthermore select the first gonotrophic 
# cycle due to some issues with the way the data were collected.

data <- data.all[which(data.all$trait.name=="TFD"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait~T, data=data)
points(trait ~ T, data = subset(data, trait2.name=="R1"), pch=16, col=2)
points(trait ~ T, data = subset(data, trait2.name=="R2"), pch=16, col=3)
points(trait ~ T, data = subset(data, trait2.name=="R3"), pch=16, col=4)
points(trait ~ T, data = subset(data, trait2.name=="R4"), pch=16, col=5)
points(trait ~ T, data = subset(data, trait2.name=="R5"), pch=16, col=6)
points(trait ~ T, data = subset(data, trait2.name=="R6"), pch=16, col=7)
points(trait ~ T, data = subset(data, trait2.name=="R7"), pch=16, col=8)
points(trait ~ T, data = subset(data, trait2.name=="R8"), pch=16, col=1)

# R indicates the oviposition cycle number
# omit R4, R5, R6, R7, R8
data <- data[which(data$trait2.name %in% c("R1", "R2", "R3", NA)),]
data = subset(data, ref != "Yee_et al_2016 JAE_in_review")
plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data,ref == "Ezeakacha_Dissertation_2015"), col = 2, pch=16)

# Given the data the Briere fuction is chosen. Jag-briere.bug contains the specifics of 
# the Briere model with the default priors.

jags <- jags.model('jags-briere.bug',
                    data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T)),
                    n.chains = n.chains, inits=list(Tm=31, T0=5, c=0.00007),
                    n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis

coda.samps <- coda.samples(jags, c('c','Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses.

samps <- make.briere.samps(coda.samps, nchains=n.chains, samp.lims=c(1, n.samps))
samps$tau <- 1/samps$sigma
TFD.samps <- samps

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
plot.hists(TFD.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", TFD.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(10, 40),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

###################################
## The next parameter is pEA, the probability a mosquito will survive from hatching to 
## maturation. 

data <- data.all[which(data.all$trait.name=="pEA"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait~T, data=data)

# Given the data the Negative Quadratic function is chosen. Jags-quad-neg.bug contains
# the specifics of the Negative Quadratic model with the default priors. 

jags <- jags.model('jags-quad-neg.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T)),
                   n.chains = n.chains,
                   inits=list(T0=5, Tm=33, n.qd=0.005), n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis

coda.samps <- coda.samples(jags, c('T0','Tm', 'qd', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 
plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps.q <- make.quad.samps(coda.samps, nchains=n.chains, 
						   samp.lims=c(1, n.samps), sig=TRUE)
samps.q$n.qd <- samps.q$qd
samps.q$tau <- 1/samps.q$sigma
e2a.samps <- samps.q

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
plot.hists(e2a.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="quad.trunc", e2a.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(5, 40), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

#########################################
# The next trait that will be fit is p, the survival probability of an adult mosquito. 
# This is fit as the longevity, lf, in lifespan-comp-script.R 

#########################################
# b, transmission probability
data <- data.all[ which(data.all$trait.name=="b"),]
# the only usable data are Xiao et al. 2014
# DENV-2 in Ae. albopictus
data = subset(data, ref=="Xiao_et_al_2014_Arch Virol")
plot(trait ~ T, data = data)

# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
b.samps <- samps

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
plot.hists(b.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", b.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(5, 40), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)


#########################################
# c, infection probability
data <- data.all[ which(data.all$trait.name=="c"),]
# the only usable data are Xiao et al. 2014
data = subset(data, ref=="Xiao_et_al_2014_Arch Virol")
plot(trait ~ T, data = data)

# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
c.samps <- samps

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
plot.hists(c.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", c.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(10, 40), ylim = c(0, 1), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)


#########################################
# PDR, parasite development rate = 1/EIP
data <- data.all[which(data.all$trait.name=="EIP"),]
data$trait <- 1/data$trait
plot(trait ~ T, data = data)

# Given the data the Briere fuction is chosen. Jag-briere.bug contains the specifics of 
# the Briere model with the default priors.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N' = length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 38, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt) 

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c','Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps <- make.briere.samps(coda.samps, nchains=n.chains, samp.lims=c(1, n.samps))
samps$tau <- 1/samps$sigma
PDR.samps <-  samps

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
plot.hists(PDR.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", PDR.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(10, 45), ylim = c(0, 0.32),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)


## This code is just save the MCMC samples for further analysis in the R0 model.
save(a.samps, MDR.samps, TFD.samps, e2a.samps, b.samps, c.samps, PDR.samps,
	 file = "Albopictus_ParameterFits_2016-03-30.Rsave")

