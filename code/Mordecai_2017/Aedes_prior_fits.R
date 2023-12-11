# Fitting models to the prior data

# Load saved versions to save time
# load("aedes_prior_samps.Rsave")
# load("aedes_prior_gamma_fits.Rsave")

## Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')
library('MASS')
library('here')

setwd("C:/Users/kd99491/Documents/GitHub/thermal-properties-mbps/code/Mordecai2017")

# This file contains tools for analysis and visualization.
source("mcmc_utils_all.R") 

# This file contains the thermal response functions and their derivatives.
source("temp_functions_all.R")  

# Creating a small constant to keep denominators from being zero.
ec<-0.000001 

## Loading the data; my data will be called aegyptiDENVmodelTempData.csv 
data.all <- read.csv("Aedes_prior_data.csv", header=TRUE)

## Load the default priors
priors1<-list()
priors1$names<-c( "T0", "Tm", "c","tau")
priors1$fun<-c( "uniform", "uniform", "gamma","gamma")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(0, 24, NA)
priors1$hyper[,2]<-c(25, 45, NA)
priors1$hyper[,3]<-c(1, 10, NA)
priors1$hyper[,4]<-c(0.0001, 0.0001, NA)

# Specifing the parameters that control the MCMC (these will be used throughout the code). 

n.chains <- 5
n.adapt <- 5000
n.samps <- 5000


######################################
## The first trait the biting rate, a.

data <- data.all[ which(data.all$trait.name=="a"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data, ref=="Callado (2002)"), col=2, pch=16)


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

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", a.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(15, 40), ylim = c(0, 0.5),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

# Use moment matching to fit a gamma distribution to each parameter estimate
gamma.fits.a = apply(a.samps, 2, function(df) fitdistr(df, "gamma")$estimate)
# warnings messages are OK

# Plot an example to show that it's working
hist(a.samps$T0)
y=dgamma(Temps, gamma.fits.a[1,1], gamma.fits.a[2,1])
lines(Temps, 10000*y/max(y))

#############################################
## The second trait to analyze is EFD, or the number of eggs laid
## per female per day.

data <- data.all[which(data.all$trait.name=="EFD"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data, ref=="Calado&Navarro-Silva_2002"), pch=16, col=2)

# The two data sources have different magnitudes but similar curves
# Just use the Joshi data since it has more data

data = subset(data, ref=="Joshi_1996")
plot(trait ~ T, data = data)

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
EFD.samps <- samps

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", EFD.samps, Temps, thinned=seq(1,n.samps, length=1000))


## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(10, 40), ylim = c(0, 18),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

# Fit a gamma distribution to each parameter
gamma.fits.EFD = apply(EFD.samps, 2, function(df) fitdistr(df, "gamma")$estimate)

################################
### Create a prior on TFD by rescaling EFD

data <- data.all[which(data.all$trait.name=="EFD"),]
# rescale EFD to have the same max value as the Aedes albopictus TFD data
data$trait <- data$trait*(77.19048/max(data$trait))

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data, ref=="Calado&Navarro-Silva_2002"), pch=16, col=2)

# The two data sources have different magnitudes but similar curves
# Just use the Joshi data since it has more data

data = subset(data, ref=="Joshi_1996")
plot(trait ~ T, data = data)

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
TFD.samps <- samps

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

# Fit a gamma distribution to each parameter
gamma.fits.TFD = apply(TFD.samps, 2, function(df) fitdistr(df, "gamma")$estimate)


########################################
## The next trait that is fitted is MDR, the mean development time for a mosquito. 

data1 <- data.all[ which(data.all$trait.name=="mdr"),]
data2 <- data.all[ which(data.all$trait.name=="1/MDR"),]
data2$trait <- 1/data2$trait
data = rbind(data1, data2)

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

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

# These plots are useful to asses model convergence and general diagnosticl information. 

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
plot(data$T, data$trait, xlim = c(5, 40), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

# Fit a gamma distribution to each parameter
gamma.fits.MDR = apply(MDR.samps, 2, function(df) fitdistr(df, "gamma")$estimate)


#####################################################
## The next parameter is pEA, the probability a mosquito will survive from hatching to 
## maturation. 

data <- data.all[which(data.all$trait.name=="e2a"),]

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


# Fit a gamma distribution to each parameter
e2a.samps.pos = e2a.samps[,c(1:4,6)]
e2a.samps.pos$qd <- -e2a.samps.pos$qd
gamma.fits.e2a = apply(e2a.samps.pos, 2, function(df) fitdistr(df, "gamma")$estimate)


######################################################
### Lifespan, lf
## Choose the response variable
data <- data.all[which(data.all$trait.name=="1/mu"),]


plot(trait ~ T, data = data)

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
lf.samps <- samps.q

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
plot.hists(e2a.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="quad", lf.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim=c(15, 35), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab="Mu (Survival Rate) - Quad w/o Informed Priors",
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

# Fit a gamma distribution to each parameter
lf.samps.pos = lf.samps[,c(1:4,6)]
lf.samps.pos$qd <- -lf.samps.pos$qd
gamma.fits.lf = apply(lf.samps.pos, 2, function(df) fitdistr(df, "gamma")$estimate)


#######################################################
### Fit bc and PDR from Lambrechts et al. 2011
data.all <- read.csv("aegyptiDENVmodelTempData_2016-03-30.csv", header=TRUE)
lambc = subset(data.all, ref=="Lambrects_et_al_2011_PNAS")

## Fit b
data = subset(lambc, trait.name=="b")

plot(trait ~ T, data = data)

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

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", b.samps, Temps, thinned=seq(1,n.samps, length=1000))

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

# Use moment matching to fit a gamma distribution to each parameter estimate
gamma.fits.b = apply(b.samps, 2, function(df) fitdistr(df, "gamma")$estimate)

#####################################
## Fit c
data = subset(lambc, trait.name=="c")

plot(trait ~ T, data = data)

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

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", c.samps, Temps, thinned=seq(1,n.samps, length=1000))


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


# Fit a gamma distribution to each parameter
gamma.fits.c = apply(c.samps, 2, function(df) fitdistr(df, "gamma")$estimate)


##############################################################
### PDR
### Data from Reisen et al from viruses in Culex mosquitoes: WNV, WEEV, SLEV
data = read.csv("EIP_priors_2015-12-04.csv", header=T)
# data = subset(data.all, trait.name=="EIP")
# data$trait <- 1/data$trait

plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data, virus=="WEEV"), pch=16, col=2)
points(trait ~ T, data = subset(data, virus=="SLEV"), pch=16, col=3)
points(trait ~ T, data = subset(data, virus=="WNV-SA"), pch=16, col=4)

# Given the data the Briere fuction is chosen. Jag-briere.bug contains the specifics of 
# the Briere model with the default priors.

jags <- jags.model('jags-briere-PDR-prior.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N' = length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
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
plot(data$T, data$trait, xlim = c(10, 40),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

# Fit a gamma distribution to each parameter
gamma.fits.PDR = apply(PDR.samps, 2, function(df) fitdistr(df, "gamma")$estimate)


##############################################################
### Save everything
save(a.samps, EFD.samps, TFD.samps, e2a.samps, MDR.samps, lf.samps, b.samps, c.samps, PDR.samps, 
      file = "aedes_prior_samps.Rsave")
save(gamma.fits.a, gamma.fits.EFD, gamma.fits.TFD, gamma.fits.e2a, gamma.fits.MDR, gamma.fits.lf, 
     gamma.fits.b, gamma.fits.c, gamma.fits.PDR, 
     file = "aedes_prior_gamma_fits.Rsave")