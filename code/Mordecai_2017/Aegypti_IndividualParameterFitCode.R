### This file will process the temperature data and create the individual parameter 
### samples for all the traits that are included in the DENV R0 model. 

setwd("C:/Users/kd99491/Documents/GitHub/thermal-properties-mbps/code/Mordecai2017")
## Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')

# This file contains tools for analysis and visualization.
source("mcmc_utils_all.R") 

# This file contains the thermal response functions
source("temp_functions_all.R") 

# Creating a small constant to keep denominators from being zero.
ec<-0.000001 

## Loading the data; my data will be called aegyptiDENVmodelTempData.csv 
data.all <- read.csv("aegyptiDENVmodelTempData_2016-03-30.csv", header=TRUE)

# Exclude the Focks & Barrera 2006 data because they're from a model
data.all = subset(data.all, ref!=as.character("Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper"))

## Now the code will choose all the temperature sensitive traits and fit either a
## Briere or Quadratic model to it using MCMC sampling.
# Specifing the parameters that control the MCMC (these will be used throughout the code). 

n.chains <- 5
n.adapt <- 5000
n.samps <- 5000

######################################
## The first trait the biting rate, a.

data <- data.all[ which(data.all$trait.name=="GCR"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

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
plot(data$T, data$trait, xlim = c(15, 40), ylim = c(0, 0.42),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

#############################################
## The second trait to analyze is EFD, or the number of eggs laid
## per female per day.

data <- data.all[which(data.all$trait.name=="EFD"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

nlfit = nls(trait ~ briere(T, c, Tm, T0), data=data,
            start=list(c=0.1, Tm=35, T0=10), algorithm="port",
            lower=c(0, max(data$T), 0))
lines(Temps, briere(Temps, coef(nlfit)[1], coef(nlfit)[2], coef(nlfit)[3]))

# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere-EFD.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 34, T0 = 15, c = 0.007),
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

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
plot.hists(EFD.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", EFD.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(15, 40), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

###########################################
## The next trait we'll be choosing is the vector competence, b*c.
## Due to the data that was collected it will be necessary
## to decompose it into its two parts, b and c.

# Choose b, the probability a human will be bitten and infected by
# an infectious mosquito (ie. transmission).

data <- data.all[ which(data.all$trait.name=="b"),]
### Can remove Lambrechts data since they are from other mosquitoes/flaviviruses
data = subset(data, ref!="Lambrects_et_al_2011_PNAS")

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data, xlim = c(15, 40), ylim = c(0, 1))
points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=2, pch=16)
points(trait ~ T, data = subset(data, ref=="Alto&Bettinardi_2013_AJTMH"), col=3, pch=16)


# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere-trunc-b.bug',
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
out1<-make.sims.temp.resp(sim="briere.trunc", b.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(5, 40), ylim = c(0,1),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)
points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=2, pch=16)
points(trait ~ T, data = subset(data, ref=="Alto&Bettinardi_2013_AJTMH"), col=3, pch=16)

#############################################
# Next, choose c, the probability that a mosquito becomes infected
# after biting an infectious human (ie. infection). 

data <- data.all[ which(data.all$trait.name=="c"),]

### Can remove Lambrechts data since they are from other mosquitoes/flaviviruses
data = subset(data, ref!="Lambrects_et_al_2011_PNAS")
data = subset(data, ref!="Alto&Bettinardi_2013_AJTMH")

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=2, pch=16)
points(trait ~ T, data = subset(data, ref=="Alto&Bettinardi_2013_AJTMH"), col=3, pch=16)


# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere-trunc.bug',
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
out1<-make.sims.temp.resp(sim="briere.trunc", c.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(5, 40), ylim = c(0,1),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

### Older version of the code, DENV_IndividualParameterFitCode.R
### fits a quadratic function to c to compare. Omitted here.

########################################
## The next trait that is fitted is MDR, the mean development time for a mosquito. 

data <- data.all[ which(data.all$trait.name=="MDR"),]

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



#####################################################
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

coda.samps <- coda.samples(jags, c('T0','Tm', 'qd'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 
plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps.q <- make.quad.samps(coda.samps, nchains=n.chains, 
                           samp.lims=c(1, n.samps), sig=FALSE)
samps.q$n.qd <- samps.q$qd
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
## The next parameter is p, the daily probability that an adult mosquito
## survives. Which will be analyzed as lf in lifespan-comp-script.R and the
## samples from that analysis will be saved and used in the overall R0
## analysis and computation. 

#########################################
## The final parameter to fit is PDR, the parasite/virus development rate.

data1 <- data.all[which(data.all$trait.name=="PDR"),]
data2 <- data.all[which(data.all$trait.name=="EIP"),]
data2$trait <- 1/data2$trait
data <- rbind(data1, data2)

# Plot the data to see which function, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data, ref=="Davis_1932_AmJEpidemiology"), col=2, pch=16)
points(trait ~ T, data = subset(data, ref=="Focks_et_al_1995_AJTMH"), col=3, pch=16)
points(trait ~ T, data = subset(data, ref=="McLean_et_al_1974_CanJMicobiol"), col=4, pch=16)
points(trait ~ T, data = subset(data, ref=="McLeah_et_al_1975_MosquitoNews"), col=5, pch=16)
points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=6, pch=16)
legend('topleft', legend=c("Davis", "Focks", "McLean 74", "McLean 75", "Watts", "Carrington"), col=c(2:6, 1), pch=c(rep(16, 5), 1))


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
plot(data$T, data$trait, xlim = c(5, 45), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)
points(trait ~ T, data = subset(data, ref=="Davis_1932_AmJEpidemiology"), col=2, pch=16)
points(trait ~ T, data = subset(data, ref=="Focks_et_al_1995_AJTMH"), col=3, pch=16)
points(trait ~ T, data = subset(data, ref=="McLean_et_al_1974_CanJMicobiol"), col=4, pch=16)
points(trait ~ T, data = subset(data, ref=="McLeah_et_al_1975_MosquitoNews"), col=5, pch=16)
points(trait ~ T, data = subset(data, ref=="Watts_et_al_1987_AJTMH"), col=6, pch=16)
legend('topleft', legend=c("Davis", "Focks", "McLean 74", "McLean 75", "Watts", "Carrington"), col=c(2:6, 1), pch=c(rep(16, 5), 1))


## This code is just save the MCMC samples for further analysis in the R0 model.
save(a.samps, b.samps, c.samps, MDR.samps, EFD.samps, 
     e2a.samps, PDR.samps,
     file = "Aegypti_DENV_ParameterFits_2016-03-30.Rsave")
save(PDR.samps, b.samps, c.samps, file = "Aegypti_DENV_b_c_PDR_samps_2016-03-30.Rsave")
