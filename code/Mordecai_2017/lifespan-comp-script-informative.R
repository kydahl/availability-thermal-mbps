## Load relevant packages and materials
library(IDPmisc)
library('rjags')
library(tidyverse)

# This file contains tools for analysis and visualization.
source("code/Mordecai_2017/mcmc_utils_all.R")

# This file contains the thermal response functions and their derivatives.
source("code/Mordecai_2017/temp_functions_all.R")
## Load data for first model fit
data.all <- read.csv("data/raw/Mordecai_2017/aegyptiDENVmodelTempData_2016-03-30.csv", header=TRUE)

load("code/Mordecai_2017/aedes_prior_gamma_fits.Rsave")

## Specify the MCMC Input Parameters
n.chains <- 5
n.adapt<-5000
n.samps<-5000

## Choose the response variable
## The Rohani data have unrealistically long lifespans; omit them
# data.p <- data.all[which(data.all$trait.name=="p"),]
data.days <- data.all[which(data.all$trait.name=="p/days"),]
data.1mu <- data.all[which(data.all$trait.name=="1/mu"),]

# Convert Yang data to lifespan
data.days$trait <- 1/data.days$trait


## Manipulating the survival data so we can examine it as mu 
ec <- 0.000001


## Rohani et al. 2009 data were the proportion surviving for 31, 27, 23 days at 26, 28, and 30C, respectively
## Convert to mortality
# data.p$trait = -log(data.p$trait)/rep(c(31,27,23),2)
# data.p$trait <- 1/data.p$trait
data <- rbind(data.days, data.1mu) #data.p,   

plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data, trait.name=="p/days"), col=2, pch=16)
points(trait ~ T, data = subset(data, trait.name=="1/mu"), col=3, pch=16)

hypers = gamma.fits.lf*0.01

## Set up the jags code, which contains the specifics of the quadratic model with the
## default priors. As well as a few different jags models that contain different things.

jags <- jags.model('code/Mordecai2017/jags-quad-neg-informative.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T), 'hypers' = hypers),
                   n.chains = n.chains,
                   inits=list(T0=5, Tm=33, n.qd=0.005), n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis

coda.samps <- coda.samples(jags, c('T0','Tm', 'qd', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnostic information. 
plot(coda.samps, ask=T)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps.q <- make.quad.samps(coda.samps, nchains=n.chains, 
                           samp.lims=c(1, n.samps), sig=TRUE)
samps.q$n.qd <- samps.q$qd
samps.q$tau = 1/samps.q$sigma
lf.DENV.samps <- samps.q


## Next we want to visualize the posterior samples of the parameters
## specifying the temp response, and compare them to the priors. We
## can look at the pair-wise joint posterior distirbution:
# ipairs(samps1[,c(1:3,4)], ztransf = function(x){x[x<1] <- 1; log2(x)})

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
priors1<-list()
priors1$names<-c( "T0", "Tm", "qd","tau")
priors1$fun<-c( "gamma", "gamma", "gamma","gamma")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(hypers[1,1], hypers[2,1], NA)
priors1$hyper[,2]<-c(hypers[1,2], hypers[2,2], NA)
priors1$hyper[,3]<-c(hypers[1,3], hypers[2,3], NA)
priors1$hyper[,4]<-c(hypers[1,4], hypers[2,4], NA)


## the plot.hists command can be used to plot histograms of posterior
## samples of parameters with overlying priors
plot.hists(cbind(lf.DENV.samps[,1:2], -lf.DENV.samps[,3], lf.DENV.samps[,4]), my.par=c(2,2), n.hists=4, priors=priors1)

## Next we want to use the parameter samples to get posterior samples
## of the temperature responses themselves
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="quad", lf.DENV.samps, Temps, thinned=seq(1,25000, by=5), trunc.num=0.0001)

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim=c(0,45), ylim = c(0,60),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab="Lifespan (days)",
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)
# Add a line plotting the mean function specified by the priors
lines(Temps, -hypers[1,3]/hypers[1,3]*(Temps- hypers[1,1]/hypers[2,1])*(Temps-hypers[1,2]/hypers[2,2]))


########################################
## Second Set of Data which will be used
data.all <- read.csv("code/Mordecai2017/albopictusCHIKVmodelTempData_2016-03-26.csv", header=TRUE)
data1 <- data.all[which(data.all$trait.name=="prop.dead"),]

## Manipulating the data so it can be used in the analysis
ec <- 0.000001
data1$trait = -log(1-data1$trait)

data2 = data.all[which(data.all$trait.name=="1/mu"),]
data2$trait <- 1/data2$trait

data = rbind(data1, data2)
data$trait <- 1/data$trait

plot(trait ~ T, data = data)
points(trait ~ T, data = subset(data, ref=="Alto_et_al_2001b"), col=2, pch=16)
points(trait ~ T, data = subset(data, ref=="Ezeakacha_Dissertation_2015" & trait2=="starved"), col = 3, pch = 16)
points(trait ~ T, data = subset(data, ref=="Ezeakacha_Dissertation_2015" & trait2=="sugar-fed"), col = 4, pch = 16)

# The starved mosquitoes had much shorter survival than all other data, so remove them
data = subset(data, trait2 %in% c("sugar-fed", NA))

## Set up the jags code, which contains the specifics of the quadratic model with the
## default priors. As well as a few different jags models that contain different things.

jags <- jags.model('code/Mordecai2017/jags-quad-neg-informative.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T), 'hypers' = hypers),
                   n.chains = n.chains,
                   inits=list(T0=5, Tm=33, n.qd=0.005), n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis

coda.samps <- coda.samples(jags, c('T0','Tm', 'qd', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnostic information. 
plot(coda.samps, ask=T)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps.q <- make.quad.samps(coda.samps, nchains=n.chains, 
                           samp.lims=c(1, n.samps), sig=TRUE)
samps.q$n.qd <- samps.q$qd
samps.q$tau = 1/samps.q$sigma
lf.CHIKV.samps <- samps.q

## Next we want to visualize the posterior samples of the parameters
## specifying the temp response, and compare them to the priors. We
## can look at the pair-wise joint posterior distirbution:
# ipairs(samps1[,c(1:3,4)], ztransf = function(x){x[x<1] <- 1; log2(x)})

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
priors1<-list()
priors1$names<-c( "T0", "Tm", "qd","tau")
priors1$fun<-c( "gamma", "gamma", "gamma","gamma")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(hypers[1,1], hypers[2,1], NA)
priors1$hyper[,2]<-c(hypers[1,2], hypers[2,2], NA)
priors1$hyper[,3]<-c(hypers[1,3], hypers[2,3], NA)
priors1$hyper[,4]<-c(hypers[1,4], hypers[2,4], NA)


## the plot.hists command can be used to plot histograms of posterior
## samples of parameters with overlying priors
plot.hists(cbind(lf.CHIKV.samps[,1:2], -lf.CHIKV.samps[,3], lf.CHIKV.samps[,4]), my.par=c(2,2), n.hists=4, priors=priors1)

## Next we want to use the parameter samples to get posterior samples
## of the temperature responses themselves
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="quad", lf.CHIKV.samps, Temps, thinned=seq(1,25000, by=5), trunc.num=0.0001)

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim=c(0,45), ylim=c(0,160),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab="Lifespan (days)",
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)
# Add a line plotting the mean function specified by the priors
lines(Temps, -hypers[1,3]/hypers[1,3]*(Temps- hypers[1,1]/hypers[2,1])*(Temps-hypers[1,2]/hypers[2,2]))


## Saving the mu samps 
save(lf.DENV.samps, lf.CHIKV.samps,
     file = "LifespanFits_2016-03-30-informative.Rsave")
