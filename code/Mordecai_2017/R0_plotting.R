# Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')
library(RColorBrewer)

# This file contains tools for analysis and visualization.
source("mcmc_utils_all.R") 

# This file contains the thermal response functions and their derivatives
source("temp_functions_all.R") 

# Temperature sequence
temp = seq(5,45,by=0.1)  
t<-length(temp)


## Loading in samples from the Rsave (CHIKV_ParameterFits.Rsave), which have the samples
## for a, b, c, MDR, TFD, e2a created in the CHIKV_IndividualParameterFitCode.R
# load("Albopictus_ParameterFits_2016-03-30.Rsave")
load("Informative_Albopictus_ParameterFits_2016-03-30.Rsave")

# load("LifespanFits_2016-03-30.Rsave")
load("LifespanFits_2016-03-30-informative.Rsave")
lf.samps <- lf.CHIKV.samps

# Length of samples, assuming they're all the same length
n = dim(a.samps)[1]
# n.mu = dim(mu.samps)[1]

# Thinned samples
thinned<-seq(1, n, by = 5) 
lthin = length(thinned)

# Creating a small constant to keep denominators from being zero.
ec<-0.000001           

### Albopictus R0 model plotting and validation
# load("Albopictus_model_outputs-uninformative.Rsave")
load("Albopictus_model_outputs-informative.Rsave")

albo.R0 = R0
albo.R0.M = R0.M
a.M.albo = a.M
b.M.albo = b.M
c.M.albo = c.M
MDR.M.albo = MDR.M
e2a.M.albo = e2a.M
PDR.M.albo = PDR.M
TFD.M.albo = TFD.M
lf.M.albo = lf.M
a.albo = a
b.albo = b
c.albo = c
MDR.albo = MDR
e2a.albo = e2a
PDR.albo = PDR
TFD.albo = TFD
lf.albo = lf
a.samps.albo = a.samps
TFD.samps.albo = TFD.samps
e2a.samps.albo = e2a.samps
MDR.samps.albo = MDR.samps
lf.samps.albo = lf.samps
b.samps.albo = b.samps
c.samps.albo = c.samps
PDR.samps.albo = PDR.samps

load("albo_R0DTR8-informative.Rsave")
# loads hourlytemps and R0DTR8
albo.R0DTR8 = R0DTR8
albo.hourlytemps = hourlytemps
albo.R0DTR8.M = colMeans(albo.R0DTR8)
albo.temps.DTR8 = colMeans(albo.hourlytemps)

# Calculate the distribution of the lower and upper limits of R0 
# and peak R0.

albo.R0.min<-albo.R0.max<-albo.R0.peak<-rep(NA, length(thinned))

# Plotting the PEAK albo.R0 distribution.
for(i in 1:length(thinned)){
  ww<-which(albo.R0[,i]==max(albo.R0[,i]))
  albo.R0.peak[i]<-temp[ww[1]]
}

# Plotting the MINIMUM albo.R0 distribution.
for(i in 1:length(thinned)){
  ww<-which(albo.R0[,i]>0)
  albo.R0.min[i]<-temp[ww[1]-1]
}

# Plotting the MAXIMUM albo.R0 distribution.

for(i in 1:length(thinned)){
  ww<-which(albo.R0[,i]>0)
  lw<-length(ww)
  albo.R0.max[i]<-temp[ww[lw]+1]
}

# Plotting the mean albo.R0 with it's quantiles, all scaled by max mean albo.R0.
pdf("Albopictus_R0_histograms-informative.pdf")
# pdf("Albopictus_R0_histograms.pdf")
layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow=T))
albo.R0.scale<-max(albo.R0.M)
albo.R0.q2<-temp.sim.quants(albo.R0, length(temp))##/albo.R0.scale
plot(temp, albo.R0.M/albo.R0.scale, type="l", col=1, lwd=3, xlim=c(10, 40), ylim=c(0, 1.7),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
add.sim.lines(temp, sim.data=NULL, q=albo.R0.q2/albo.R0.scale, mycol=2)

hist(albo.R0.min, xlab=expression(paste("Temp. of min. ", R[0], sep="")), freq=TRUE, main="")
hist(albo.R0.peak, xlab=expression(paste("Temp. of peak ", R[0], sep="")), freq=TRUE, main="")
hist(albo.R0.max, xlab=expression(paste("Temp. of max. ", R[0], sep="")), freq=TRUE, main="")
par(mfrow=c(1,1), bty="n")
dev.off()

albo.results = t(matrix(c(
  # median(albo.R0.peak)
  mean(albo.R0.peak),
  HPDinterval(mcmc(albo.R0.peak)),
  # median(albo.R0.min)
  mean(albo.R0.min),
  HPDinterval(mcmc(albo.R0.min)),
  # median(albo.R0.max)
  mean(albo.R0.max),
  HPDinterval(mcmc(albo.R0.max))
), 3, 3))
rownames(albo.results) = c("peak", "min", "max")
colnames(albo.results) = c("mean", "lower 95%", "upper 95%")
albo.results

# Calculate the distribution of the lower and upper limits of R0 
# and peak R0 for 8DTR

albo.R0DTR8.min<-albo.R0DTR8.max<-albo.R0DTR8.peak<-rep(NA, length(thinned))

# Plotting the PEAK albo.R0DTR8 distribution.
for(i in 1:length(thinned)){
  ww<-which(albo.R0DTR8[i,]==max(albo.R0DTR8[i,]))
  albo.R0DTR8.peak[i]<-albo.temps.DTR8[ww[1]]
}

# Plotting the MINIMUM albo.R0DTR8 distribution.
for(i in 1:length(thinned)){
  ww<-which(albo.R0DTR8[i,]>0)
  albo.R0DTR8.min[i]<-albo.temps.DTR8[ww[1]-1]
}

# Plotting the MAXIMUM albo.R0DTR8 distribution.

for(i in 1:length(thinned)){
  ww<-which(albo.R0DTR8[i,]>0)
  lw<-length(ww)
  albo.R0DTR8.max[i]<-albo.temps.DTR8[ww[lw]+1]
}

albo.DTR8.results = t(matrix(c(
  # median(albo.R0DTR8.peak)
  mean(albo.R0DTR8.peak),
  HPDinterval(mcmc(albo.R0DTR8.peak)),
  # median(albo.R0DTR8.min)
  mean(albo.R0DTR8.min),
  HPDinterval(mcmc(albo.R0DTR8.min)),
  # median(albo.R0DTR8.max)
  mean(albo.R0DTR8.max),
  HPDinterval(mcmc(albo.R0DTR8.max))
), 3, 3))
rownames(albo.DTR8.results) = c("peak", "min", "max")
colnames(albo.DTR8.results) = c("mean", "lower 95%", "upper 95%")
albo.DTR8.results


## Plotting the data against the thermal responses fits.

# Subsetting the data.
data.albo <- read.csv("albopictusCHIKVmodelTempData_2016-03-26.csv", header=TRUE)

a.data.albo = subset(data.albo, trait.name=="GCD")
a.data.albo$trait <- 1/a.data.albo$trait
TFD.data.albo = subset(data.albo, trait.name=="TFD" & trait2.name %in% c("R1", "R2", "R3", NA) & ref != "Yee_et al_2016 JAE_in_review")
e2a.data.albo = subset(data.albo, trait.name=="pEA")
data1 <- data.albo[which(data.albo$trait.name=="MDR"),]
data2 <- data.albo[which(data.albo$trait.name=="1/MDR"),]
data2$trait <- 1/data2$trait
MDR.data.albo = rbind(data1, data2)
data1 <- data.albo[which(data.albo$trait.name=="prop.dead"),]

## Manipulating the data so it can be used in the analysis
data1$trait = -log(1-data1$trait)

data2 = data.albo[which(data.albo$trait.name=="1/mu"),]
data2$trait <- 1/data2$trait

lf.data.albo = rbind(data1, data2)
lf.data.albo$trait <- 1/lf.data.albo$trait

b.data.albo <- data.albo[ which(data.albo$trait.name=="b" & data.albo$ref=="Xiao_et_al_2014_Arch Virol"),]
c.data.albo <- data.albo[ which(data.albo$trait.name=="c" & data.albo$ref=="Xiao_et_al_2014_Arch Virol"),]
PDR.data.albo <- data.albo[which(data.albo$trait.name=="EIP"),]
PDR.data.albo$trait <- 1/PDR.data.albo$trait


# # Getting the HPD intervals.
# 
# a.int.albo = HPDinterval(mcmc(t(a.albo)))
# b.int.albo = HPDinterval(mcmc(t(b.albo)))
# c.int.albo = HPDinterval(mcmc(t(c.albo)))
# TFD.int.albo = HPDinterval(mcmc(t(TFD.albo)))
# e2a.int.albo = HPDinterval(mcmc(t(e2a.albo)))
# MDR.int.albo = HPDinterval(mcmc(t(MDR.albo)))
# lf.int.albo = HPDinterval(mcmc(t(lf.albo)))

# Use the built-in quantile function because the HPD intervals look funky
quant.fun = function(x, probs = c(0.025, 0.975)){
  y = matrix(NA, nrow(x), 2)
  for (i in 1:nrow(x)){
    y[i,] = quantile(x[i,], probs, na.rm=T)
  }
  y
}
a.int.albo = quant.fun(a.albo)
# b.int.albo = quant.fun(b.albo)
# c.int.albo = quant.fun(c.albo)
TFD.int.albo = quant.fun(TFD.albo)
e2a.int.albo = quant.fun(e2a.albo)
MDR.int.albo = quant.fun(MDR.albo)
# PDR.int.albo = quant.fun(PDR.albo)
lf.int.albo = quant.fun(lf.albo)
b.int.albo = quant.fun(b.albo)
c.int.albo = quant.fun(c.albo)
PDR.int.albo = quant.fun(PDR.albo)

# Look at the Tpk for each trait
pk.fun = function(df) temp[which(df==max(df, na.rm=T))]
pk.fun(a.M.albo)
pk.fun(TFD.M.albo)
pk.fun(e2a.M.albo)
pk.fun(MDR.M.albo)
pk.fun(lf.M.albo)
pk.fun(b.M.albo)
pk.fun(c.M.albo)
pk.fun(PDR.M.albo)

# Setting up the plot area.
# pdf("Albopictus_trait_responses-informative.pdf")

pdf("Albopictus_trait_responses.pdf", width = 8, height = 8)
par(mfrow=c(3,3))
plot(trait~T, data=a.data.albo, main="Biting Rate", xlim=c(10,40), ylim=c(0,0.4), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)")
lines(temp, a.M.albo)
lines(temp, a.int.albo[,1], lty=2, col=2)
lines(temp, a.int.albo[,2], lty=2, col=2)

plot(trait~T, data=TFD.data.albo, main="Fecundity", xlim=c(10,40), ylim = c(0, 83), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Eggs per female per cycle")
lines(temp, TFD.M.albo)
lines(temp, TFD.int.albo[,1], lty=2, col=2)
lines(temp, TFD.int.albo[,2], lty=2, col=2)

plot(trait~T, data=e2a.data.albo, main="Egg-to-Adult Survival", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability")
lines(temp, e2a.M.albo)
lines(temp, e2a.int.albo[,1], lty=2, col=2)
lines(temp, e2a.int.albo[,2], lty=2, col=2)

plot(trait~T, data=MDR.data.albo, main="Mosquito Development Rate", xlim=c(10,40), ylim=c(0,0.18), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)")
lines(temp, MDR.M.albo)
lines(temp, MDR.int.albo[,1], lty=2, col=2)
lines(temp, MDR.int.albo[,2], lty=2, col=2)

plot(trait~T, data=lf.data.albo, main="Adult Lifespan", xlim=c(10,40), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Lifespan (days)")
lines(temp, lf.M.albo)
lines(temp, lf.int.albo[,1], lty=2, col=2)
lines(temp, lf.int.albo[,2], lty=2, col=2)

plot(trait~T, data=b.data.albo, main="Transmission Probability", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability")
lines(temp, b.M.albo)
lines(temp, b.int.albo[,1], lty=2, col=2)
lines(temp, b.int.albo[,2], lty=2, col=2)

plot(trait~T, data=c.data.albo, main="Infection Probability", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability")
lines(temp, c.M.albo)
lines(temp, c.int.albo[,1], lty=2, col=2)
lines(temp, c.int.albo[,2], lty=2, col=2)

plot(trait~T, data=PDR.data.albo, main="Extrinsic Incubation Rate", xlim=c(10,40), ylim=c(0,0.35), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)")
lines(temp, PDR.M.albo)
lines(temp, PDR.int.albo[,1], lty=2, col=2)
lines(temp, PDR.int.albo[,2], lty=2, col=2)


par(mfrow=c(1,1))

dev.off()

### Create a version of trait plots with means and error bars
require(plyr)
require(plotrix)
trait.fun = function(df){
  mean = mean(df$trait, na.rm=T)
  se = std.error(df$trait, na.rm=T)
  data.frame(mean, se)
}


pdf("Albopictus_trait_responses_errorbars.pdf", width = 8, height = 8)
par(mfrow=c(3,3))
a.data.bars = ddply(a.data.albo, .(T), trait.fun)
plot(mean~T, data=a.data.bars, main="Biting Rate", xlim=c(10,40), ylim=c(0,0.35), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)", pch = 16)
for (i in 1:nrow(a.data.bars)) {lines(rep(a.data.bars$T[i],2), c(a.data.bars$mean[i] - a.data.bars$se[i], a.data.bars$mean[i] + a.data.bars$se[i]))}
lines(temp, a.M.albo)
lines(temp, a.int.albo[,1], lty=2, col=2)
lines(temp, a.int.albo[,2], lty=2, col=2)

TFD.data.bars = ddply(TFD.data.albo, .(T), trait.fun)
plot(mean~T, data=TFD.data.bars, main="Fecundity", xlim=c(10,40), ylim = c(0, 90), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Eggs per female per cycle", pch = 16)
for (i in 1:nrow(TFD.data.bars)) {lines(rep(TFD.data.bars$T[i],2), c(TFD.data.bars$mean[i] - TFD.data.bars$se[i], TFD.data.bars$mean[i] + TFD.data.bars$se[i]))}
lines(temp, TFD.M.albo)
lines(temp, TFD.int.albo[,1], lty=2, col=2)
lines(temp, TFD.int.albo[,2], lty=2, col=2)

e2a.data.bars = ddply(e2a.data.albo, .(T), trait.fun)
plot(mean~T, data=e2a.data.bars, main="Egg-to-Adult Survival", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability", pch = 16)
for (i in 1:nrow(e2a.data.bars)) {lines(rep(e2a.data.bars$T[i],2), c(e2a.data.bars$mean[i] - e2a.data.bars$se[i], e2a.data.bars$mean[i] + e2a.data.bars$se[i]))}
lines(temp, e2a.M.albo)
lines(temp, e2a.int.albo[,1], lty=2, col=2)
lines(temp, e2a.int.albo[,2], lty=2, col=2)

MDR.data.bars = ddply(MDR.data.albo, .(T), trait.fun)
plot(mean~T, data=MDR.data.bars, main="Mosquito Development Rate", xlim=c(10,40), ylim=c(0,0.2), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)", pch = 16)
for (i in 1:nrow(MDR.data.bars)) {lines(rep(MDR.data.bars$T[i],2), c(MDR.data.bars$mean[i] - MDR.data.bars$se[i], MDR.data.bars$mean[i] + MDR.data.bars$se[i]))}
lines(temp, MDR.M.albo)
lines(temp, MDR.int.albo[,1], lty=2, col=2)
lines(temp, MDR.int.albo[,2], lty=2, col=2)

lf.data.bars = ddply(lf.data.albo, .(T), trait.fun)
plot(mean~T, data=lf.data.bars, main="Adult Lifespan", xlim=c(10,40), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Lifespan (days)", pch = 16, ylim = c(15, 145))
for (i in 1:nrow(lf.data.bars)) {lines(rep(lf.data.bars$T[i],2), c(lf.data.bars$mean[i] - lf.data.bars$se[i], lf.data.bars$mean[i] + lf.data.bars$se[i]))}
lines(temp, lf.M.albo)
lines(temp, lf.int.albo[,1], lty=2, col=2)
lines(temp, lf.int.albo[,2], lty=2, col=2)

b.data.bars = ddply(b.data.albo, .(T), trait.fun)
plot(mean~T, data=b.data.bars, main="Transmission Probability", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability", pch = 16)
for (i in 1:nrow(b.data.bars)) {lines(rep(b.data.bars$T[i],2), c(b.data.bars$mean[i] - b.data.bars$se[i], b.data.bars$mean[i] + b.data.bars$se[i]))}
lines(temp, b.M.albo)
lines(temp, b.int.albo[,1], lty=2, col=2)
lines(temp, b.int.albo[,2], lty=2, col=2)

c.data.bars = ddply(c.data.albo, .(T), trait.fun)
plot(mean~T, data=c.data.bars, main="Infection Probability", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability", pch = 16)
for (i in 1:nrow(c.data.bars)) {lines(rep(c.data.bars$T[i],2), c(c.data.bars$mean[i] - c.data.bars$se[i], c.data.bars$mean[i] + c.data.bars$se[i]))}
lines(temp, c.M.albo)
lines(temp, c.int.albo[,1], lty=2, col=2)
lines(temp, c.int.albo[,2], lty=2, col=2)

PDR.data.bars = ddply(PDR.data.albo, .(T), trait.fun)
plot(mean~T, data=PDR.data.bars, main="Extrinsic Incubation Rate", xlim=c(10,40), ylim=c(0,0.35), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)", pch = 16)
for (i in 1:nrow(PDR.data.bars)) {lines(rep(PDR.data.bars$T[i],2), c(PDR.data.bars$mean[i] - PDR.data.bars$se[i], PDR.data.bars$mean[i] + PDR.data.bars$se[i]))}
lines(temp, PDR.M.albo)
lines(temp, PDR.int.albo[,1], lty=2, col=2)
lines(temp, PDR.int.albo[,2], lty=2, col=2)

par(mfrow=c(1,1))

dev.off()




### Create data source and fit table
albo.traits = c("a", "TFD", "e2a", "MDR", "lf", "b", "c", "PDR")
albo.refs = c(paste(unique(a.data.albo$ref), collapse=", "), paste(unique(TFD.data.albo$ref), collapse=", "), paste(unique(e2a.data.albo$ref), collapse=", "), paste(unique(MDR.data.albo$ref), collapse=", "), paste(unique(lf.data.albo$ref), collapse=", "), paste(unique(b.data.albo$ref), collapse=", "), paste(unique(c.data.albo$ref), collapse=", "), paste(unique(PDR.data.albo$ref), collapse=", "))
### summarize the posterior samples for each trait fit
fit.fun = function(df){
  mean = colMeans(df)[1:3]
  tmp = HPDinterval(mcmc(df))
  lower = tmp[1:3,1]
  upper = tmp[1:3,2]
  out = c(mean, lower, upper)
  names(out) = c("T0.mean", "Tm.mean", "c.mean", "T0.lower", "Tm.lower", "c.lower", "T0.upper", "Tm.upper", "c.upper")
  out
}
fits = t(sapply(list(a.samps.albo, TFD.samps.albo, e2a.samps.albo, MDR.samps.albo, lf.samps.albo, b.samps.albo, c.samps.albo, PDR.samps.albo), fit.fun))
albo.table = data.frame('trait' = albo.traits, 'refs' = albo.refs, fits)
write.csv(albo.table, file = "Albopictus_trait_data_fits.csv", row.names = F)


### Plot R0 alone
# pdf("Albopictus_R0_alone-informative.pdf", width = 8.5, height = 7)
# pdf("Albopictus_R0_alone.pdf", width = 8.5, height = 7)
# pdf("Albopictus_R0_with_DTR8-informative.pdf", width = 8.5, height = 7)
albo.R0.scale = max(albo.R0.M)
par(mfrow = c(1,1), mar=c(6, 4, 4, 5))
plot(temp, albo.R0.M/albo.R0.scale, ylim=c(0, 1.05), xlim = c(15, 35), type="l", col="dodgerblue", lwd=3, xlab="Temperature (C)", ylab=expression(paste("Relative ", R[0], " (line)", sep="")))

lines(albo.temps.DTR8, albo.R0DTR8.M/max(albo.R0DTR8.M), col = "gray40", lwd = 3)
legend('topleft', c("DTR 0 C", "DTR 8 C"), col = c("dodgerblue", "gray40"), lty = 1, bty = 'n', lwd = 3)
# add.sim.lines(temp, sim.data=NULL, q=R0.q2/R0.scale, mycol="gray")

# dev.off()

###############################################
### Aegypti DENV model
# Load previous results
# load("Aegypti_DENV_model_outputs-uninformative.Rsave")
load("Aegypti_DENV_model_outputs-informative.Rsave")

# load("Aegypti_DENV_ParameterFits_2016-03-30.Rsave")
load("Informative_Aegypti_DENV_ParameterFits_2016-03-30.Rsave")

# load("LifespanFits_2016-03-30.Rsave")
load("LifespanFits_2016-03-30-informative.Rsave")
lf.samps <- lf.DENV.samps

aegy.R0 = R0
aegy.R0.M = R0.M
a.M.aegy = a.M
b.M.aegy = b.M
c.M.aegy = c.M
MDR.M.aegy = MDR.M
e2a.M.aegy = e2a.M
PDR.M.aegy = PDR.M
EFD.M.aegy = EFD.M
lf.M.aegy = lf.M
a.aegy = a
b.aegy = b
c.aegy = c
MDR.aegy = MDR
e2a.aegy = e2a
PDR.aegy = PDR
EFD.aegy = EFD
lf.aegy = lf
a.samps.aegy = a.samps
b.samps.aegy = b.samps
c.samps.aegy = c.samps
MDR.samps.aegy = MDR.samps
e2a.samps.aegy = e2a.samps
EFD.samps.aegy = EFD.samps
lf.samps.aegy = lf.samps
PDR.samps.aegy = PDR.samps

load("aegy_R0DTR8-informative.Rsave")
# loads hourlytemps and R0DTR8
aegy.R0DTR8 = R0DTR8
aegy.hourlytemps = hourlytemps
aegy.R0DTR8.M = colMeans(aegy.R0DTR8)
aegy.temps.DTR8 = colMeans(aegy.hourlytemps)


# Calculate the distribution of the lower and upper limits of R0 
# and peak R0.
aegy.R0col = ncol(aegy.R0)
aegy.R0.min<-aegy.R0.max<-aegy.R0.peak<-rep(NA, aegy.R0col)

# Plotting the PEAK aegy.R0 distribution.
for(i in 1:aegy.R0col){
  ww<-which(aegy.R0[,i]==max(aegy.R0[,i]))
  aegy.R0.peak[i]<-temp[ww[1]]
}

# Plotting the MINIMUM aegy.R0 distribution.
for(i in 1:aegy.R0col){
  ww<-which(aegy.R0[,i]>0)
  aegy.R0.min[i]<-temp[ww[1]-1]
}

# Plotting the MAXIMUM aegy.R0 distribution.

for(i in 1:aegy.R0col){
  ww<-which(aegy.R0[,i]>0)
  lw<-length(ww)
  aegy.R0.max[i]<-temp[ww[lw]+1]
}

# Plotting the mean aegy.R0 with its quantiles, all scaled by max mean aegy.R0
pdf("Aegypti_R0_histograms-informative.pdf")
# pdf("Aegypti_R0_histograms.pdf")
layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow=T))
aegy.R0.scale<-max(aegy.R0.M)
aegy.R0.q2<-temp.sim.quants(aegy.R0, length(temp))##/aegy.R0.scale
plot(temp, aegy.R0.M/aegy.R0.scale, type="l", col=1, lwd=3, xlim=c(10, 40), ylim=c(0, 1.5),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
add.sim.lines(temp, sim.data=NULL, q=aegy.R0.q2/aegy.R0.scale, mycol=2)

hist(aegy.R0.min, xlab=expression(paste("Temp. of min. ", R[0], sep="")), freq=TRUE, main="")
hist(aegy.R0.peak, xlab=expression(paste("Temp. of peak ", R[0], sep="")), freq=TRUE, main="")
hist(aegy.R0.max, xlab=expression(paste("Temp. of max. ", R[0], sep="")), freq=TRUE, main="")
par(mfrow=c(1,1), bty="n")
dev.off()

aegy.results = t(matrix(c(
  # median(aegy.R0.peak)
  mean(aegy.R0.peak),
  HPDinterval(mcmc(aegy.R0.peak)),
  # median(aegy.R0.min)
  mean(aegy.R0.min),
  HPDinterval(mcmc(aegy.R0.min)),
  # median(aegy.R0.max)
  mean(aegy.R0.max),
  HPDinterval(mcmc(aegy.R0.max))
), 3, 3))
rownames(aegy.results) = c("peak", "min", "max")
colnames(aegy.results) = c("mean", "lower 95%", "upper 95%")
aegy.results

# Calculate the distribution of the lower and upper limits of R0 
# and peak R0 for 8DTR

aegy.R0DTR8.min<-aegy.R0DTR8.max<-aegy.R0DTR8.peak<-rep(NA, length(thinned))

# Plotting the PEAK aegy.R0DTR8 distribution.
for(i in 1:length(thinned)){
  ww<-which(aegy.R0DTR8[i,]==max(aegy.R0DTR8[i,]))
  aegy.R0DTR8.peak[i]<-aegy.temps.DTR8[ww[1]]
}

# Plotting the MINIMUM aegy.R0DTR8 distribution.
for(i in 1:length(thinned)){
  ww<-which(aegy.R0DTR8[i,]>0)
  aegy.R0DTR8.min[i]<-aegy.temps.DTR8[ww[1]-1]
}

# Plotting the MAXIMUM aegy.R0DTR8 distribution.

for(i in 1:length(thinned)){
  ww<-which(aegy.R0DTR8[i,]>0)
  lw<-length(ww)
  aegy.R0DTR8.max[i]<-aegy.temps.DTR8[ww[lw]+1]
}

aegy.DTR8.results = t(matrix(c(
  # median(aegy.R0DTR8.peak)
  mean(aegy.R0DTR8.peak),
  HPDinterval(mcmc(aegy.R0DTR8.peak)),
  # median(aegy.R0DTR8.min)
  mean(aegy.R0DTR8.min),
  HPDinterval(mcmc(aegy.R0DTR8.min)),
  # median(aegy.R0DTR8.max)
  mean(aegy.R0DTR8.max),
  HPDinterval(mcmc(aegy.R0DTR8.max))
), 3, 3))
rownames(aegy.DTR8.results) = c("peak", "min", "max")
colnames(aegy.DTR8.results) = c("mean", "lower 95%", "upper 95%")
aegy.DTR8.results

## Plotting the data against the thermal responses fits.

# Subsetting the data.
data.aegy <- read.csv("aegyptiDENVmodelTempData_2016-03-30.csv", header=TRUE)
data.aegy = subset(data.aegy, ref!=as.character("Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper"))
a.data.aegy = subset(data.aegy, trait.name=="GCR")
b.data.aegy = subset(data.aegy, trait.name=="b")
b.data.aegy = subset(b.data.aegy, ref!="Lambrects_et_al_2011_PNAS")
c.data.aegy = subset(data.aegy, trait.name=="c")
c.data.aegy = subset(c.data.aegy, ref!="Lambrects_et_al_2011_PNAS")
c.data.aegy = subset(c.data.aegy, ref!="Alto&Bettinardi_2013_AJTMH")
EFD.data.aegy = subset(data.aegy, trait.name=="EFD")
e2a.data.aegy = subset(data.aegy, trait.name=="pEA")
MDR.data.aegy = subset(data.aegy, trait.name=="MDR")
PDR.data.aegy1 <- data.aegy[which(data.aegy$trait.name=="PDR"),]
PDR.data.aegy2 <- data.aegy[which(data.aegy$trait.name=="EIP"),]
PDR.data.aegy2$trait <- 1/PDR.data.aegy2$trait
PDR.data.aegy <- rbind(PDR.data.aegy1, PDR.data.aegy2)
# data.p <- data.aegy[which(data.aegy$trait.name=="p"),]
data.days <- data.aegy[which(data.aegy$trait.name=="p/days"),]
data.days$trait = 1/data.days$trait
# data.p$trait = -log(data.p$trait)/rep(c(31,27,23),2)
# mu.data.aegy <- rbind(data.days, data.p)
data.1mu <- data.aegy[which(data.aegy$trait.name=="1/mu"),]
lf.data.aegy <- rbind(data.days, data.1mu)


# # Getting the HPD intervals.
# 
# a.int.aegy = HPDinterval(mcmc(t(a.aegy)))
# b.int.aegy = HPDinterval(mcmc(t(b.aegy)))
# c.int.aegy = HPDinterval(mcmc(t(c.aegy)))
# EFD.int.aegy = HPDinterval(mcmc(t(EFD.aegy)))
# e2a.int.aegy = HPDinterval(mcmc(t(e2a.aegy)))
# MDR.int.aegy = HPDinterval(mcmc(t(MDR.aegy)))
# PDR.int.aegy = HPDinterval(mcmc(t(PDR.aegy)))
# mu.int.aegy = HPDinterval(mcmc(t(mu.aegy)))

# Use the built-in quantile function because the HPD intervals look funky
quant.fun = function(x, probs = c(0.025, 0.975)){
  y = matrix(NA, nrow(x), 2)
  for (i in 1:nrow(x)){
    y[i,] = quantile(x[i,], probs, na.rm=T)
  }
  y
}
a.int.aegy = quant.fun(a.aegy)
b.int.aegy = quant.fun(b.aegy)
c.int.aegy = quant.fun(c.aegy)
EFD.int.aegy = quant.fun(EFD.aegy)
e2a.int.aegy = quant.fun(e2a.aegy)
MDR.int.aegy = quant.fun(MDR.aegy)
PDR.int.aegy = quant.fun(PDR.aegy)
lf.int.aegy = quant.fun(lf.aegy)

# Look at the Tpk for each trait
pk.fun = function(df) temp[which(df==max(df, na.rm=T))]
pk.fun(a.M.aegy)
pk.fun(b.M.aegy)
pk.fun(c.M.aegy)
pk.fun(EFD.M.aegy)
pk.fun(e2a.M.aegy)
pk.fun(MDR.M.aegy)
pk.fun(PDR.M.aegy)
pk.fun(lf.M.aegy)

# Setting up the plot area.
# pdf("Aegypti_DENV_trait_responses-informative.pdf")

pdf("Aegypti_DENV_trait_responses.pdf", width=8, height=8)
par(mfrow=c(3,3))
plot(trait~T, data=a.data.aegy, main="Biting Rate", xlim=c(10,40), ylim=c(0,0.4), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)")
lines(temp, a.M.aegy)
lines(temp, a.int.aegy[,1], lty=2, col=2)
lines(temp, a.int.aegy[,2], lty=2, col=2)

plot(trait~T, data=EFD.data.aegy, main="Fecundity", xlim=c(10,40), ylim=c(0,12), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Eggs per female per day")
lines(temp, EFD.M.aegy)
lines(temp, EFD.int.aegy[,1], lty=2, col=2)
lines(temp, EFD.int.aegy[,2], lty=2, col=2)

plot(trait~T, data=e2a.data.aegy, main="Egg-to-Adult Survival", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability")
lines(temp, e2a.M.aegy)
lines(temp, e2a.int.aegy[,1], lty=2, col=2)
lines(temp, e2a.int.aegy[,2], lty=2, col=2)

plot(trait~T, data=MDR.data.aegy, main="Mosquito Development Rate", xlim=c(10,40), ylim=c(0,0.18), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)")
lines(temp, MDR.M.aegy)
lines(temp, MDR.int.aegy[,1], lty=2, col=2)
lines(temp, MDR.int.aegy[,2], lty=2, col=2)

# plot(trait~T, data=mu.data.aegy, main="mu", xlim=c(10,40))
# lines(temp, mu.M.aegy)
# lines(temp, mu.int.aegy[,1], lty=2, col=2)
# lines(temp, mu.int.aegy[,2], lty=2, col=2)

plot(trait~T, data=lf.data.aegy, main="Adult Lifespan", xlim=c(10,40), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Lifespan (days)")
lines(temp, lf.M.aegy)
lines(temp, lf.int.aegy[,1], lty=2, col=2)
lines(temp, lf.int.aegy[,2], lty=2, col=2)

plot(trait~T, data=b.data.aegy, main="Transmission Probability", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability")
lines(temp, b.M.aegy)
lines(temp, b.int.aegy[,1], lty=2, col=2)
lines(temp, b.int.aegy[,2], lty=2, col=2)

plot(trait~T, data=c.data.aegy, main="Infection Probability", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability")
lines(temp, c.M.aegy)
lines(temp, c.int.aegy[,1], lty=2, col=2)
lines(temp, c.int.aegy[,2], lty=2, col=2)

plot(trait~T, data=PDR.data.aegy, main="Extrinsic Incubation Rate", xlim=c(10,40), ylim=c(0,0.25), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)")
lines(temp, PDR.M.aegy)
lines(temp, PDR.int.aegy[,1], lty=2, col=2)
lines(temp, PDR.int.aegy[,2], lty=2, col=2)

par(mfrow=c(1,1))

dev.off()

### Create a version of trait plots with means and errorbars
pdf("Aegypti_DENV_trait_responses_errorbars.pdf", width=8, height=8)
par(mfrow=c(3,3))
a.data.b = ddply(a.data.aegy, .(T), trait.fun)
plot(mean~T, data=a.data.b, main="Biting Rate", xlim=c(10,40), ylim=c(0,0.4), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)", pch = 16)
for (i in 1:nrow(a.data.b)) {lines(rep(a.data.b$T[i],2), c(a.data.b$mean[i] - a.data.b$se[i], a.data.b$mean[i] + a.data.b$se[i]))}
lines(temp, a.M.aegy)
lines(temp, a.int.aegy[,1], lty=2, col=2)
lines(temp, a.int.aegy[,2], lty=2, col=2)

EFD.aegy.round = EFD.data.aegy
EFD.aegy.round$T = round(EFD.aegy.round$T)
EFD.data.b = ddply(EFD.aegy.round, .(T), trait.fun)
plot(mean~T, data=EFD.data.b, main="Fecundity", xlim=c(10,40), ylim=c(0,14), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Eggs per female per day", pch = 16)
for (i in 1:nrow(EFD.data.b)) {lines(rep(EFD.data.b$T[i],2), c(EFD.data.b$mean[i] - EFD.data.b$se[i], EFD.data.b$mean[i] + EFD.data.b$se[i]))}
lines(temp, EFD.M.aegy)
lines(temp, EFD.int.aegy[,1], lty=2, col=2)
lines(temp, EFD.int.aegy[,2], lty=2, col=2)

e2a.data.b = ddply(e2a.data.aegy, .(T), trait.fun)
plot(mean~T, data=e2a.data.b, main="Egg-to-Adult Survival", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability", pch = 16)
for (i in 1:nrow(e2a.data.b)) {lines(rep(e2a.data.b$T[i],2), c(e2a.data.b$mean[i] - e2a.data.b$se[i], e2a.data.b$mean[i] + e2a.data.b$se[i]))}
lines(temp, e2a.M.aegy)
lines(temp, e2a.int.aegy[,1], lty=2, col=2)
lines(temp, e2a.int.aegy[,2], lty=2, col=2)

MDR.data.b = ddply(MDR.data.aegy, .(T), trait.fun)
plot(mean~T, data=MDR.data.b, main="Mosquito Development Rate", xlim=c(10,40), ylim=c(0,0.18), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)", pch = 16)
for (i in 1:nrow(MDR.data.b)) {lines(rep(MDR.data.b$T[i],2), c(MDR.data.b$mean[i] - MDR.data.b$se[i], MDR.data.b$mean[i] + MDR.data.b$se[i]))}
lines(temp, MDR.M.aegy)
lines(temp, MDR.int.aegy[,1], lty=2, col=2)
lines(temp, MDR.int.aegy[,2], lty=2, col=2)

# plot(trait~T, data=mu.data.aegy, main="mu", xlim=c(10,40))
# lines(temp, mu.M.aegy)
# lines(temp, mu.int.aegy[,1], lty=2, col=2)
# lines(temp, mu.int.aegy[,2], lty=2, col=2)

lf.aegy.round = lf.data.aegy
lf.aegy.round$T = round(lf.aegy.round$T)
lf.data.b = ddply(lf.aegy.round, .(T), trait.fun)
plot(mean~T, data=lf.data.b, main="Adult Lifespan", xlim=c(10,40), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Lifespan (days)", pch = 16, ylim = c(10, 38))
for (i in 1:nrow(lf.data.b)) {lines(rep(lf.data.b$T[i],2), c(lf.data.b$mean[i] - lf.data.b$se[i], lf.data.b$mean[i] + lf.data.b$se[i]))}
lines(temp, lf.M.aegy)
lines(temp, lf.int.aegy[,1], lty=2, col=2)
lines(temp, lf.int.aegy[,2], lty=2, col=2)

b.data.b = ddply(b.data.aegy, .(T), trait.fun)
plot(mean~T, data=b.data.b, main="Transmission Probability", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability", pch = 16)
for (i in 1:nrow(b.data.b)) {lines(rep(b.data.b$T[i],2), c(b.data.b$mean[i] - b.data.b$se[i], b.data.b$mean[i] + b.data.b$se[i]))}
lines(temp, b.M.aegy)
lines(temp, b.int.aegy[,1], lty=2, col=2)
lines(temp, b.int.aegy[,2], lty=2, col=2)

c.data.b = ddply(c.data.aegy, .(T), trait.fun)
plot(mean~T, data=c.data.b, main="Infection Probability", xlim=c(10,40), ylim=c(0,1), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Probability", pch = 16)
for (i in 1:nrow(c.data.b)) {lines(rep(c.data.b$T[i],2), c(c.data.b$mean[i] - c.data.b$se[i], c.data.b$mean[i] + c.data.b$se[i]))}
lines(temp, c.M.aegy)
lines(temp, c.int.aegy[,1], lty=2, col=2)
lines(temp, c.int.aegy[,2], lty=2, col=2)

PDR.aegy.round = PDR.data.aegy
PDR.aegy.round$T = round(PDR.aegy.round$T)
PDR.data.b = ddply(PDR.aegy.round, .(T), trait.fun)
plot(mean~T, data=PDR.data.b, main="Extrinsic Incubation Rate", xlim=c(10,40), ylim=c(0,0.25), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/day)", pch = 16)
for (i in 1:nrow(PDR.data.b)) {lines(rep(PDR.data.b$T[i],2), c(PDR.data.b$mean[i] - PDR.data.b$se[i], PDR.data.b$mean[i] + PDR.data.b$se[i]))}
lines(temp, PDR.M.aegy)
lines(temp, PDR.int.aegy[,1], lty=2, col=2)
lines(temp, PDR.int.aegy[,2], lty=2, col=2)

par(mfrow=c(1,1))

dev.off()



### Create data source and fit table
aegy.traits = c("a", "EFD", "e2a", "MDR", "lf", "b", "c", "PDR")
aegy.refs = c(paste(unique(a.data.aegy$ref), collapse=", "), paste(unique(EFD.data.aegy$ref), collapse=", "), paste(unique(e2a.data.aegy$ref), collapse=", "), paste(unique(MDR.data.aegy$ref), collapse=", "), paste(unique(lf.data.aegy$ref), collapse=", "), paste(unique(b.data.aegy$ref), collapse = ", "), paste(unique(c.data.aegy$ref), collapse = ", "), paste(unique(PDR.data.aegy$ref), collapse = ", "))
### summarize the posterior samples for each trait fit
fit.fun = function(df){
  mean = colMeans(df)[1:3]
  tmp = HPDinterval(mcmc(df))
  lower = tmp[1:3,1]
  upper = tmp[1:3,2]
  out = c(mean, lower, upper)
  names(out) = c("T0.mean", "Tm.mean", "c.mean", "T0.lower", "Tm.lower", "c.lower", "T0.upper", "Tm.upper", "c.upper")
  out
}
fits = t(sapply(list(a.samps.aegy, EFD.samps.aegy, e2a.samps.aegy, MDR.samps.aegy, lf.samps.aegy, b.samps.aegy, c.samps.aegy, PDR.samps.aegy), fit.fun))
aegy.table = data.frame('trait' = aegy.traits, 'refs' = aegy.refs, fits)
write.csv(aegy.table, file = "Aegypti_trait_data_fits.csv", row.names = F)

# pdf("Aegypti_R0_with_DTR8-informative.pdf", width = 8.5, height = 7)
aegy.R0.scale = max(aegy.R0.M)
par(mfrow = c(1,1), mar=c(6, 4, 4, 5))
plot(temp, aegy.R0.M/aegy.R0.scale, ylim=c(0, 1.05), xlim = c(15, 35), type="l", col="dodgerblue", lwd=3, xlab="Temperature (C)", ylab=expression(paste("Relative ", R[0], " (line)", sep="")))

lines(aegy.temps.DTR8, aegy.R0DTR8.M/max(aegy.R0DTR8.M), col = "gray40", lwd = 3)
legend('topleft', c("DTR 0 C", "DTR 8 C"), col = c("dodgerblue", "gray40"), lty = 1, bty = 'n', lwd = 3)

# dev.off()

### Start plotting everything

bin.fun = function(df, tem, inc, scale, by = 10, col = 1, pch = 16){
  maxes = c()
  # in a sequence of data, find the order by temperature
  ord = order(df[,tem], na.last=NA)
  # reorder the data frame
  df2 = df[ord,]
  # keep track of the original positions
  df2$ord = ord
  # find the max within each bin of 10
  n = length(ord)
  bins = c(seq(1, n, by = by), n+1)
  for (i in 2:length(bins)){
    tmp = df2[c(bins[i-1]:(bins[i]-1)), ]
    tmp2 =  subset(tmp, inc==max(inc, na.rm=T))
    maxes[i] = tmp2$ord[1]
  }
  maxes
}


# Load R0 functions from two previous papers
source("Morin2015_temp_functions.R")
source("Wesolowski_etal_2015_temp_functions.R")
source("liu-helmersson.R")
source("johansson.R")  
caminade = read.csv("caminade.csv", header = T)

# pdf("Albo_Aegy_R0_alone-informative.pdf", width = 8.5, height = 7)
pdf("Albo_Aegy_R0_alone.pdf", width = 8.5, height = 7)
plot(temp, aegy.R0.M/aegy.R0.scale, type="l", col="darkblue", lwd=3, xlab="Temperature (C)", ylab=expression(paste("Relative ", R[0], " (lines)", sep="")), ylim=c(0, 1.05), xlim=c(15, 36))
lines(temp, albo.R0.M/albo.R0.scale, col="dodgerblue", lwd=3)
legend('topleft', c(expression(italic('Ae. albopictus')), expression(italic('Ae. aegypti'))), col = c("dodgerblue", "darkblue"), lty = 1, lwd=3, bty = 'n')
dev.off()

# add 95% CI to the R0 fit

aegy.R0.scale<-max(aegy.R0.M)
aegy.R0.q2<-temp.sim.quants(aegy.R0, length(temp))##/aegy.R0.scale
albo.R0.scale<-max(albo.R0.M)
albo.R0.q2<-temp.sim.quants(albo.R0, length(temp))##/albo.R0.scale


pdf("Albo_Aegy_R0_CI.pdf", width = 8.5, height = 7)
plot(temp, aegy.R0.M/aegy.R0.scale, type="l", col="darkblue", lwd=3, xlab="Temperature (C)", ylab=expression(paste("Relative ", R[0], " (lines)", sep="")), ylim=c(0, 1.05), xlim=c(15, 36))
lines(temp, albo.R0.M/albo.R0.scale, col="dodgerblue", lwd=3)
lines(temp, albo.R0.q2[1,]/max(albo.R0.q2[1,]), col="dodgerblue", lty = 2, lwd = 1)
lines(temp, albo.R0.q2[2,]/max(albo.R0.q2[2,]), col="dodgerblue", lty = 2, lwd = 1)
lines(temp, aegy.R0.q2[1,]/max(aegy.R0.q2[1,]), col="darkblue", lty = 2, lwd = 1)
lines(temp, aegy.R0.q2[2,]/max(aegy.R0.q2[2,]), col="darkblue", lty = 2, lwd = 1)
legend('topleft', c(expression(italic('Ae. albopictus')), expression(italic('Ae. aegypti'))), col = c("dodgerblue", "darkblue"), lty = 1, lwd=3, bty = 'n')
dev.off()

pdf("Albo_Aegy_R0_CI_hists.pdf", width = 8.5, height = 8.5)
layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow=T))
plot(temp, aegy.R0.M/aegy.R0.scale, type="l", col="darkblue", lwd=3, xlab="Temperature (C)", ylab=expression(paste("Relative ", R[0], " (lines)", sep="")), ylim=c(0, 1.05), xlim=c(15, 36))
lines(temp, albo.R0.M/albo.R0.scale, col="dodgerblue", lwd=3)
lines(temp, albo.R0.q2[1,]/max(albo.R0.q2[1,]), col="dodgerblue", lty = 2, lwd = 1)
lines(temp, albo.R0.q2[2,]/max(albo.R0.q2[2,]), col="dodgerblue", lty = 2, lwd = 1)
lines(temp, aegy.R0.q2[1,]/max(aegy.R0.q2[1,]), col="darkblue", lty = 2, lwd = 1)
lines(temp, aegy.R0.q2[2,]/max(aegy.R0.q2[2,]), col="darkblue", lty = 2, lwd = 1)
legend('topleft', c(expression(italic('Ae. albopictus')), expression(italic('Ae. aegypti'))), col = c("dodgerblue", "darkblue"), lty = 1, lwd=3, bty = 'n')

hist(albo.R0.min, xlab=expression(paste("Temp. of min. ", R[0], sep="")), freq=FALSE, main="", border = "dodgerblue", col = rgb(30/255,144/255,255/255,0.3))
hist(aegy.R0.min, xlab=expression(paste("Temp. of min. ", R[0], sep="")), freq=FALSE, main="", border = "darkblue", col = rgb(0,0,1,0.3), add = T)
hist(aegy.R0.peak, xlab=expression(paste("Temp. of peak ", R[0], sep="")), freq=FALSE, main="", border = "darkblue", col = rgb(0,0,1,0.3), add = F, xlim = c(24, 31), breaks = 8)
hist(albo.R0.peak, freq=FALSE, main="", add = T, border = "dodgerblue", col = rgb(30/255,144/255,255/255,0.3))
hist(aegy.R0.max, xlab=expression(paste("Temp. of max. ", R[0], sep="")), freq=FALSE, main="", border = "darkblue", col = rgb(0,0,1,0.3), add = F, xlim = c(28, 37), breaks = 8)
hist(albo.R0.max, freq=FALSE, main="", border = "dodgerblue", col = rgb(30/255,144/255,255/255,0.3), add = T)
legend('topleft', c(expression(italic('Ae. albopictus')), expression(italic('Ae. aegypti'))), fill = c(rgb(30/255,144/255,255/255,0.3), rgb(0,0,1,0.3)), bty = 'n', border = c("dodgerblue", "darkblue"))

dev.off()


pdf("Albo_Aegy_R0_CI_poster.pdf", width = 10, height = 8)
par(mar=c(5.1,5.1,4.1,2.1))
plot(temp, aegy.R0.M/aegy.R0.scale, type="l", col="darkblue", lwd=3, xlab="Temperature (C)", ylab=expression(paste("Relative ", R[0], sep="")), ylim=c(0, 1.05), xlim=c(15, 36), cex.axis = 1.5, cex.lab = 1.5)
lines(temp, albo.R0.M/albo.R0.scale, col="dodgerblue", lwd=3)
lines(temp, albo.R0.q2[1,]/max(albo.R0.q2[1,]), col="dodgerblue", lty = 2, lwd = 2)
lines(temp, albo.R0.q2[2,]/max(albo.R0.q2[2,]), col="dodgerblue", lty = 2, lwd = 2)
lines(temp, aegy.R0.q2[1,]/max(aegy.R0.q2[1,]), col="darkblue", lty = 2, lwd = 2)
lines(temp, aegy.R0.q2[2,]/max(aegy.R0.q2[2,]), col="darkblue", lty = 2, lwd = 2)
legend('topleft', c(expression(italic('Ae. albopictus')), expression(italic('Ae. aegypti'))), col = c("dodgerblue", "darkblue"), lty = 1, lwd=3, bty = 'n', cex = 1.5)
dev.off()

pdf("Albo_Aegy_R0_Morin_Wesolowski_Liu-Helmersson_Johansson_Caminade-informative.pdf", width = 8.5, height = 7)
plot(temp, aegy.R0.M/aegy.R0.scale, type="l", col=1, lwd=2, xlab="Temperature (C)", ylab=expression(paste("Relative ", R[0], " (lines)", sep="")), ylim=c(0, 1.05), xlim=c(15, 40))
lines(temp, albo.R0.M/albo.R0.scale, col=1, lwd=2, lty = 2)
cols = brewer.pal(6, "Spectral")
morin.R0 = morin(temp)
lines(temp, morin.R0/max(morin.R0), col=cols[1], lwd=2, lty=1)
wesolowski.R0 = wesolowski(temp)
lines(temp, wesolowski.R0/max(wesolowski.R0), col=cols[2], lwd=2, lty=1)
liu_helmersson.R0 = liu_helmersson(temp)
lines(temp, liu_helmersson.R0/max(liu_helmersson.R0), col=cols[3], lwd=2, lty=1)
johansson.R0 = johansson(temp)
lines(temp, johansson.R0/max(johansson.R0), col=cols[4], lwd=2, lty=1)
c.ae = subset(caminade, mosquito=="Aeae")
c.al = subset(caminade, mosquito=="Aealbo")
lines(c.ae$T, c.ae$R0/max(c.ae$R0), col = cols[5], lwd = 2, lty = 1)
lines(c.al$T, c.al$R0/max(c.al$R0), col = cols[6], lwd = 2, lty = 2)
lines(temp, aegy.R0.M/aegy.R0.scale, col=1, lwd=2)
legend('topleft', c(expression(italic('Ae. aegypti')), expression(italic('Ae. albopictus')), 'Morin', 'Wesolowski', 'Liu-Helmersson', 'Johansson', 'Caminade aegy.', 'Caminade albo.'), col = c(1,1,cols), lty = c(1, 2, rep(1,5),2), lwd=2, bty = 'n')
dev.off()
### The shape of the Wesolowski curve depends strongly on the assumed lower bound on mu
### From their Fig. S3 it looks like they chose 0.01 as the lower bound

### Plot R0 and R0^2
pdf("Aegy_R0_squared.pdf", width = 8.5, height = 5)
layout(matrix(c(1,1,2), 1, 3, byrow=T))
plot(temp, aegy.R0.M^2+1, lty = 2, type="l", col="darkblue", lwd=2, xlab="Temperature (C)", ylab=expression(paste("log ", R[0], " (lines)", sep="")), xlim=c(10, 40), yaxt = 'n', log = "y")
lines(temp, aegy.R0.M+1, col = "darkblue", lwd = 2, lty = 1)

minT = temp[which(aegy.R0.M^2>0)[1]]
maxT = temp[tail(which(aegy.R0.M^2>0), 1)]
pkT = temp[which(aegy.R0.M^2==max(aegy.R0.M^2))]
abline(v = minT, col = "gray", lwd = 1)
abline(v = maxT, col = "gray", lwd = 1)
abline(v = pkT, col = "gray", lwd = 1)
abline(h = 2, col = "gray", lwd = 1)
legend(12, 400, legend = c(expression("R"[0]), expression("R"[0]^2)), col = "darkblue", lty = c(1,2), bty = 'n', lwd = 2, cex = 1.2)

plot(aegy.R0.M, aegy.R0.M^2, xlab = expression("R"[0]), ylab = expression("square of R"[0]), type = "l")
dev.off()

###################################
### Create an Aedes spp. prior data and fit table
### Create data source and fit table
# Subsetting the data.
data.aedes <- read.csv("Aedes_prior_data.csv", header=TRUE)
lambc = subset(data.aegy, ref=="Lambrects_et_al_2011_PNAS")
data.aedes = subset(data.aedes, ref!=as.character("Focks_Barrera_2006_Research&TrainingTropicalDis_Geneva_Paper"))

a.data.aedes = subset(data.aedes, trait.name=="a")
a.data.aedes = subset(a.data.aedes, ref=="Callado (2002)")
b.data.aedes = subset(lambc, trait.name=="b")
c.data.aedes = subset(lambc, trait.name=="c")
EFD.data.aedes = subset(data.aedes, trait.name=="EFD")
EFD.data.aedes = subset(EFD.data.aedes, ref=="Joshi_1996")
e2a.data.aedes = subset(data.aedes, trait.name=="e2a")
data1 <- data.aedes[ which(data.aedes$trait.name=="mdr"),]
data2 <- data.aedes[ which(data.aedes$trait.name=="1/MDR"),]
data2$trait <- 1/data2$trait
MDR.data.aedes = rbind(data1, data2)
PDR.data.aedes <- read.csv("EIP_priors_2015-12-04.csv", header=T)
lf.data.aedes <- data.aedes[which(data.aedes$trait.name=="1/mu"),]

load("aedes_prior_samps.Rsave")
a.samps.aedes = a.samps
b.samps.aedes = b.samps
c.samps.aedes = c.samps
EFD.samps.aedes = EFD.samps
MDR.samps.aedes = MDR.samps
e2a.samps.aedes = e2a.samps
lf.samps.aedes = lf.samps
PDR.samps.aedes = PDR.samps


aedes.traits = c("a", "EFD", "e2a", "MDR", "lf", "b", "c", "PDR")
aedes.refs = c(paste(unique(a.data.aedes$ref), collapse=", "), paste(unique(EFD.data.aedes$ref), collapse=", "), paste(unique(e2a.data.aedes$ref), collapse=", "), paste(unique(MDR.data.aedes$ref), collapse=", "), paste(unique(lf.data.aedes$ref), collapse=", "), paste(unique(b.data.aedes$ref), collapse = ", "), paste(unique(c.data.aedes$ref), collapse = ", "), paste(unique(PDR.data.aedes$ref), collapse = ", "))
### summarize the posterior samples for each trait fit
fit.fun = function(df){
  mean = colMeans(df)[1:3]
  tmp = HPDinterval(mcmc(df))
  lower = tmp[1:3,1]
  upper = tmp[1:3,2]
  out = c(mean, lower, upper)
  names(out) = c("T0.mean", "Tm.mean", "c.mean", "T0.lower", "Tm.lower", "c.lower", "T0.upper", "Tm.upper", "c.upper")
  out
}
fits = t(sapply(list(a.samps.aedes, EFD.samps.aedes, e2a.samps.aedes, MDR.samps.aedes, lf.samps.aedes, b.samps.aedes, c.samps.aedes, PDR.samps.aedes), fit.fun))
aedes.table = data.frame('trait' = aedes.traits, 'refs' = aedes.refs, fits)
write.csv(aedes.table, file = "aedes_prior_trait_data_fits.csv", row.names = F)


