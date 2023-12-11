### This file will take the individual parameter samples for all the traits and use
### them to build the CHIKV R0 model. 

# Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')
setwd("code/Mordecai2017/")
# This file contains tools for analysis and visualization.
source("mcmc_utils_all.R") 

# This file contains the thermal response functions and their derivatives
source("temp_functions_all.R") 

## Loading in samples from the Rsave (CHIKV_ParameterFits.Rsave), which have the samples
## for a, b, c, MDR, TFD, e2a created in the CHIKV_IndividualParameterFitCode.R
# load("Albopictus_ParameterFits_2016-03-30.Rsave")
load("Informative_Albopictus_ParameterFits_2016-03-30.Rsave")

# load("LifespanFits_2016-03-30.Rsave")
load("LifespanFits_2016-03-30-informative.Rsave")
lf.samps <- lf.CHIKV.samps


## Load previously computed and saved model outputs to save time
# load("Albopictus_model_outputs-uninformative.Rsave")
load("Albopictus_model_outputs-informative.Rsave")

# Temperature sequence
temp = seq(5,45,by=0.1)  
t<-length(temp)

# Length of samples, assuming they're all the same length
n = dim(a.samps)[1]
# n.mu = dim(mu.samps)[1]

# Thinned samples
thinned<-seq(1, n, by = 5) 
lthin = length(thinned)

# Creating a small constant to keep denominators from being zero.
ec<-0.000001           

## Creating the function encoding the value of R0 as a function of the parameters
## Substituting EFD = TFD/gonotrophic cycle length = TFD*a

myR0<-function(a, b, c, PDR, MDR, TFD, e2a, lf){
  mu = 1/(lf + ec)
  ((a^3*b*c*(TFD*e2a*MDR/(mu)^2)*exp((-(mu)/(PDR+ec))))/(mu))^0.5
}


## The following code runs through the samples and calculates the
## posterior trajectories (i.e. as a function of temperatures) of each
## component and of R0, and sets up a matrix to store them (each
## column is a trajectory).

R0<-matrix(NA,t,lthin)
a<-b<-c<-PDR<-MDR<-TFD<-e2a<-lf<-matrix(NA,t,lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\n")
  # calculate parameter trajectories
  i<-thinned[j]

  a[,j] = briere(temp, a.samps[i,3], a.samps[i,2], a.samps[i,1])
  PDR[,j] = briere(temp, PDR.samps[i,3], PDR.samps[i,2], PDR.samps[i,1])
  MDR[,j] = briere(temp, MDR.samps[i,3], MDR.samps[i,2], MDR.samps[i,1])
  TFD[,j] = briere(temp, TFD.samps[i,3], TFD.samps[i,2], TFD.samps[i,1])
  e2a[,j] = quad.2.trunc(temp, e2a.samps[i,1], e2a.samps[i,2], e2a.samps[i,3])
  b[,j] = briere.trunc(temp, b.samps[i,3], b.samps[i,2], b.samps[i,1])
  c[,j] = briere.trunc(temp, c.samps[i,3], c.samps[i,2], c.samps[i,1])
  lf[,j] = quad.2(temp, lf.samps[i,1], lf.samps[i,2], lf.samps[i,3])
  
# Calculate Ro equation
R0[,j]<-myR0(a[,j], b[,j], c[,j], PDR[,j], MDR[,j], TFD[,j], e2a[,j], lf[,j])  

}


## Next we calculate the posterior mean trajectory of each component of
## R0 and R0 itself. These will be used as part of the uncertainty
## analysis.

a.M<-rowMeans(a)
b.M<-rowMeans(b)
c.M<-rowMeans(c)
PDR.M<-rowMeans(PDR)
MDR.M<-rowMeans(MDR)
TFD.M<-rowMeans(TFD)
e2a.M<-rowMeans(e2a)
lf.M<-rowMeans(lf)
R0.M<-rowMeans(R0)

albo = R0.M
# save(albo, file="albo_R0.Rsave")
# save(albo, file="albo_R0_informative.Rsave")

## Reality check
plot(temp,a.M)
plot(temp,b.M)
plot(temp,c.M)
plot(temp,PDR.M)
plot(temp,MDR.M)
plot(temp,TFD.M)
plot(temp,e2a.M)
plot(temp,lf.M)
plot(temp,R0.M)

# Build matrices to hold results

R0.a<-R0.b<-R0.c<-R0.TFD<-R0.e2a<-R0.MDR<-R0.lf<-R0.PDR<-matrix(NA,t,lthin)

## For uncertainty analysis: calculate posterior samples for R0 with
## all but a single component fixed the posterior mean.

for (j in 1:lthin){
  if(j%%100==0) cat("iteration ", j, "\n")
  # Calculating trajectories
  i<-thinned[j]

R0.a[,j] = myR0(a[,j], b.M, c.M, PDR.M, MDR.M, TFD.M, e2a.M, lf.M)
R0.b[,j] = myR0(a.M, b[,j], c.M, PDR.M, MDR.M, TFD.M, e2a.M, lf.M)
R0.c[,j] = myR0(a.M, b.M, c[,j], PDR.M, MDR.M, TFD.M, e2a.M, lf.M)
R0.TFD[,j] = myR0(a.M, b.M, c.M, PDR.M, MDR.M, TFD[,j], e2a.M, lf.M)
R0.e2a[,j] = myR0(a.M, b.M, c.M, PDR.M, MDR.M, TFD.M, e2a[,j], lf.M)
R0.MDR[,j] = myR0(a.M, b.M, c.M, PDR.M, MDR[,j], TFD.M, e2a.M, lf.M)
R0.lf[,j] = myR0(a.M, b.M, c.M, PDR.M, MDR.M, TFD.M, e2a.M, lf[,j])
R0.PDR[,j] = myR0(a.M, b.M, c.M, PDR[,j], MDR.M, TFD.M, e2a.M, lf.M)

}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.

R0.q<-  apply(R0, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0, 1, FUN=quantile, probs=0.025, na.rm=T)

a.q<-  apply(R0.a, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.a, 1, FUN=quantile, probs=0.025, na.rm=T)
b.q<- apply(R0.b, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.b, 1, FUN=quantile, probs=0.025, na.rm=T)
c.q<- apply(R0.c, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.c, 1, FUN=quantile, probs=0.025, na.rm=T)
TFD.q<- apply(R0.TFD, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.TFD, 1, FUN=quantile, probs=0.025, na.rm=T)
e2a.q<-apply(R0.e2a, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.e2a, 1, FUN=quantile, probs=0.025, na.rm=T)
MDR.q<-  apply(R0.MDR, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.MDR, 1, FUN=quantile, probs=0.025, na.rm=T)
# mu.q <-  apply(R0.mu, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.mu, 1, FUN=quantile, probs=0.025, na.rm=T)
lf.q <-  apply(R0.lf, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.lf, 1, FUN=quantile, probs=0.025, na.rm=T)
PDR.q<- apply(R0.PDR, 1, FUN=quantile, probs=0.925, na.rm=T)- apply(R0.PDR, 1, FUN=quantile, probs=0.025, na.rm=T)

## Next plot relative width of quantiles 

# Creating a small constant to keep denominators from being zero.
ec<-0.2 
pdf("Albopictus_uncertainty-informative.pdf")
par(mfrow=c(1,1))
plot(temp, a.q/(R0.q +ec), col=2, type="l", lwd=2,
     xlab="Temperature (C)", ylab="Relative width of quantiles", xlim=c(17,37), ylim = c(0, 1.05))
lines(temp, b.q/(R0.q +ec), col=3, lwd=2)
lines(temp, c.q/(R0.q +ec), col="orange", lwd=2)
lines(temp, TFD.q/(R0.q +ec), col=4, lwd=2)
lines(temp, e2a.q/(R0.q +ec), col=5, lwd=2)
lines(temp, MDR.q/(R0.q +ec), col=6, lwd=2)
lines(temp, lf.q/(R0.q +ec), col=7, lwd=3)
lines(temp, PDR.q/(R0.q +ec), col=8, lwd=3)

leg.text<-c("a", "b", "c", "TFD", expression(p[EA]), "MDR", "lf", "PDR")
leg.col<-c(2:3, "orange", 4:8)
legend('topright',  leg.text, col=leg.col, lwd=2, bty='n')
dev.off()


plot(temp, R0.q, col=1, type="l", lwd=2, 
     xlab="Temperature (C)", ylab="Width of quantiles", xlim=c(12,37))
lines(temp, b.q, col=3, lwd=2)
lines(temp, c.q, col="orange", lwd=2)
lines(temp, TFD.q, col=4, lwd=2)
lines(temp, e2a.q, col=5, lwd=2)
lines(temp, MDR.q, col=6, lwd=2)
lines(temp, lf.q, col=7, lwd=3)
lines(temp, PDR.q, col=8, lwd=3)
lines(temp, a.q, col=2, lwd=2)

# Adding a legend to the plot.

leg.text<-c("a", "b", "c", "TFD", "e2a", "MDR", "LF", "PDR")
leg.col<-c(2:3, "orange", 4:8)
legend('topright',  leg.text, col=leg.col, lwd=2, bty='n')


# Save everything
save(
  temp, t,
  a, PDR, MDR, TFD, e2a, b, c, lf, R0,
  a.M, PDR.M, MDR.M, TFD.M, e2a.M, b.M, c.M, lf.M, R0.M,
  R0.a, R0.PDR, R0.MDR, R0.TFD, R0.e2a, R0.b, R0.c, R0.lf,
  a.q, PDR.q, MDR.q, TFD.q, e2a.q, b.q, c.q, lf.q, R0.q,
  file = "Albopictus_model_outputs-informative.Rsave"
  )

# save(
#   temp, t,
#   a, PDR, MDR, TFD, e2a, b, c, lf, R0,
#   a.M, PDR.M, MDR.M, TFD.M, e2a.M, b.M, c.M, lf.M, R0.M,
#   R0.a, R0.PDR, R0.MDR, R0.TFD, R0.e2a, R0.b, R0.c, R0.lf,
#   a.q, PDR.q, MDR.q, TFD.q, e2a.q, b.q, c.q, lf.q, R0.q,
#   file = "Albopictus_model_outputs-uninformative.Rsave"
# )
# 
# save(
#   temp, t,
#   a.M, PDR.M, MDR.M, TFD.M, e2a.M, b.M, c.M, lf.M, R0.M,
#   file = "Albopictus_trait_means_uninformative.Rsave"
# )

save(
  temp, t,
  a.M, PDR.M, MDR.M, TFD.M, e2a.M, b.M, c.M, lf.M, R0.M,
  file = "Albopictus_trait_means_informative.Rsave"
)

