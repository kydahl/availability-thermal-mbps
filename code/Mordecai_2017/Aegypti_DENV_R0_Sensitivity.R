### This file will take resulting analysis from the Aegypti DENV model and run 
### through a bit of sensitvity analysis. 

# Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')
setwd("code/Mordecai2017/")
# This file contains tools for analysis and visualization.
source("mcmc_utils_all.R") 

# This file contains the thermal response functions and their derivatives
source("temp_functions_all.R") 

## Loading in samples from the Rsave (DENV_ParameterFits.Rsave), which have the samples
## for a, b, c, MDR, EFD, e2a, and PDR created in the DENV_IndividualParameterFitCode.R.

# load("Aegypti_DENV_ParameterFits_2016-03-30.Rsave")
load("Informative_Aegypti_DENV_ParameterFits_2016-03-30.Rsave")

# load("LifespanFits_2016-03-30.Rsave")
load("LifespanFits_2016-03-30-informative.Rsave")
lf.samps <- lf.DENV.samps

## Load the outputs from Aegypti_DENV_R0_Analysis.R
# load("Aegypti_DENV_model_outputs-uninformative.Rsave")
# load("Aegypti_DENV_trait_means_uninformative.Rsave")
load("Aegypti_DENV_model_outputs-informative.Rsave")
load("Aegypti_DENV_trait_means_informative.Rsave")

## Load previously calculated and saved derivatives to save time
# load("aegy_derivative_outputs-uninformative.Rsave")
# load("aegy_derivative_outputs-informative.Rsave")

## Next we set up the temperatures over which we will be evaluating
## R0, etc, as well as the thinning interval for the samples.

temp = seq(5,45,by=0.1)  ###temperature sequence
t<-length(temp)

n = dim(a.samps)[1]    	### length of samples

thinned<-seq(1, n, by=5)### thinned samples

lthin<-length(thinned)  ### number of of thinned samples

ec<-0.000001            ### small constant used to keep denominators
### from being numerically zero

## Creating the function encoding the value of R0 as a function of the parameters
myR0<-function(a, b, c, PDR, MDR, EFD, e2a, lf){
  mu = 1/(lf + ec)
  bc = (b*c)
  ((a^2*bc*(EFD*e2a*MDR/(mu)^2)*exp((-mu/(PDR+ec))))/(mu))^0.5
}


## Build matrices to hold derivatives, etc.

dR0.dT= dR0.da = dR0.db = dR0.dc= dR0.dEFD= dR0.de2a= dR0.dMDR= dR0.dlf= dR0.dPDR= dR0.R0dT<-matrix(NA,t,lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\n")
  ## calculate derivative trajectories
  i<-thinned[j]
  da.dT = d_briere(temp,a.samps[i,3],a.samps[i,2],a.samps[i,1])  
  db.dT = d_briere_trunc(temp,b.samps[i,3], b.samps[i,2], b.samps[i,1])  
  dc.dT = d_briere_trunc(temp,c.samps[i,3], c.samps[i,2], c.samps[i,1]) 
  dPDR.dT = d_briere(temp,PDR.samps[i,3],PDR.samps[i,2],PDR.samps[i,1])
  dMDR.dT = d_briere(temp,MDR.samps[i,3],MDR.samps[i,2],MDR.samps[i,1])
  dEFD.dT = d_briere(temp,EFD.samps[i,3],EFD.samps[i,2],EFD.samps[i,1])
  de2a.dT = d_quad.2.trunc(temp,e2a.samps[i,1],e2a.samps[i,2],e2a.samps[i,3])
  dlf.dT = d_quad.2(temp,lf.samps[i,1],lf.samps[i,2],lf.samps[i,3])
  
  ## Calculate derivative equations
  ec<-0.000001
  
  dR0.da[,j] = 1/2*(myR0(a[,j], b.M, c.M, PDR.M, MDR.M, EFD.M, e2a.M, lf.M)*2/(a[,j]+ec) * da.dT)
  dR0.db[,j] = 1/2*(myR0(a.M, b[,j], c.M, PDR.M, MDR.M, EFD.M, e2a.M, lf.M)/(b[,j]+ec) * db.dT)
  dR0.dc[,j] = 1/2*(myR0(a.M, b.M, c[,j], PDR.M, MDR.M, EFD.M, e2a.M, lf.M)/(c[,j]+ec) * dc.dT)
  dR0.dEFD[,j] = 1/2*(myR0(a.M,  b.M, c.M, PDR.M, MDR.M, EFD[,j], e2a.M,lf.M)/(EFD[,j]+ec) * dEFD.dT)
  dR0.de2a[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR.M, MDR.M, EFD.M, e2a[,j], lf.M)/(e2a[,j]+ec) * de2a.dT)
  dR0.dMDR[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR.M, MDR[,j], EFD.M, e2a.M, lf.M)/(MDR[,j]+ec) * dMDR.dT)
### Double check that this one's right
  dR0.dlf[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR.M, MDR.M, EFD.M, e2a.M, lf[,j])*(-3*lf[,j]
                                                                          -1/(PDR.M+ec))*(-1/(lf[,j]+ec)^2) * dlf.dT)
  dR0.dPDR[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR[,j], MDR.M, EFD.M, e2a.M, lf.M)*lf.M/(PDR[,j]+ec)^2 * dPDR.dT)
  dR0.dT[,j] =  dR0.da[,j] + dR0.db[,j]+ dR0.dc[,j] + dR0.dEFD[,j] + dR0.de2a[,j] + dR0.dMDR[,j] + dR0.dlf[,j] + dR0.dPDR[,j]
}

## Relative sensitivities in R0 (different from mordecai paper, just
## normalized by max R0 for each curve) with regard to temperature,
## broken down into individual parameters' contributions:

R0.Med<-apply(R0, 1, FUN=median, na.rm=FALSE)

dR0.R0da=dR0.da/max(R0.M)
dR0.R0db=dR0.db/max(R0.M)
dR0.R0dc=dR0.dc/max(R0.M)
dR0.R0dEFD=dR0.dEFD/max(R0.M)
dR0.R0de2a=dR0.de2a/max(R0.M)
dR0.R0dMDR=dR0.dMDR/max(R0.M)
dR0.R0dlf=dR0.dlf/max(R0.M)
dR0.R0dPDR=dR0.dPDR/max(R0.M)
dR0.R0dT = dR0.dT/max(R0.M)

## Width of the quantiles of these normalized sensitivities
dR0.q<-  apply(dR0.R0dT, 1, FUN=quantile, probs=0.925, na.rm=F) - apply(dR0.R0dT, 1, FUN=quantile, probs=0.025, na.rm=F)
dR0da.q<-  apply(dR0.R0da, 1, FUN=quantile, probs=0.925)- apply(dR0.R0da, 1, FUN=quantile, probs=0.025)
dR0db.q<- apply(dR0.R0db, 1, FUN=quantile, probs=0.925)- apply(dR0.R0db, 1, FUN=quantile, probs=0.025)
dR0dc.q<- apply(dR0.R0dc, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dc, 1, FUN=quantile, probs=0.025)
dR0dEFD.q<- apply(dR0.R0dEFD, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dEFD, 1, FUN=quantile, probs=0.025)
dR0de2a.q<-apply(dR0.R0de2a, 1, FUN=quantile, probs=0.925)- apply(dR0.R0de2a, 1, FUN=quantile, probs=0.025)
dR0dMDR.q<-  apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.025)
dR0dlf.q <-  apply(dR0.R0dlf, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dlf, 1, FUN=quantile, probs=0.025)
dR0dPDR.q<- apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.025)

plot(temp, dR0.q, type="l", xlim=c(15,35), lwd=2, xlab="Temperature (C)", ylab="width of HPD interval of dR0/dT")
lines(temp, dR0da.q, col=2, lwd=2)
lines(temp, dR0db.q, col=3, lwd=2)
lines(temp, dR0dc.q, col="orange", lwd=2)
lines(temp, dR0dEFD.q, col=4, lwd=2)
lines(temp, dR0de2a.q, col=5, lwd=2)
lines(temp, dR0dMDR.q, col=6, lwd=2)
lines(temp, dR0dlf.q, col=7, lwd=2)
lines(temp, dR0dPDR.q, col=8, lwd=2)

leg.text<-c("R0", "a", "b", "c", "EFD", "e2a", "MDR", "LF", "PDR")
leg.col<-c(1:3, "orange", 4:8)
legend('topleft',  leg.text, col=leg.col, lwd=2)


## Now plot the relative width of the quantiles of the
## sensitivities. We again include a small offset to keep the
## denominator from being numerically zero. Figure 3(b)
ec=10^(-6)
plot(temp, dR0da.q/(dR0.q +ec), col=2, type="l", ylim=c(0,1), lwd=2, xlab="Temperature (C)", ylab="Relative width of HPD intervals", xlim=c(10, 40))
lines(temp, dR0db.q/(dR0.q +ec), col=3, lwd=2)
lines(temp, dR0dc.q/(dR0.q +ec), col="orange", lwd=2)
lines(temp, dR0dEFD.q/(dR0.q +ec), col=4, lwd=2)
lines(temp, dR0de2a.q/(dR0.q +ec), col=5, lwd=2)
lines(temp, dR0dMDR.q/(dR0.q +ec), col=6, lwd=2)
lines(temp, dR0dlf.q/(dR0.q +ec), col=7, lwd=3)
lines(temp, dR0dPDR.q/(dR0.q +ec), col=8, lwd=3)
## add legend
leg.text<-c("a", "b", "c", "EFD", "e2a", "MDR", "LF", "PDR")
leg.col<-c(2, 3, "orange", 4:8)
legend(9, 0.7,  leg.text, col=leg.col, lwd=2, bty='n')


## Plot the normalized sensitivities vs. temperature
pdf("Aegypti_rel_sensitivity.pdf")
plot(temp, rowMeans(dR0.R0dPDR), type="l", xlim=c(15,35), lwd=2, xlab="Temperature (C)", ylab="Normalized dR0/dX")
lines(temp, rowMeans(dR0.R0da), col=2, lwd=2)
lines(temp, rowMeans(dR0.R0db), col=3, lwd=2)
lines(temp, rowMeans(dR0.R0dc), col="orange", lwd=2)
lines(temp, rowMeans(dR0.R0dEFD), col=4, lwd=2)
lines(temp, rowMeans(dR0.R0de2a), col=5, lwd=2)
lines(temp, rowMeans(dR0.R0dMDR), col=6, lwd=2)
lines(temp, rowMeans(dR0.R0dlf), col=7, lwd=2)
lines(temp, rowMeans(dR0.R0dPDR), col=8, lwd=2)

leg.text<-c("a", "b", "c", "EFD", expression(p[EA]), "MDR", "lf", "PDR")
leg.col<-c(2:3, "orange", 4:8)
legend('topleft',  leg.text, col=leg.col, lwd=2)
dev.off()

## Plot the average parameter fits together to look at their optima
pdf("Aegypti_thermal_responses_together.pdf")
plot(temp, a.M/max(a.M), col=2, xlab = "Temperature (C)", ylab = "Standardized trait response", type="l", lwd = 2)
lines(temp, b.M/max(b.M), col=3, lwd=2)
lines(temp, c.M/max(c.M), col="orange", lwd=2)
lines(temp, EFD.M/max(EFD.M), col=4, lwd=2)
lines(temp, e2a.M/max(e2a.M), col=5, lwd=2)
lines(temp, MDR.M/max(MDR.M), col=6, lwd = 2)
lines(temp, lf.M/max(lf.M), col=7, lwd = 2)
lines(temp, PDR.M/max(PDR.M), col=8, lwd=2)
lines(temp, R0.M/max(R0.M), col=1, lwd=3)

leg.text<-c("R0", "a", "b", "c", "EFD", "e2a", "MDR", "LF", "PDR")
leg.col<-c(1:3, "orange", 4:8)
legend('topleft',  leg.text, col=leg.col, lwd=2)
dev.off()

save(
  temp, t,
  dR0.dT, dR0.da, dR0.db, dR0.dc, dR0.dEFD, dR0.de2a, dR0.dMDR, dR0.dlf, dR0.dPDR, dR0.R0dT,
  file = "aegy_derivative_outputs-informative.Rsave"
  )