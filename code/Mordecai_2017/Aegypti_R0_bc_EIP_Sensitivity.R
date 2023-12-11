# Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')

####
# Load previous model runs to save time
# load("aegy_model_outputs-sensitivity.Rsave")

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


## Next we set up the temperatures over which we will be evaluating
## R0, etc, as well as the thinning interval for the samples.

temp = seq(5,45,by=0.1)  ###temperature sequence
t<-length(temp)

n = dim(a.samps)[1]      ### length of samples

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


##############################
#### Create a function to calculate R0 for each parameter set
#### Then plot the results of each R0 set

R0calcfun = function(tempsPDR = temp, tempsb = temp, tempsc = temp, a.s = a.samps, PDR.s = PDR.samps, MDR.s = MDR.samps, TFD.s = TFD.samps, e2a.s = e2a.samps, b.s = b.samps, c.s = c.samps, lf.s = lf.samps){
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
    
    a[,j] = briere(temp, a.s[i,3], a.s[i,2], a.s[i,1])
    PDR[,j] = briere(tempsPDR, PDR.s[i,3], PDR.s[i,2], PDR.s[i,1])
    MDR[,j] = briere(temp, MDR.s[i,3], MDR.s[i,2], MDR.s[i,1])
    TFD[,j] = briere(temp, TFD.s[i,3], TFD.s[i,2], TFD.s[i,1])
    e2a[,j] = quad.2.trunc(temp, e2a.s[i,1], e2a.s[i,2], e2a.s[i,3])
    # For DENV, b and c are Briere trunc
    b[,j] = briere.trunc(tempsb, b.s[i,3], b.s[i,2], b.s[i,1])
    c[,j] = briere.trunc(tempsc, c.s[i,3], c.s[i,2], c.s[i,1])
    lf[,j] = quad.2(temp, lf.s[i,1], lf.s[i,2], lf.s[i,3])
    
    # Calculate Ro equation
    R0[,j]<-myR0(a[,j], b[,j], c[,j], PDR[,j], MDR[,j], TFD[,j], e2a[,j], lf[,j])  
    
  }
  list('R0' = R0, 'a' = a, 'PDR' = PDR, 'MDR' = MDR, 'TFD' = TFD, 'e2a' = e2a, 'b' = b, 'c' = c, 'lf' = lf)
}

plot.R0 = function(R0){
  # Calculate the distribution of the lower and upper limits of R0 
  # and peak R0.
  
  R0.min<-R0.max<-R0.peak<-rep(NA, length(thinned))
  
  maxes = c()
  for (i in 1:ncol(R0)) maxes[i] = max(R0[,i])
  ind = which(maxes>0)
  R0.sub = R0[,ind]
  
  # Plotting the PEAK R0 distribution.
  for(i in 1:ncol(R0.sub)){
    ww<-which(R0.sub[,i]==max(R0.sub[,i]))
    R0.peak[i]<-temp[ww[1]]
  }
  
  # Plotting the MINIMUM R0 distribution.
  for(i in 1:ncol(R0.sub)){
    ww<-which(R0.sub[,i]>0)
    R0.min[i]<-temp[ww[1]-1]
  }
  
  # Plotting the MAXIMUM R0 distribution.
  
  for(i in 1:ncol(R0.sub)){
    ww<-which(R0.sub[,i]>0)
    lw<-length(ww)
    R0.max[i]<-temp[ww[lw]+1]
  }
  
  # Plotting the mean R0 with it's quantiles, all scaled by max mean R0.
  R0.M = rowMeans(R0.sub)
  
  layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
  R0.scale<-max(R0.M)
  R0.q2<-temp.sim.quants(R0.sub, length(temp))##/R0.scale
  plot(temp, R0.M/R0.scale, type="l", col=1, lwd=3, xlim=c(10, 40), ylim=c(0, 1.5),
       ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
  add.sim.lines(temp, sim.data=NULL, q=R0.q2/R0.scale, mycol=2)
  
  hist(R0.min, xlab="Temp of min R0", freq=TRUE, main="")
  hist(R0.peak, xlab="Temp of peak R0", freq=TRUE, main="")
  hist(R0.max, xlab="Temp of max R0", freq=TRUE, main="")
  par(mfrow=c(1,1))
}

####################################
## Calculate R0 for the base model with DENV b, c, and PDR fits
R0.base = R0calcfun()
plot.R0(R0.base$R0)

## Shift b, c, PDR +3C
R0.plus3C = R0calcfun(tempsPDR = temp - 3, tempsb = temp - 3, tempsc = temp - 3)
plot.R0(R0.plus3C$R0)

## Shift b, c, PDR +5C
R0.plus5C = R0calcfun(tempsPDR = temp - 5, tempsb = temp - 5, tempsc = temp - 5)
plot.R0(R0.plus5C$R0)

## Shift b, c, PDR -5C
R0.minus5C = R0calcfun(tempsPDR = temp + 5, tempsb = temp + 5, tempsc = temp + 5)
plot.R0(R0.minus5C$R0)

## Shift b, c, PDR -3C
R0.minus3C = R0calcfun(tempsPDR = temp + 3, tempsb = temp + 3, tempsc = temp + 3)
plot.R0(R0.minus3C$R0)

## Shift b +3C
R0.bplus3C = R0calcfun(tempsb = temp - 3)
plot.R0(R0.bplus3C$R0)

## Shift c +3C
R0.cplus3C = R0calcfun(tempsc = temp - 3)
plot.R0(R0.cplus3C$R0)

## Shift PDR +3C
R0.PDRplus3C = R0calcfun(tempsPDR = temp - 3)
plot.R0(R0.PDRplus3C$R0)

## Shift b -3C
R0.bminus3C = R0calcfun(tempsb = temp + 3)
plot.R0(R0.bminus3C$R0)

## Shift c -3C
R0.cminus3C = R0calcfun(tempsc = temp + 3)
plot.R0(R0.cminus3C$R0)

## Shift PDR -3C
R0.PDRminus3C = R0calcfun(tempsPDR = temp + 3)
plot.R0(R0.PDRminus3C$R0)

## Shift b, c, PDR Tmin -1.5C, Tmax +1.5C
rangeshift = function(df, Tminshift, Tmaxshift){
  tmp = df
  tmp$T0 = df$T0 + Tminshift
  tmp$Tm = df$Tm + Tmaxshift
  tmp
}
PDR.wider = rangeshift(PDR.samps, -1.5, 1.5)
b.wider = rangeshift(b.samps, -1.5, 1.5)
c.wider = rangeshift(c.samps, -1.5, 1.5)
R0.3Cwider = R0calcfun(b.s = b.wider, c.s = c.wider, PDR.s = PDR.wider)
plot.R0(R0.3Cwider$R0)

## Shift b, c, PDR Tmin +1.5C, Tmax -1.5C
# check for problematic cases where curves are narrower than 3C
length(which(b.samps$Tm - b.samps$T0 < 3))
length(which(c.samps$Tm - c.samps$T0 < 3))
length(which(PDR.samps$Tm - PDR.samps$T0 < 3))

PDR.narrower = rangeshift(PDR.samps, 1.5, -1.5)
b.narrower = rangeshift(b.samps, 1.5, -1.5)
c.narrower = rangeshift(c.samps, 1.5, -1.5)
R0.3Cnarrower = R0calcfun(b.s = b.narrower, c.s = c.narrower, PDR.s = PDR.narrower)
plot.R0(R0.3Cnarrower$R0)


## Plot the R0 curves together
pdf("Aegypti_bc_EIP_sensitivity_curves.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
plot(temp, rowMeans(R0.base$R0)/max(rowMeans(R0.base$R0)), type = "l", col=1, lwd=3, xlim=c(15, 37), ylim=c(0, 1),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
lines(temp, rowMeans(R0.plus3C$R0)/max(rowMeans(R0.plus3C$R0)), col=2, lwd=3)
lines(temp, rowMeans(R0.minus3C$R0)/max(rowMeans(R0.minus3C$R0)), col=4, lwd=3)
lines(temp, rowMeans(R0.bplus3C$R0)/max(rowMeans(R0.bplus3C$R0)), col=2, lwd=3, lty=2)
lines(temp, rowMeans(R0.cplus3C$R0)/max(rowMeans(R0.cplus3C$R0)), col=2, lwd=3, lty=3)
lines(temp, rowMeans(R0.PDRplus3C$R0)/max(rowMeans(R0.PDRplus3C$R0)), col=2, lwd=3, lty=4)
lines(temp, rowMeans(R0.bminus3C$R0)/max(rowMeans(R0.bminus3C$R0)), col=4, lwd=3, lty=2)
lines(temp, rowMeans(R0.cminus3C$R0)/max(rowMeans(R0.cminus3C$R0)), col=4, lwd=3, lty=3)
lines(temp, rowMeans(R0.PDRminus3C$R0)/max(rowMeans(R0.PDRminus3C$R0)), col=4, lwd=3, lty=4)
lines(temp, rowMeans(R0.base$R0)/max(rowMeans(R0.base$R0)), col=1, lwd=3)
legend('topleft', bty='n', legend = c("baseline", "all +3C", "all -3C", "b +3C", "c +3C", "PDR +3C", "b -3C", "c -3C", "PDR -3C"), col = c(1, 2, 4, 2, 2, 2, 4, 4, 4), lwd=3, lty = c(1, 1, 1, 2, 3, 4, 2, 3, 4))


plot(temp, rowMeans(R0.base$R0)/max(rowMeans(R0.base$R0)), type = "l", col=1, lwd=3, xlim=c(15, 37), ylim=c(0, 1),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
lines(temp, rowMeans(R0.plus3C$R0)/max(rowMeans(R0.plus3C$R0)), col=2, lwd=3)
lines(temp, rowMeans(R0.minus3C$R0)/max(rowMeans(R0.minus3C$R0)), col=3, lwd=3)
lines(temp, rowMeans(R0.3Cwider$R0)/max(rowMeans(R0.3Cwider$R0)), col=4, lwd=3)
lines(temp, rowMeans(R0.3Cnarrower$R0)/max(rowMeans(R0.3Cnarrower$R0)), col=5, lwd=3)
lines(temp, rowMeans(R0.plus5C$R0)/max(rowMeans(R0.plus5C$R0)), col = 6, lwd = 3)
lines(temp, rowMeans(R0.minus5C$R0)/max(rowMeans(R0.minus5C$R0)), col = 7, lwd = 3)
lines(temp, rowMeans(R0.base$R0)/max(rowMeans(R0.base$R0)), col=1, lwd=3)
legend('topleft', bty='n', legend = c("baseline", "all +3C", "all -3C", "all +5C", "all -5C", "all 3C wider", "all 3C narrower"), col = c(1:3, 6, 7, 4, 5), lwd=3)
par(mfrow=c(1,1))
dev.off()

##############################
### Comparing the mean and 95% CI for T0, Tpk, and Tm across models
compare.R0 = function(R0){
  # Calculate the distribution of the lower and upper limits of R0 
  # and peak R0.
  
  R0.min<-R0.max<-R0.peak<-rep(NA, length(thinned))
  
  maxes = c()
  for (i in 1:ncol(R0)) maxes[i] = max(R0[,i])
  ind = which(maxes>0)
  R0.sub = R0[,ind]
  
  # Plotting the PEAK R0 distribution.
  for(i in 1:ncol(R0.sub)){
    ww<-which(R0.sub[,i]==max(R0.sub[,i]))
    R0.peak[i]<-temp[ww[1]]
  }
  
  # Plotting the MINIMUM R0 distribution.
  for(i in 1:ncol(R0.sub)){
    ww<-which(R0.sub[,i]>0)
    R0.min[i]<-temp[ww[1]-1]
  }
  
  # Plotting the MAXIMUM R0 distribution.
  
  for(i in 1:ncol(R0.sub)){
    ww<-which(R0.sub[,i]>0)
    lw<-length(ww)
    R0.max[i]<-temp[ww[lw]+1]
  }

  out = c(mean(R0.peak, na.rm=T), HPDinterval(mcmc(R0.peak)),
          mean(R0.min, na.rm = T), HPDinterval(mcmc(R0.min)),
          mean(R0.max, na.rm = T), HPDinterval(mcmc(R0.max))
    )
  names(out) = c("Tpk.mean", "Tpk.lower95%", "Tpk.upper95%", "Tmin.mean", "Tmin.lower95%", "Tmin.upper95%", "Tmax.mean", "Tmax.lower95%", "Tmax.upper95%")
  out
}  


results = sapply(list(R0.base$R0, R0.plus3C$R0, R0.minus3C$R0, R0.plus5C$R0, R0.minus5C$R0, R0.3Cwider$R0, R0.3Cnarrower$R0), function(df) compare.R0(df))
colnames(results) = c("base", "+3C", "-3C", "+5C", "-5C", "3C wider", "3C narrower")
results

### Summary: The most extreme changes--shifting all three curves up or down by 3C or 5C, resulted in less than one degree change in Tpk.
### Tmin shifted approximately proportionally to the shift in b, c, and PDR (i.e., Tmin increased by ~5C for the +5C model).
### Tmax was essentially unaffected by all shifts in b, c, and PDR.

##############################
### Save the results
save(R0.base, R0.plus3C, R0.minus3C, R0.plus5C, R0.minus5C, R0.3Cwider, R0.3Cnarrower, R0.bplus3C, R0.cplus3C, R0.PDRplus3C, R0.bminus3C, R0.cminus3C, R0.PDRminus3C,
     file = "aegy_model_outputs-sensitivity.Rsave")
