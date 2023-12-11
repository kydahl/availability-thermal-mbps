# temp_functions_all.R
AIC.c<-function(myfit){##, ndata){
  k<-(length(coef(myfit))+1)
  lLfit<-logLik(myfit)
  
  n<-attr(lLfit, "nobs")
  cort<-2*k*(k+1)/(n-k-1)

  AICc <- -2*lLfit + 2*k + cort

  return(AICc)

}

linear<-function(T, inter, slope){
  x = inter+slope*T
  x[which(x<0)] <- 0
  x
}

linear.trunc<-function(T, inter, slope){
  x = inter+slope*T
  x[which(x>1)] <- 1
  x[which(x<0)] <- 0
  x
}


quad.trunc<-function(T, T0, Tm, qd){
  1*((qd*(T-T0)*(T-Tm))>1) + (qd*(T-T0)*(T-Tm)*(T<=Tm)*(T>=T0))*((qd*(T-T0)*(T-Tm))<1)
}

briere<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if(t[i]>T0 && t[i]<Tm){  b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))  }
    else {b[i]<-0}
  }
  b
}

briere.trunc<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if(t[i]>T0 && t[i]<Tm){  b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))  }
    else {b[i]<-0}
  }
  b[which(b>1)]<-1
  b
}



d_briere<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-(c*(-5*(t[i]^2) + 3*t[i]*T0 + 4*t[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-t[i]))
    )}
    else {b[i]<-0}
  }
  b
}

d_briere_trunc<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-(c*(-5*(t[i]^2) + 3*t[i]*T0 + 4*t[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-t[i]))
    )}
    else {b[i]<-0}
    if (briere.trunc(t[i], c, Tm, T0)==1) b[i] <- 0
  }
  b
}



quad<-function(t, inter, n.slope, qd, lim=0.0001){
  b=c()
  for (i in 1:length(t)){
    if (inter-n.slope*t[i]+qd*t[i]^2>lim) {b[i]<-inter-n.slope*t[i]+qd*t[i]^2}
    else {b[i]<-lim}
  }
  b
}

quad.2<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if(t[i]>T0 && t[i]<Tm) {b[i]<-qd*(t[i]-T0)*(t[i]-Tm)}
    else {b[i]<-0}
  }
  b
}

d_quad.2<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if (t[i]>T0 && t[i]<Tm) {b[i]<-qd*(2*t[i]-T0-Tm)}
    else {b[i]<-0}
  }
  b
}

quad.2.trunc<-function(t, T0, Tm, qd, lim=0.0001){
  b=c()
  for (i in 1:length(t)){
    if(t[i]>T0 && t[i]<Tm) {b[i]<-qd*(t[i]-T0)*(t[i]-Tm)}
    else {b[i]<-lim}
    if(b[i]>1){ b[i] <-1 }
  }
  b
}

d_quad.2.trunc<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if (t[i]>T0 && t[i]<Tm) {b[i]<-qd*(2*t[i]-T0-Tm)}
    else {b[i]<-0}
    if(qd*(t[i]-T0)*(t[i]-Tm)>1) b[i]<-0
  }
  b
}

