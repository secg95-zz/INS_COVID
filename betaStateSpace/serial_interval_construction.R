weightConstructionIncubacionInfeccion=function(lags=7, parameterINC1=3.169434, parameterINC2=5.163921, parameterINF1=24.206087, parameterINF2=2.984198 ){
  #Distribution Incubation: Weibull
  #Parameter INC 1: Shape Weibull
  #Parameter INC 2: Scale Weibull
  #Distribution Infection: Gamma
  #Parameter INF 1: Shape Gamma
  #Parameter INF 2: rate Gamma  

  prob=diff(pweibull(1:(lags+1),parameterINC1, parameterINC2 ))
  omegaLag=NULL
  for(Om in 1:lags){
    t=sum(prob[1:Om]*rev(pgamma(1:Om, parameterINF1, parameterINF2, lower.tail = FALSE)))
    omegaLag=c(omegaLag,t)
  }
  omegaLag=(parameterINF1/parameterINF2)*omegaLag/sum(omegaLag)
  omegaLag=rev(omegaLag)
  return(omegaLag)
}
weightConstructionIncubacion=function(lags=7, parameterINC1=3.169434, parameterINC2=5.163921 ){
  #Distribution Incubation: Weibull
  #Parameter INC 1: Shape Weibull
  #Parameter INC 2: Scale Weibull
  #Distribution Infection: Gamma
  #Parameter INF 1: Shape Gamma
  #Parameter INF 2: rate Gamma  
  
  omegaLag=diff(pweibull(1:(lags+1),parameterINC1, parameterINC2 ))
  omegaLag=omegaLag/sum(omegaLag)
  return(omegaLag)
}