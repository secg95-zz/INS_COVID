
weightConstructionIncubacionInfeccion=function(lags=7, par1_inf=24.206087, par2_inf=2.984198, distribucion_inf ){
  texto=paste0("p",distribucion_inf,"(1:lags, par1_inf, par2_inf, lower.tail = FALSE)")
  omegaLag=rev(eval(parse(text = texto)))
  omegaLag <- omegaLag / (sum(omegaLag) * par2_inf)
  return(omegaLag)
}

R_caso=function(betas, omegas){
  lags=length(omegas)-1
  initialCases=as.matrix(rep(betas[1], lags))
  cases=betas
  initialVector=initialCases
  
  #Generate expectation Matrix
  longCases=c(initialVector,cases)
  casesMatrix_T=matrix(0,length(longCases)-lags,lags+1)
  for(i in 1:dim(casesMatrix_T)[1]){
    casesMatrix_T[i,]=longCases[i:(i+lags)]
  }
  t=as.vector(casesMatrix_T%*%omegas)
  return(t)
}  