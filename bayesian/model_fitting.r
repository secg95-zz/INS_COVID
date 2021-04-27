library(nloptr)

fit = function(cases, window=7, lags=7, omega) {
  #Contructs de R with thompson method... (Cori et al)
  priorRate=2
  priorShape=1
  initialVector=cases[1:lags]
  
  omega_matrix=matrix(0,length(cases)-lags,lags)
  for(i in 1:length(cases)-lags){
    for(j in 1:lags){
      omega_matrix[i,j] = omega(i, j)
    }
  }
  
  infectPotential_T=matrix(0,length(cases)-lags,1)
  for(i in 1:dim(infectPotential_T)[1]){
    infectPotential_T[i]=sum(cases[i:(i+lags-1)]*rev(omega_matrix[i,]))
  }
  #Seriel interval times observed cases, there are lags values
  # computing denominatortofequation 4?
  casesMatrix=matrix(0,length(infectPotential_T)-(window-1),window)
  for(i in 1:dim(casesMatrix)[1]){
    casesMatrix[i,]=infectPotential_T[i:(i+window-1)]
  }
  #doble suma ecuacion 7
  infectPotential=casesMatrix%*%matrix(1,window,1)
  
  #Generate observation Matrix
  postcasesMatrix=matrix(0,length(cases)-(window-1),window)
  for(i in 1:dim(postcasesMatrix)[1]){
    postcasesMatrix[i,]=cases[i:(i+window-1)]
  }
  
  #suma iw(t) se supone esto es Poisson
  infectObservados=postcasesMatrix%*%matrix(1,window,1)
  
  betaParameteres=cbind(priorShape, priorRate)
  
  for(i in 1:dim(casesMatrix)[1]){
    shape1=priorShape+infectObservados[i]
    rate1=priorRate+infectPotential[i]
    betaParameteres=rbind(betaParameteres,cbind(shape1, rate1))
  }
  shape1=betaParameteres[-1,1]
  rate1=betaParameteres[-1,2]
  
  lb=qgamma(0.05, shape=shape1, rate=rate1)
  ub=qgamma(0.95, shape=shape1, rate=rate1)
  beta=((shape1-1)/rate1)
  for(i in 1:window){
    lb=c(lb[1],lb)
    ub=c(ub[1],ub)
    beta=c(beta[1],beta)
  }
  
  data=cbind(lb, beta,ub)
  data=as.data.frame(data)
  data$day=1:length(lb)
  print(data$beta)
  data$R=data$beta*sum(omega_matrix[1,])
  data$Rlb=data$lb*sum(omega_matrix[1,])
  data$Rub=data$ub*sum(omega_matrix[1,]) 
  
  #R_Caso
  lb=data$lb
  ub=data$ub
  beta=data$beta
  longCaseslb=c(rep(lb[1], length=lags-1),lb)
  longCasesub=c(rep(ub[1], length=lags-1),ub)
  longCasesbeta=c(rep(beta[1], length=lags-1),beta)
  
  data$R_c=data$beta*sum(omega_matrix[1,])
  data$R_clb=data$lb*sum(omega_matrix[1,])
  data$R_cub=data$ub*sum(omega_matrix[1,])    
  
  for(i in 1:length(lb)){
    data$R_c[i]=sum(longCasesbeta[i:(i+lags-1)]*omega_matrix[1,])
    data$R_clb[i]=sum(longCaseslb[i:(i+lags-1)]*omega_matrix[1,])      
    data$R_cub[i]=sum(longCasesub[i:(i+lags-1)]*omega_matrix[1,])
  }
  return(data)
}

fit2 = function(I, window=7, prior_shape, prior_rate, omega) {
  "
  Unvectorized model fitting.
  "
  steps = length(I)
  mode = rep(NaN, window)
  lb = rep(NaN, window)
  ub = rep(NaN, window)
  for (t in window:steps) {
    shape = prior_shape + sum(I[(t - window):t])
    # solve rate iteratively
    rate = prior_rate
    for (t_prime in (t - window):t) {
      for (tau in 1:(t_prime - 1)) {
        rate = rate + (omega[tau] * I[t_prime - tau])
      }
    }
    mode = c(mode, (shape - 1)/rate)
    lb = c(lb, qgamma(0.05, shape=shape, rate=rate))
    ub = c(ub, qgamma(0.95, shape=shape, rate=rate))
  }
  R = list("mode"=mode, "lb"=lb, "ub"=lb)
  return(R)
}

fit3 = function(I, window=7, prior_shape, prior_rate, omega) {
  "
  Like fit2, but omega is a function of (t, tau).
  "
  steps = length(I)
  mode = rep(NaN, window)
  for (t in window:steps) {
    shape = prior_shape + sum(I[(t - window):t])
    # solve rate iteratively
    rate = prior_rate
    for (t_prime in (t - window):t) {
      for (tau in 1:(t_prime - 1)) {
        rate = rate + (omega(t, tau) * I[t_prime - tau])
      }
    }
    mode = c(mode, (shape - 1)/rate)
  }
  R = list("mode"=mode)
  return(R)
}