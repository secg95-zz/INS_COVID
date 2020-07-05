library(nloptr)
library(rlist)
library(Rlab)
library(bssm)

fit = function(observed_I, f_inc, f_inf, prior_R){
  "
  Parameters
  ----------
  observed_I : numeric vector
    Daily observed number of new cases.
  f_inc : numeric vector
    Incubation  time distribution. Should be the same length as 'observed_I'.
  f_inf : numeric vector
    Infectious  time distribution. Should be the same length as 'observed_I'.

  Returns
  -------
  model : list
    model betas , R and R_c.
  "
  omega = NULL
  F_inf = cumsum(f_inf)
  right_tail_inf = c(1, 1 - F_inf)
  for (tau in 1:length(observed_I)) {
    omega = c(omega, 0)
    for (tau_prime in 1:tau) {
      omega[tau] = omega[tau] + f_inc[tau_prime] * right_tail_inf[tau - tau_prime + 1]
    }
  }
  priorMean= prior_R / sum(omega)
  #Dependiente
  lags = length(observed_I)
  incubPotential_T = 1:length(lags)
  
  for(i in 1:(lags-1)){
    infected_t = 0
    for(j in 1:(lags-i)){
      infected_t = infected_t + observed_I[i + j]*f_inc[j]
    }
    incubPotential_T[i] = infected_t
  }
  
  
  #Independiente
  #longCases=c(initialVector,cases)
  infectPotential_T=observed_I
  
  for(i in 2:lags){
    caseMatrix_t=0
    for(j in 1:(i-1)){
      caseMatrix_t = caseMatrix_t + observed_I[i-j]*omega[j]
    }     
    infectPotential_T[i]=caseMatrix_t
  }
  #infectPotential_T=casesMatrix_T%*%omega
  infectPotential_T=infectPotential_T[1:length(incubPotential_T)]
  if (sum(infectPotential_T[1:2]) ==0){
    infectPotential_T[1:2] = 1
  }
  
  pars=c(0,0)
  
  obj1 <- function(pars) {
    tryCatch({
      pars2=pars
      pars2=exp(pars2)
      x=log(infectPotential_T) 
      y=as.matrix(incubPotential_T)
      Z=t(cbind(as.matrix(x),1,0))
      n=length(incubPotential_T)
      dim(Z)=c(1,3,n)
      Z=as.array(Z)
      Ti=rep(c(1,0,0,0,1,0,0,1,1),n)
      dim(Ti)=c(3,3,n)
      
      R=rep(c(0,pars2[1],0,0,0,pars2[2]),n)
      dim(R)=c(3,2,n)
      
      a1=c(1,priorMean,priorMean)
      P1=diag(c(0,10^9,10^9))
      distribution="poisson"
      u=0
      model=ssm_mng(y, Z, Ti, R, a1, P1, distribution, state_names=c("Base", "Beta", "Beta_P"))
      
      objectiveFunction=-logLik(model, nsim = 1000, nsim_states = 0)
      if(objectiveFunction==-Inf){
        objectiveFunction=-objectiveFunction
      }
    }, error = function(e) {
      objectiveFunction=Inf
    } )
    return(objectiveFunction)
  }
  
  lb=c(-10,-10)
  ub=c(10,10)
  x0 = c(0,0)
  
  local_opts = list( "algorithm" = "NLOPT_LN_BOBYQA",
                     "xtol_rel" = 1.0e-7 )
  opts = list( "algorithm" = "NLOPT_LN_BOBYQA",
               "xtol_rel" = 1.0e-7,
               "maxeval" = 100000,
               "local_opts" = local_opts )
  
  res = nloptr( x0=x0, eval_f=obj1, eval_grad_f=NULL,lb=lb, ub=ub, opts=opts)
  pars=res$solution
  pars2=pars
  pars2=exp(pars2)
  x=log(infectPotential_T) 
  y=as.matrix(incubPotential_T)
  Z=t(cbind(as.matrix(x),1,0))
  n=length(incubPotential_T)
  dim(Z)=c(1,3,n)
  Z=as.array(Z)
  Ti=rep(c(1,0,0,0,1,0,0,1,1),n)
  dim(Ti)=c(3,3,n)
  
  R=rep(c(0,pars2[1],0,0,0,pars2[2]),n)
  dim(R)=c(3,2,n)
  
  a1=c(1,priorMean,priorMean)
  P1=diag(c(0,10^9,10^9))
  distribution="poisson"
  u=0
  model=ssm_mng(y, Z, Ti, R, a1, P1, distribution, state_names=c("Base", "Beta", "Beta_P"))
  betaT=sim_smoother(model, nsim = 1)[-1,2,1] 
  for(h in 2:1000){
    betaT=cbind(betaT,sim_smoother(model, nsim = 1)[-1,2,1])  
  }
  
  beta=exp(betaT)
  
  simuQL=function(simu){
    return(quantile(simu,c(0.05,0.95)))
  }
  
  data=cbind(t(apply(beta, 1, simuQL)),apply(beta, 1, mean))
  colnames(data)=c("lb","ub", "beta")
  
  # lb=qgamma(0.05, shape=shape1, rate=rate1)
  # ub=qgamma(0.95, shape=shape1, rate=rate1)
  # beta=((shape1-1)/rate1)
  # 
  # data=cbind(betaParameteres,lb, beta,ub)
  data=as.data.frame(data)
  data$R=data$beta*sum(omega)
  data$Rlb=data$lb*sum(omega)
  data$Rub=data$ub*sum(omega)
  
  
  #R_Caso
  lb=data$lb
  ub=data$ub
  beta=data$beta
  longCaseslb=c(rep(lb[1], length=lags-1),lb)
  longCasesub=c(rep(ub[1], length=lags-1),ub)
  longCasesbeta=c(rep(beta[1], length=lags-1),beta)
  
  data$R_c=data$beta*sum(omega)
  data$R_clb=data$lb*sum(omega)
  data$R_cub=data$ub*sum(omega)    
  
  for(i in 1:dim(betaT)[1]){
    for(j in 1:i){
      data$R_c[i]=data$R_c[i] + (longCasesbeta[i]*omega[j])
      data$R_clb[i]=data$R_clb[i] + (longCaseslb[i]*omega[j])      
      data$R_cub[i]= data$R_cub[i] + (longCasesub[i]*omega[j])      
    }
  }
  return(list(data=data))
}