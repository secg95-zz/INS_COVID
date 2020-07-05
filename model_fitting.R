fit_likelihood_exp = function(observed_I, beta0, beta_min, beta_max, tau10, tau1_min,
               tau1_max, tau20, tau2_min, tau2_max, N00, N0_min, N0_max, A00,
               A0_min, A0_max, lambda, ignore_beta_diff) {
  "
  Fits the model to observed daily new case counts.
  
  Parameters
  ----------
  observed_I : numeric vector
    Daily observed number of new cases.
  beta0 : numeric vector
    Initial value for beta; the expected number of cases stemming from a single
    person in a single day.
  beta_min : numeric
    Beta lower bound.
  beta_max : numeric
    Beta upper bound.
  ignore_beta_diff : numeric vector
    List of beta indices for which differences should be ignored while
    calculating the loss function. This amounts to moments in time in which we
    allow the beta series to be discontinuous.
  
  Returns
  -------
  beta : numeric vector
    Expected number of cases stemming from a single person in a single day.
  "
  # concatenate initial parameters into a single vector for optimization
  x0 = c(beta0, tau10, tau20, N00, A00)
  steps = length(observed_I)
  lb = c(rep(beta_min, steps), tau1_min, tau2_min, N0_min, A0_min)
  ub = c(rep(beta_max, steps), tau1_max, tau2_max, N0_max, A0_max)
  # set optimization parameters
  opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-7, "maxeval" = 10000000)
  # define loss function
  loss = function(x) {
    # unpack values
    beta = x[1:steps]
    tau1 = x[steps + 1]
    tau2 = x[steps + 2]
    N0 = x[steps + 3]
    A0 = x[steps + 4]
    # calculate loss
    expected = get_expected(beta, tau1, tau2, N0, A0)
    regularization =  (diff(beta) / beta[1:(steps - 1)]) ^ 2
    regularization = regularization[!1:length(regularization) %in% ignore_beta_diff]
    observed_I = round(observed_I, 0)
    factSum = function(x) {
      return(sum(log(1:round(x, 0))))
    }
    loglikelihood = observed_I*log(expected$I)-expected$I-sapply(as.matrix(observed_I), FUN=factSum)
    loglikelihood[observed_I==0] = -expected$I[observed_I==0]
    return(-mean(loglikelihood) + lambda*mean(regularization))
  }
  result = nloptr(x0=x0, eval_f=loss, lb=lb, ub=ub, opts=opts)
  model = list(
    "beta" = result$solution[1:steps],
    "tau1" = result$solution[steps + 1],
    "tau2" = result$solution[steps + 2],
    "N0" = result$solution[steps + 3],
    "A0" = result$solution[steps + 4],
    "loss" = result$objective
  )
  regularization =  (diff(model$beta) / model$beta[1:(steps - 1)]) ^ 2
  regularization = regularization[!1:(steps - 1) %in% ignore_beta_diff]
  model$likelihood = (model$loss - lambda * mean(regularization)) * steps
  return(model)
}


fit_StateSpace=function(observed_I, f_inc, f_inf){
  "
  Fits the model to observed daily new case counts.
  
  Parameters
  ----------
  observed_I : numeric vector
    Daily observed number of new cases.
  f_inc : numeric vector
    incubation  probaility, a vector of size I_observed is expected.
  omega : numeric vector
   omega from the state space model can be constructed with the function 
   weightConstructionIncubacionInfeccion, a vector of size I_observed is expected
  Returns
  -------
  model : list
    model betas , R and R_c.
  "
  omega = "unos fores"
  priorMean=3/sum(omega)
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


fit_bayesian = function(I, window=7, prior_shape, prior_rate, omega) {
  "
  Unvectorized model fitting.
  "
  steps = length(I)
  mode = rep(NaN, window)
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
  }
  R = list("mode"=mode)
  return(R)
}