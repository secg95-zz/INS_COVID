library(nloptr)
library(rlist)
source("bootstrapping.r")

get_expected_I = function(beta, N0, lambda) {
  "
  Calculates expected new case counts for each step in the time series.
  
  Parameters
  ----------
  beta : numeric vector
    Expected number of cases stemming from a single person in a single day.
    Lagged one step (beta[1] = beta(0) in our notation).
  N0 : numeric
    Initial number of cases.
  lambda : numeric
    (1 / mean) number of days an individual will continue to be infectious.
  
  Returns
  -------
  expected_I : numeric vector
    Daily expected number of new cases.
  "
  moments = length(beta)
  expected_I = rep(0, moments)
  expected_N = rep(0, moments)
  N_previous = N0
  for(t in 1:moments){
    expected_N[t] = N_previous * (beta[t] + pexp(1, rate=lambda, lower.tail=FALSE))
    expected_I[t] = N_previous * beta[t]
    N_previous = expected_N[t]
  }
  return(expected_I)
}

optimize_lambda = function(observed_I, beta, N0, lambda0, lambda_min, lambda_max) {
  "
  Find lambda that optimizes the quartic error of estimation of E[I(t)].
  
  Parameters
  ----------
  observed_I : numeric vector
    Daily observed number of new cases.
  beta : numeric vector
    Expected number of cases stemming from a single person in a single day.
    Lagged one step (beta[1] = beta(0) in our notation).
  N0 : numeric
    Initial number of cases.
  lambda0 : numeric
    Initial value for lambda.
  
  Returns
  -------
  lambda : numeric
    Optimal value for lambda.
  "
  # define objective function: mean quartic error of estimation  of new cases
  mean_quartic_error = function(lambda) {
    expected_I = get_expected_I(beta, N0, lambda)
    return(mean((observed_I - expected_I) ^ 4))
  }
  # optimize
  local_opts <- list( "algorithm" = "NLOPT_LN_BOBYQA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LN_BOBYQA",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 100000,
                "local_opts" = local_opts )
  result <- nloptr(x0=lambda0, eval_f=mean_quartic_error, eval_grad_f=NULL,
                   lb=lambda_min, ub=lambda_max, opts=opts)
  return(result$solution)
}

gen_simulated=function(beta, N0, lambda){
  casos_Serie=rep(0,length=length(beta)+40)
  casos_Serie[1:gen_casosDisc(lambda)]=N0
  casos_Nuevo=N0
  casos_SerieV=N0
  for(i in 1:length(beta)){
    nCasos=rpois(1,beta[i]*casos_Serie[i])
    casos_SerieV[i]=casos_Serie[i+1]
    for(j in 1: nCasos){
      diasInf=gen_casosDisc(lambda)
      casos_Serie[(i+1):min(i+1+diasInf,length(casos_Serie))]=casos_Serie[(i+1):min(i+1+diasInf,length(casos_Serie))]+1
    }
    casos_Nuevo=c(casos_Nuevo,nCasos)
  }
  return(casos_Nuevo)
}

gen_casosDisc=function(lambda){
  casos=0
  r=1
  while(rbern(1,pexp(1, rate=lambda, lower.tail=FALSE))==1){
    casos=casos+1
  }
  return(casos)
}  

fit_log_likelihood = function(observed_I, beta0, beta_min, beta_max, lambda0, lambda_min,
               lambda_max, N00, N0_min, N0_max, Castigo, ignore_beta_diff) {
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
  x0 = c(beta0, lambda0, N00)
  steps = length(observed_I)
  lb = c(rep(beta_min, steps), lambda_min, N0_min)
  ub = c(rep(beta_max, steps), lambda_max, N0_max)
  # set optimization parameters
  opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-7, "maxeval" = 10000000)
  # define loss function
  loss = function(x) {
    # unpack values
    beta = x[1:steps]
    lambda = x[steps + 1]
    N0 = x[steps + 2]
    # calculate loss
    expected = get_expected(beta, N0, lambda)
    beta_diff = diff(beta)
    regularization =  (beta_diff/beta[1:length(beta_diff)]) ^ 2
    regularization = regularization[!1:length(regularization) %in% ignore_beta_diff]
    observed_I=round(observed_I,0)
    factSum=function(x){
      return(sum(log(1:round(x,0))))
    }
    loglikelihood=observed_I*log(expected$I)-expected$I-sapply(as.matrix(observed_I), FUN=factSum)
    loglikelihood[observed_I==0]=-expected$I[observed_I==0]
    loss = -mean(loglikelihood)+ Castigo*mean(regularization)
    #loss=sum((observed_I-expected$I)^2)+ sum(regularization)
    
  }
  result = nloptr(x0=x0, eval_f=loss, eval_grad_f=NULL, lb=lb, ub=ub, opts=opts)
  model = list()
  model$beta = result$solution[1:steps]
  model$lambda = result$solution[steps + 1]
  model$N0 = result$solution[steps + 2]
  model$loss = result$objective
  beta=model$beta
  beta_diff = diff(beta)
  regularization =  (beta_diff/beta[1:length(beta_diff)]) ^ 2
  regularization = regularization[!1:length(regularization) %in% ignore_beta_diff]
  model$likelihood=(model$loss-Castigo*mean(regularization))*length(beta)
  return(model)
}


