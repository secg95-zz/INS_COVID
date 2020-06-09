library(nloptr)
library(rlist)

get_expected = function(beta, N0, lambda) {
  "
  Returns
  -------
  expected_I : numeric vector
    Daily expected number of new cases.
  expected_N : numeric vector
    Daily expected number of infectious cases.
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
  expected_N = c(N0, expected_N)
  expected = list("I"=expected_I, "N"=expected_N)
  return(expected)
}

fit = function(observed_I, beta0, beta_min, beta_max, lambda0, lambda_min,
               lambda_max, N00, N0_min, N0_max, alpha, ignore_beta_diff) {
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
  opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-7,
              "maxeval" = 10000000)
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
    return(-mean(loglikelihood)+ alpha*mean(regularization))
  }
  result = nloptr(x0=x0, eval_f=loss, eval_grad_f=NULL, lb=lb, ub=ub, opts=opts)
  model = list()
  model$beta = result$solution[1:steps]
  model$lambda = result$solution[steps + 1]
  model$N0 = result$solution[steps + 2]
  model$loss = result$objective
  return(model)
}