library(nloptr)
library(rlist)

get_expected_I = function(N0, C0, beta, gamma, lambda) {
  "
  Calculates expected new case counts for each step in the time series.
  
  Parameters
  ----------
  N0 : numeric
    Initial number of contagious patients.
  C0 : numeric
    Initial number of infected but non-contagious patients.
  beta : numeric vector
    Expected number of cases stemming from a single person in a single day.
    Lagged one step (beta[1] = beta(0) in our notation).
  gamma : numeric
    Probability that an infected, non-contagious patient will turn infectious in
    the next time step.
  lambda : numeric
    (1 / mean) number of days an individual will continue to be infectious.
  
  Returns
  -------
  expected_I : numeric vector
    Daily expected number of new cases.
  "
  steps = length(beta)
  expected_I = rep(0, steps)
  expected_N = rep(0, steps)
  expected_C = rep(0, steps)
  N_previous = N0
  C_previous = C0
  for(t in 1:steps) {
    expected_I[t] = gamma * C_previous
    expected_N[t] = N_previous * pexp(1, rate=lambda, lower.tail=FALSE) + expected_I[t]
    expected_C[t] = (1 - gamma) * C_previous + beta[t] * N_previous
    N_previous = expected_N[t]
    C_previous = expected_C[t]
  }
  return(expected_I)
}

fit = function(observed_I, beta0, beta_min, beta_max, gamma0, gamma_min,
               gamma_max, lambda0, lambda_min, lambda_max, N00, N0_min, N0_max,
               C00, C0_min, C0_max, alpha, ignore_beta_diff) {
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
  x0 = c(beta0, gamma0, lambda0, N00, C00)
  steps = length(observed_I)
  lb = c(rep(beta_min, steps), gamma_min, lambda_min, N0_min, C0_min)
  ub = c(rep(beta_max, steps), gamma_max, lambda_max, N0_max, C0_max)
  # set optimization parameters
  opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-7,
              "maxeval" = 10000000)
  # define loss function
  loss = function(x) {
    # unpack values
    beta = x[1:steps]
    gamma = x[steps + 1]
    lambda = x[steps + 2]
    N0 = x[steps + 3]
    C0 = x[steps + 4]
    # calculate loss
    expected_I = get_expected_I(N0, C0, beta, gamma, lambda)
    beta_diff = diff(beta)
    beta_diff = beta_diff[!1:length(beta_diff) %in% ignore_beta_diff]
    loss = mean((expected_I - observed_I) ^ 4) + alpha * mean(beta_diff ^ 2)
  }
  result = nloptr(x0=x0, eval_f=loss, eval_grad_f=NULL, lb=lb, ub=ub, opts=opts)
  model = list()
  model$beta = result$solution[1:steps]
  model$gamma = result$solution[steps + 1]
  model$lambda = result$solution[steps + 2]
  model$N0 = result$solution[steps + 3]
  model$C0 = result$solution[steps + 4]
  model$loss = result$objective
  return(model)
}