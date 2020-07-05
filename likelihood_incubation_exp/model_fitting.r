library(nloptr)

get_expected = function(beta, tau1, tau2, N0, A0) {
  "
  Calculates expected new case counts for each step in the time series.
  
  Parameters
  ----------
  N0 : numeric
    Initial number of contagious patients.
  A0 : numeric
    Initial number of asymptomatic patients.
  beta : numeric vector
    Expected number of cases stemming from a single person in a single day.
    Lagged one step (beta[1] = beta(0) in our notation).
  tau1 : numeric
    Probability that an asymptomatic patient will turn contagious in the next
    time step.
  tau : numeric
    (1 / mean) number of days an individual will continue to be infectious.
  
  Returns
  -------
  expected : list
    Daily expected number of infectious (N), new asymptomatic (A), and new
    infectious (I) patients.
  "
  steps = length(beta)
  expected_I = rep(0, steps)
  expected_N = rep(0, steps)
  expected_A = rep(0, steps)
  N_previous = N0
  A_previous = A0
  for(t in 1:steps) {
    expected_I[t] = (1 - exp(-tau1)) * A_previous
    expected_A[t] = (exp(-tau1) * A_previous) + (beta[t] * N_previous)
    expected_N[t] = (exp(-tau2) * N_previous) + expected_I[t]
    N_previous = expected_N[t]
    A_previous = expected_A[t]
  }
  return(list("N"=c(N0, expected_N), "A"=c(A0, expected_A), "I"=expected_I))
}

fit = function(observed_I, beta0, beta_min, beta_max, tau10, tau1_min,
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

fit_robust = function(observed_I, beta_min, beta_max, tau1_min, tau1_max,
                      tau2_min, tau2_max, N0_min, N0_max, A0_min, A0_max, lambda,
                      ignore_beta_diff, n_iter) {
  best_loss = Inf
  best_model = NULL
  for (i in 1:n_iter) {
    beta0 = runif(length(observed_I), beta_min, beta_max)
    tau10 = runif(1, tau1_min, tau1_max)
    tau20 = runif(1, tau2_min, tau2_max)
    N00 = runif(1, N0_min, N0_max)
    A00 = runif(1, A0_min, A0_max)
    model = fit(observed_I, beta0=beta0, beta_min=beta_min, beta_max=beta_max,
                tau10=tau10, tau1_min=tau1_min, tau1_max=tau1_max, tau20=tau20,
                tau2_min=tau2_min, tau2_max=tau2_max, N00=N00, N0_min=N0_min,
                N0_max=N0_max, A00=A00, A0_min=A0_min, A0_max=A0_max,
                lambda=lambda, ignore_beta_diff=ignore_beta_diff)
    if (model$loss < best_loss) {
      best_model = model
      best_loss = model$loss
    }
  }
  return(best_model)
}