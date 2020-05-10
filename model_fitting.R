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
  expected_I = matrix(0, 1, moments)
  expected_N = matrix(0, 1, moments)
  N_previous = N0
  for(t in 1:moments){
    expected_N[t] = N_previous * (beta[t] + pexp(1, rate=lambda, lower.tail=FALSE))
    expected_I[t] = N_previous * beta[t]
    N_previous = expected_N[t]
  }
  return(expected_I)
}

bootstrap_samples = function( expected_I, observed_I, number_samples, window_size=3){
  "
  Once the expected_I vector is computed is posible to generate the bootsrapped
  samples for  a confidence interval cosntruction.
  
  Parameters
  ----------
  Expected_I : numeric vector
    Number of sintomatic people on a given day (estimated).
  Observed_I : numeric vector
    Number of sintomatic people on a given day.
  number_samples : numeric
    number of samples.
  window_size : numeric
    smoothing factor.
  
  Returns
  -------
  standarized_residual : data_frame
    bootstrapped sample
  "
  standarized_residual <- standarized_residuals(expected_I, observed_I, window_size)
  samples = list()
  for( i in 1:number_samples)
  {
      sample <- sample(standarized_residual$SR, length(expected_I), replace = TRUE)
      sample = sample*standarized_residual$DE
      sample_I = observed_I - sample
      samples <- list.append(samples, i=sample_I)
  }

  df <- data.frame(matrix(unlist(samples), nrow=length(samples), byrow=T))
  print(df)
  return(df)
}

standarized_residuals = function(expected_I, observed_I, window_size=3){
"
  Residuals normalization for an umbaised bootstrap.
  
  Parameters
  ----------
  Expected_I : numeric vector
    Number of sintomatic people on a given day (estimated).
  Observed_I : numeric vector
    Number of sintomatic people on a given day.
  window_size : numeric
    smoothing factor.
  
  Returns
  -------
  standarized_residual : list
    SR : residuals
    DE : standard deviation
  "
  standarized_residual <- observed_I - expected_I
  DES = 1:length(expected_I)
  for( i in seq(1, length(expected_I) - (length(expected_I) %% window_size), window_size))
  {
    DE = sd(observed_I[i:(i + window_size-1)] - expected_I[i:(i + window_size-1)])
    standarized_residual[i:(i + window_size-1)] = standarized_residual[i:(i + window_size-1)]/DE 
    DES[i:(i + window_size-1)] = DE 
  }

  if (length(expected_I) %% window_size != 0)
  {
    aux = length(expected_I) %% window_size
    DE = sd (observed_I[length(expected_I) - aux:length(expected_I)] - expected_I[length(expected_I) - aux:length(expected_I)])
    DES[length(expected_I) - aux:length(expected_I)] = DE 
    standarized_residual[length(expected_I) - aux:length(expected_I)] = standarized_residual[(length(expected_I) - aux):length(expected_I)]/DE
  }
  return(list("SR"=standarized_residual, "DE"=DES))
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

fit = function(observed_I, beta0, beta_min, beta_max, lambda0, lambda_min, lambda_max, alpha) {
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
  
  Returns
  -------
  beta : numeric vector
    Expected number of cases stemming from a single person in a single day.
  "
  # set optimization parameters
  moments = length(observed_I)
  beta_min =rep(beta_min, times=moments)
  beta_max =rep(beta_max, times=moments)
  local_opts <- list( "algorithm" = "NLOPT_LN_BOBYQA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LN_BOBYQA",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 100000,
                "local_opts" = local_opts )
  # optimize over several N0
  best = list()
  best$loss = Inf
  for (N0 in 1:5) {
    loss = function(beta) {
      lambda = optimize_lambda(observed_I, beta, N0, lambda0, lambda_min,
                               lambda_max)
      expected_I = get_expected_I(beta, N0, lambda)
      loss = mean((expected_I - observed_I) ^ 4) + alpha * mean(diff(beta) ^ 2)
      return(loss)
    }
    result <- nloptr(x0=beta0, eval_f=loss, eval_grad_f=NULL, lb=beta_min,
                     ub=beta_max, opts=opts)
    if(result$objective < best$loss) {
      best$loss = result$objective
      best$N0 = N0
      best$beta = result$solution
    }
  }
  best$lambda = optimize_lambda(observed_I, best$beta, best$N0, lambda0,
                                lambda_min, lambda_max)
  return(best)
}