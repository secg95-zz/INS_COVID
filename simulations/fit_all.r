library(RJSONIO)

# import each model
bayesian = new.env(); source("bayesian/model_fitting.r", local=bayesian)
ss = new.env(); source("state_space/model_fitting.r", local=ss)
poisson = new.env(); source("likelihood_incubation_exp/model_fitting.r", local=poisson)
# R series comparison method
mape = function(truth, estimator) {
  mean(abs((truth - estimator) / truth), na.rm=TRUE)
}

fit_all = function(simulation, out_dir, ignore_beta_diff=NULL, lambda) {
  dir.create(out_dir, recursive=TRUE)
  # assign priors by adding a small amount of noise to theoretical values
  prior_R = mean(simulation$R) + rnorm(1, sd=0.5)
  prior_rate = 1
  prior_shape = (prior_R * prior_rate) + 1
  tau1 = simulation$tau1
  tau2 = simulation$tau2
  steps = simulation$steps
  lags = 9
  window = 7
  n_nloptr_iter = 10000
  n_robust_iter = 20
  bayesian_fit = list(
    "prior_rate"=prior_rate, "prior_shape"=prior_shape, "window"=window,
    "omega"=simulation$omega
  )
  ss_fit = list(
    "lags"=lags, "prior_R"=prior_R, "f_inc"=c(1e-30, simulation$f_inc[1:(steps - 1)]),
    "f_inf"=c(0, simulation$f_inf[1:(steps - 1)])
  )
  poisson_fit = list()
  
  # fit models
  bayesian_fit[["R"]] = bayesian$fit2(
    simulation$I, window=bayesian_fit$window, prior_shape=bayesian_fit$prior_shape,
    prior_rate=bayesian_fit$prior_rate, omega=bayesian_fit$omega
  )$mode
  ss_fit_df = ss$fit(
    simulation$I, f_inc=ss_fit$f_inc, f_inf=ss_fit$f_inf, prior_R=ss_fit$prior_R,
    lags=ss_fit$lags
  )$data
  for (col in colnames(ss_fit_df)) {
    ss_fit[[col]] = ss_fit_df[[col]]
  }
  poisson_fit = list()
  for (i in 1:length(lambda)) {
    poisson_fit[[i]] = poisson$fit_robust(
      simulation$I, beta_min=0.05, beta_max=2, tau1_min=1/tau1, tau1_max=1/tau1,
      tau2_min=1/tau2, tau2_max=1/tau2, N0_min=1, N0_max=5, A0_min=0, A0_max=5,
      lambda=lambda[i], ignore_beta_diff=ignore_beta_diff, use_history=FALSE,
      n_nloptr_iter=n_nloptr_iter, n_robust_iter=n_robust_iter
    )
  }
  #store models
  write(toJSON(bayesian_fit), paste(out_dir, "bayesian_fit.json", sep="/"))
  write(toJSON(ss_fit), paste(out_dir, "ss_fit.json", sep="/"))
  write(toJSON(poisson_fit), paste(out_dir, "poisson_fit.json", sep="/"))
  
  # store R MAPE with respect to theoretical
  relevant = round(simulation$steps * c(0.1, 0.9))
  relevant = relevant[1]:relevant[2]
  R_mape = list(
    "bayesian"=mape(simulation$R[relevant], bayesian_fit$R[relevant]),
    "ss"=mape(simulation$R[relevant], ss_fit$R[relevant])
  )
  for (i in 1:length(lambda)) {
    R_mape[[paste("Poisson", i)]] = mape(simulation$R[relevant], poisson_fit[[i]]$R[relevant])
  }
  write(toJSON(R_mape), paste(out_dir, "R_mape.json", sep="/"))
  
  # store simulation
  write(toJSON(simulation), paste(out_dir, "simulation.json", sep="/"))
}