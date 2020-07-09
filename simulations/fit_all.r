library(RJSONIO)

# import each model
bayesian = new.env(); source("bayesian/model_fitting.r", local=bayesian)
ss = new.env(); source("state_space/model_fitting.r", local=ss)
poisson = new.env(); source("likelihood_incubation_exp/model_fitting.r", local=poisson)
# R series comparison method
smape = function(x, y) {
  mean(abs(x - y) / ((abs(x) + abs(y)) / 2), na.rm=TRUE)
}

fit_all = function(simulation, name, ignore_beta_diff=NULL, lambda1=2^10, lambda2=2^20) {
  out_dir = paste("simulations", name, sep="/")
  dir.create(out_dir, recursive=TRUE)
  png(paste(out_dir, "I.png", sep="/"))
  plot(simulation$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(simulation$I, simulation$expected_I))))
  lines(simulation$expected_I, col="blue")
  legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
  dev.off()
  # assign priors by adding a small amount of noise to theoretical values
  prior_R = mean(simulation$R) + rnorm(1, sd=0.5)
  prior_rate = 1
  prior_shape = (prior_R * prior_rate) + 1
  tau1 = simulation$tau1
  tau2 = simulation$tau2
  # fit models
  bayesian_fit = bayesian$fit2(
    simulation$I, window=7, prior_shape=prior_shape, prior_rate=prior_rate,
    omega=simulation$omega
  )$mode
  ss_fit = ss$fit(
    simulation$I, f_inc=simulation$f_inc, f_inf=simulation$f_inf, prior_R=3
  )# $data$R
  poisson_fit = poisson$fit_robust(
    simulation$I, beta_min=0.05, beta_max=2, tau1_min=1/tau1, tau1_max=1/tau1,
    tau2_min=1/tau2, tau2_max=1/tau2, N0_min=0, N0_max=5, A0_min=0, A0_max=5,
    lambda=lambda1, ignore_beta_diff=ignore_beta_diff, use_history=FALSE, n_iter=10
  )
  poisson2_fit = poisson$fit_robust(
    simulation$I, beta_min=0.05, beta_max=2, tau1_min=1/3, tau1_max=1/3,
    tau2_min=1/7, tau2_max=1/7, N0_min=0, N0_max=5, A0_min=0, A0_max=5,
    lambda=lambda2, ignore_beta_diff=ignore_beta_diff, use_history=TRUE, n_iter=10
  )
  # store R graphical comparison
  png(paste(out_dir, "R.png", sep="/"))
  plot(simulation$R, type="l", ylab="R", xlab="t", ylim=c(0, max(
    c(bayesian_fit, ss_fit$data$R, poisson_fit$R, poisson2_fit$R, simulation$R),
    na.rm=TRUE)))
  lines(bayesian_fit, col="blue")
  lines(ss_fit$data$R, col="green")
  lines(poisson_fit$R, col="red")
  lines(poisson2_fit$R, col="purple")
  legend("topleft", legend=c("Teórico", "Bayesiano", "EE", "Poisson", "Poisson 2"),
         col=c("black", "blue", "green", "red", "purple"), pch=20)
  dev.off()
  # store R SMAPE with respect to theoretical
  relevant = round(simulation$steps * c(0.1, 0.9))
  relevant = relevant[1]:relevant[2]
  R_smape = list(
    "bayesian"=smape(bayesian_fit[relevant], simulation$R[relevant]),
    "ss"=smape(ss_fit$data$R[relevant], simulation$R[relevant]),
    "poisson"=smape(poisson_fit$R[relevant], simulation$R[relevant]),
    "poisson2"=smape(poisson2_fit$R[relevant], simulation$R[relevant])
  )
  write(toJSON(R_smape), paste(out_dir, "R_smape.json", sep="/"))
  # store fitted models
  write(toJSON(poisson_fit), paste(out_dir, "poisson_fit.json", sep="/"))
  write(toJSON(poisson2_fit), paste(out_dir, "poisson2_fit.json", sep="/"))
  # store priors
  priors = list("prior_R"=prior_R, "prior_rate"=prior_rate,
                "prior_shape"=prior_shape, "tau1"=tau1, "tau2"=tau2)
  write(toJSON(priors), paste(out_dir, "priors.json", sep="/"))
}