library(RJSONIO)
library(Metrics)

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
    lambda=lambda, ignore_beta_diff=ignore_beta_diff, use_history=FALSE, n_iter=10
  )
  poisson2_fit = poisson$fit_robust(
    simulation$I, beta_min=0.05, beta_max=2, tau1_min=1/3, tau1_max=1/3,
    tau2_min=1/7, tau2_max=1/7, N0_min=0, N0_max=5, A0_min=0, A0_max=5,
    lambda=lambda, ignore_beta_diff=ignore_beta_diff, use_history=TRUE, n_iter=10
  )
  # store R graphical comparison
  png(paste(out_dir, "R.png", sep="/"))
  plot(simulation$R, type="l", ylab="R", xlab="t", ylim=c(0, max(c(bayesian_R, ss_R, poisson_R, poisson_R2, simulation$R), na.rm=TRUE)))
  lines(bayesian_R, col="blue")
  lines(ss_R, col="green")
  lines(poisson_R, col="red")
  lines(poisson_R2, col="purple")
  legend("topleft", legend=c("Teórico", "Bayesiano", "EE", "General", "General+"),
         col=c("black", "blue", "green", "red", "purple"), pch=20)
  dev.off()
  # store R RMSE with respect to theoretical
  rmse = list(
    "bayesian"=mean((bayesian_R - simulation$R) ^ 2, na.rm=TRUE),
    "ss"=mean((ss_R - simulation$R) ^ 2, na.rm=TRUE),
    "general"=mean((poisson_R - simulation$R) ^ 2, na.rm=TRUE),
    "general+"=mean((poisson_R2 - simulation$R) ^ 2, na.rm=TRUE)
  )
  write(toJSON(rmse), paste(out_dir, "R_rmse.json", sep="/"))
}