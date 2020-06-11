library(nloptr)
library(rlist)
library(Rlab)
library(RJSONIO)

source("likelihood_incubation_exp/model_fitting.r")
# load data and initialize parameters
load("CasosAjustados_Diagnostico.Rda")
params = list(
  "lambda" = 2 ^ 20,
  "beta_min" = 0.1,
  "beta_max" = 0.7,
  "tau1_min" = 1 / 3,
  "tau1_max" = 10,
  "tau2_min" = 1 / 40,
  "tau2_max" = 1 / 4,
  "N0_min" = 0,
  "N0_max" = 5,
  "A0_min" = 0,
  "A0_max" = 5,
  "t0" = 1
)

observed_I = data$`Estimados Totales`
observed_I = c(observed_I[1], 0, tail(observed_I, -1))
observed_I = observed_I[params$t0:length(observed_I)]

for (i in 1:1000) {
  params$beta0 = runif(length(observed_I), params$beta_min, params$beta_max)
  params$tau10 = runif(1, params$tau1_min, params$tau1_max)
  params$tau20 = runif(1, params$tau2_min, params$tau2_max)
  params$N00 = runif(1, params$N0_min, params$N0_max)
  params$A00 = runif(1, params$A0_min, params$A0_max)
  params$observed_I = observed_I
  # fit Nicolas's model
  model = fit(observed_I, beta0=params$beta0, beta_min=params$beta_min,
              beta_max=params$beta_max, tau10=params$tau10,
              tau1_min=params$tau1_min, tau1_max=params$tau1_max,
              tau20=params$tau20, tau2_min=params$tau2_min,
              tau2_max=params$tau2_max, N00=params$N00, N0_min=params$N0_min,
              N0_max=params$N0_max, A00=params$A00, A0_min=params$A0_min,
              A0_max=params$A0_max, lambda=params$lambda, ignore_beta_diff=26)
  # save fitted model and parameters
  timestamp = format(Sys.time(), "%Y%m%d%H%M")
  out_dir = paste("likelihood_incubation_exp", "tuned", timestamp, sep="/")
  dir.create(out_dir, recursive=TRUE)
  write(toJSON(params), paste(out_dir, "params.json", sep="/"))
  write(toJSON(model), paste(out_dir, "model.json", sep="/"))
  # plot reproduction number series
  R = model$beta / model$tau2
  png(paste(out_dir, "R(t).png", sep="/"))
  plot(R, xlab="t", ylab="R(t)", type="l")
  dev.off()
  # plot expected new cases with confidence intervals
  expected_I = get_expected(model$beta, model$tau1, model$tau2, model$N0,
                            model$A0)$I
  png(paste(out_dir, "I(t).png", sep="/"))
  plot(1:length(expected_I), expected_I, type="l")
  lines(1:length(observed_I), observed_I, col="red")
  dev.off()
}