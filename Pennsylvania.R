"
Compare the results of Nicolas's model to Thompson's in one of Thompson's
datasets.
"

source("model_fitting.R")
load("Flu2009.rda")
observed_I = Flu2009$incidence$I
# initialize beta
raw_N = cumsum(observed_I)
lagged_N = lag(raw_N)
lagged_N[1] = 1
beta0 = rep(0.3, length(observed_I)) # rep(mean(observed_I / lagged_N), length(observed_I))
beta_min = min(observed_I / lagged_N)
beta_max = max(observed_I / lagged_N)
# fit Nicolas's model at multiple levels of regularization
colors = c("black", "blue", "green", "red", "orange")
alphas = 2 ^ (1:length(colors))
for (alpha in alphas) {
  model = fit(observed_I, beta0=beta0, beta_min=beta_min, beta_max=beta_max,
              lambda0=1/6, lambda_min=1/8, lambda_max=1/4, alpha=alpha)
  # plot reproduction number series
  R = model$beta / model$lambda
  if (alpha == alphas[1]) {
    plot(R, xlab="t", ylab="R(t)", type="l")
  } else {
    lines(R, xlab="t", ylab="R(t)", type="l")
  }
}
