"
Fit Nicolas's model to the Colombia data.
"
library(RJSONIO)
library(foreach)
library(doParallel)
source("model_fitting.r")
source("visualization.r")

params = list(
  "alpha" = 6.4e9,
  "lambda0" = 1/5.6,
  "lambda_min" = 1/40,
  "lambda_max" = 1/4,
  "beta_min" = 0.1,
  "beta_max" = 0.7
)
bootstrap_params = list(
  "number_samples" = 1000,
  "window_size" = 10,
  "confidence" = 0.95
)
# load data and initial beta
load("CasosAjustados_Diagnostico.Rda")
observed_I = data$`Estimados Totales`
observed_I = c(observed_I[1], 0, tail(observed_I, -1))
beta0 = read.table("BetaInitials.txt", quote="\"", comment.char="")$V1
params$beta0 = beta0
params$observed_I = observed_I
# fit Nicolas's model
model = fit(observed_I, beta0=beta0, beta_min=params$beta_min,
            beta_max=params$beta_max, lambda0=params$lambda0,
            lambda_min=params$lambda_min, lambda_max=params$lambda_max,
            alpha=params$alpha)
# save fitted model and parameters
timestamp = format(Sys.time(), "%Y%m%d%H%M")
out_dir = paste("tuned", timestamp, sep="/")
dir.create(out_dir, recursive=TRUE)
write(toJSON(params), paste(out_dir, "params.json", sep="/"))
write(toJSON(model), paste(out_dir, "model.json", sep="/"))
# plot reproduction number series
R = model$beta / model$lambda
png(paste(out_dir, "R(t).png", sep="/"))
plot(R, xlab="t", ylab="R(t)", type="l")
dev.off()
# plot expected new cases with confidence intervals
expected_I = get_expected_I(model$beta, model$N0, model$lambda)
b_samples = bootstrap_samples(expected_I, observed_I,
                              bootstrap_params$number_samples,
                              bootstrap_params$window_size,
                              use_wa = TRUE)
write(toJSON(bootstrap_params), paste(out_dir, "bootstrap_params.json", sep="/"))
png(paste(out_dir, "I(t).png", sep="/"))
plot_I_intervals(expected_I, observed_I, b_samples, confidence=bootstrap_params$confidence)
dev.off()
