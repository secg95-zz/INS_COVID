"
Fit Nicolas's model to the Colombia data.
"
library(RJSONIO)
source("original/model_fitting.r")
source("visualization.r")

# load data and initialize parameters
load("CasosAjustados_Diagnostico.Rda")
aux = rep(1,length(data$`Estimados Totales`))
aux[c(27)] = 0
params = list(
  "regularization_weights" = aux,
  "lambda_min" = 1/40,
  "lambda_max" = 1/4,
  "beta_min" = 0.1,
  "beta_max" = 0.7,
  "N0_min" = 0,
  "N0_max" = 5,
  "t0" = 1
)
bootstrap_params = list(
  "number_samples" = 500,
  "window_size" = 10,
  "confidence" = 0.95
)

observed_I = data$`Estimados Totales`
observed_I = c(observed_I[1], 0, tail(observed_I, -1))
observed_I = observed_I[params$t0:length(observed_I)]
# params$beta0 = read.table("BetaInitials.txt", quote="\"", comment.char="")$V1
params$beta0 = runif(length(observed_I), params$beta_min, params$beta_max)
params$lambda0 = runif(1, params$lambda_min, params$lambda_max)
params$N00 = runif(1, params$N0_min, params$N0_max)
params$observed_I = observed_I
# fit Nicolas's model
model = fit2(observed_I, beta0=params$beta0, beta_min=params$beta_min,
             beta_max=params$beta_max, lambda0=params$lambda0, N00=params$N00,
             lambda_min=params$lambda_min, lambda_max=params$lambda_max,
             N0_min=params$N0_min, N0_max=params$N0_max,regularization=6.4e7 , regularization_weights=params$regularization_weights,
             ignore_beta_diff=NULL)
# save fitted model and parameters
timestamp = format(Sys.time(), "%Y%m%d%H%M")
out_dir = paste("original", "tuned", timestamp, sep="/")
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
                              bootstrap_params$window_size)

plot_I_intervals(expected_I, observed_I, b_samples, confidence=bootstrap_params$confidence, out_dir)

