library(RJSONIO)
source("model_fitting.r")
source("visualization.r")

timestamp = "202005161202"
bootstrap_params = list(
  "number_samples" = 1000,
  "window_size" = 10,
  "confidence" = 0.95
)

# obtain bootsrap samples for I
model = fromJSON(paste("tuned", timestamp, "model.json", sep="/"))
load("CasosAjustados_Diagnostico.Rda")
observed_I = data$`Estimados Totales`
observed_I = c(observed_I[1], 0, tail(observed_I, -1))
expected_I = get_expected_I(model$beta, model$N0, model$lambda)
b_samples = bootstrap_samples(expected_I, observed_I,
                              bootstrap_params$number_samples,
                              bootstrap_params$window_size)

# estimate beta at each
params = fromJSON(paste("tuned", timestamp, "params.json", sep="/"))
beta0 = read.table("BetaInitials.txt", quote="\"", comment.char="")$V1
params$beta0 = beta0
out_file = paste("tuned", timestamp, "bootstrapped_R.csv", sep="/")
for (i in 1:nrow(b_samples)) {
  bootstrapped_I = as.numeric(b_samples[i,])
  boot_model = fit(bootstrapped_I, beta0=beta0, beta_min=params$beta_min,
                   beta_max=params$beta_max, lambda0=params$lambda0,
                   lambda_min=params$lambda_min, lambda_max=params$lambda_max,
                   alpha=params$alpha)
  row = data.frame(t(boot_model$beta / boot_model$lambda))
  row["lambda"] = boot_model$lambda
  write.table(row, file=out_file, append = TRUE, row.names = FALSE)
}

# plot the R series with confidence intervals
b_samples = read.table(out_file, header=TRUE)
b_samples$lambda = NULL
# plot_I_intervals(model$beta, model$beta, b_samples, confidence=bootstrap_params$confidence)
plot_I_intervals(model$beta / model$lambda, model$beta / model$lambda, b_samples, confidence=bootstrap_params$confidence)