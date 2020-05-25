library(RJSONIO)
source("model_fitting.r")
source("visualization.r")

timestamp = "202005161202"
bootstrap_params = list(
  "number_samples" = 1000,
  "window_size" = 10,
  "confidence" = 0.95
)

model = fromJSON(paste("tuned", timestamp, "model.json", sep="/"))
load("CasosAjustados_Diagnostico.Rda")
observed_I = data$`Estimados Totales`
observed_I = c(observed_I[1], 0, tail(observed_I, -1))
expected_I = get_expected_I(model$beta, model$N0, model$lambda)
b_samples = bootstrap_samples(expected_I, observed_I,
                              bootstrap_params$number_samples,
                              bootstrap_params$window_size)
plot_I_intervals(expected_I, observed_I, b_samples, confidence=bootstrap_params$confidence)