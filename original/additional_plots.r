source("original/model_fitting.r")
source("visualization.r")
library(RJSONIO)

timestamp = "202005262246"
bootstrap_params = list(
  "number_samples" = 1000,
  "window_size" = 10,
  "confidence" = 0.95
)

model = fromJSON(paste("original", "tuned", timestamp, "model.json", sep="/"))

# plot retrospective R
steps = length(model$beta)
cum_exp = pexp(0:steps, model$lambda, lower.tail=FALSE)
R = NULL
for (t in 1:steps) {
  R = c(R, sum(model$beta[1:t] * rev(cum_exp[1:t])))
}
plot(1:steps, R, main="Retrospective R", type="l")
# plot prospective R
plot(1:steps, model$beta / model$lambda, main="Prospective R", type="l")


# plot retrospective R
steps = length(model$beta)
cum_exp = pexp(0:steps, model$lambda, lower.tail=FALSE)
R = NULL
for (t in 1:49) {
  R = c(R, sum(model$beta[1:t] * rev(cum_exp[1:t])))
}
plot(1:49, R, main="Retrospective R", type="l")
# plot prospective R
plot(1:steps, model$beta / model$lambda, main="Prospective R", type="l")

load("CasosAjustados_Diagnostico.Rda")
observed_I = data$`Estimados Totales`
observed_I = c(observed_I[1], 0, tail(observed_I, -1))
expected_I = get_expected_I(model$beta, model$N0, model$lambda)
b_samples = bootstrap_samples(expected_I, observed_I,
                              bootstrap_params$number_samples,
                              bootstrap_params$window_size)
plot_I_intervals(expected_I, observed_I, b_samples, confidence=bootstrap_params$confidence)