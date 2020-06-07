library(nloptr)
library(rlist)
library(Rlab)
library(RJSONIO)

source("likelihood_method//model_fitting.r")
source("visualization.r")
# load data and initialize parameters
load("CasosAjustados_Diagnostico.Rda")
aux = rep(1,length(data$`Estimados Totales`))
aux[26] = 0
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

LossFunction=NULL
castigos=0:80
for(C in castigos){
  print(C)
  params$beta0 = runif(length(observed_I), params$beta_min, params$beta_max)
  params$lambda0 = runif(1, params$lambda_min, params$lambda_max)
  params$N00 = runif(1, params$N0_min, params$N0_max)
  sol=fit(observed_I, params$beta0[1:length(observed_I)], params$beta_min, params$beta_max, params$lambda0, params$lambda_min,
          params$lambda_max, params$N00, params$N0_min, params$N0_max, Castigo=C, ignore_beta_diff)
  LossFunction=c(LossFunction,sol$likelihood)
  beta0=sol$beta
  timestamp = format(Sys.time(), "%Y%m%d%H%M")
  out_dir = paste("likelihood_method", "tuned", timestamp, sep="/")
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

}

#Calcular Likelihood Ratio Test

LRTest=2*(LossFunction-LossFunction[1])

prob=pchisq(LRTest,length(beta), lower.tail=TRUE)
plot(castigos,prob)
lines(castigos,castigos*0+0.95, type="l", col="red")
#Elegir Parametro de suavizamiento
C_opt=min(castigos[prob>0.95])


sol=fit(observed_I, beta0, beta_min, beta_max, lambda0, lambda_min,
        lambda_max, N00, N0_min, N0_max, Castigo=C, ignore_beta_diff)

serie=1:length(observed_I)
#Beta estimado vs. Beta teorico
plot(serie, beta, type="l")
lines(serie, sol$beta, type="l", col="red")
#Ajuste
plot(serie, observed_I, type="l")
lines(serie, get_expected(sol$beta, N0, lambda)$I, type="l", col="red")


b_samples = bootstrap_samples(get_expected(sol$beta, N0, lambda)$I, observed_I,
                              bootstrap_params$number_samples,
                              bootstrap_params$window_size)
timestamp = format(Sys.time(), "%Y%m%d%H%M")
out_dir = paste("original", "tuned", timestamp, sep="/")
dir.create(out_dir, recursive=TRUE)
plot_I_intervals(get_expected(sol$beta, N0, lambda)$I, observed_I, b_samples, confidence=bootstrap_params$confidence, out_dir)

