"
run simulations using bootstrapped data.
"
source("model_fitting.R")
load("Flu2009.rda")
load("CasosAjustados_Diagnostico.Rda")
library(nloptr)
library(rlist)
library(ggplot2)
library(reshape2)

observed_I = Flu2009$incidence$I
# initialize beta
num_simulations= 3
raw_N = cumsum(observed_I)
lagged_N = lag(raw_N)
lagged_N[1] = 1
beta0 = rep(0.3, length(observed_I)) # rep(mean(observed_I / lagged_N), length(observed_I))
beta_min = min(observed_I / lagged_N)
beta_max = max(observed_I / lagged_N)
model = fit(observed_I, beta0=beta0, beta_min=beta_min, beta_max=beta_max,
                lambda0=1/6, lambda_min=1/8, lambda_max=1/4, alpha=6.4e+13)
simulados = get_expected_I( model$beta, model$N0, model$lambda)
bootstraped_betas = list()
# fit Nicolas's model
for( i in 1:num_simulations){
    # plot reproduction number series
    # R = model$beta / model$lambda
    print(i)
    bootstrapped_residual = bootstrap_samples(simulados, observed_I, 10, 3)
}    


df <- data.frame(matrix(unlist(bootstraped_betas), nrow=length(bootstraped_betas), byrow=T))
df <- t(df)
colnames(df) <- c( 1:num_simulations)
rownames(df) <- NULL
df <- melt(df ,  id.vars = "time", variable.name = 'series')
colnames(df) <- c('time', 'series', 'value')
ggplot(df, aes(time,value)) + geom_line(aes(colour = series))
