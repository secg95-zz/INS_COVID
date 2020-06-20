source("simulations/generators.r")
source("betaStateSpace/serial_interval_construction.R")

nicolas = source("likelihood_incubation_exp/model_fitting.r")

# Bayesian model
sim = simulate_bayesian(50, 2, 1)
plot(sim$I, type="l", col="red", ylab="Count")
lines(sim$expected, col="blue")
lines(sim$expected_conditional)
legend("topleft", legend=c("I", "E[I]", "E[I | previous]"), col=c("red", "blue", "black"), pch=20)

# Almost Bayesian model with time-dependent serial interval
sim = simulate_inseparable_beta(50, 2, 10, 1)
lags = 21
initialCases = sim$I[1:lags]
base = sim$I[(lags + 1): length(sim$I)]
result = betaThompson_non_constant_SI(base=base, window=7, initialCases=initialCases, lags=lags,  update=1, parameter1=2.235567, parameter2=5.419495, omega=sim$omega)
plot(sim$I, type="l", col="red", ylab="Count")
lines(sim$expected, col="blue")
lines(sim$expected_conditional)
legend("topleft", legend=c("I", "E[I]", "E[I | previous]"), col=c("red", "blue", "black"), pch=20)
