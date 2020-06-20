source("simulations/generators.r")
nicolas = source("likelihood_incubation_exp/model_fitting.r")

# Bayesian model
sim = simulate_bayesian(50, 2, 1)
plot(sim$I, type="l", col="red", ylab="Count")
lines(sim$expected, col="blue")
lines(sim$expected_conditional)
legend("topleft", legend=c("I", "E[I]", "E[I | previous]"), col=c("red", "blue", "black"), pch=20)

# Almost Bayesian model with time-dependent serial interval
sim = simulate_inseparable_beta(50, 2, 10, 1)
plot(sim$I, type="l", col="red", ylab="Count")
lines(sim$expected, col="blue")
lines(sim$expected_conditional)
legend("topleft", legend=c("I", "E[I]", "E[I | previous]"), col=c("red", "blue", "black"), pch=20)
