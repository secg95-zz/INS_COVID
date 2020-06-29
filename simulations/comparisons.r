source("simulations/generators.r")

# import each model
bayesian = new.env(); source("bayesian/model_fitting.r", local=bayesian)
nicolas = new.env(); source("likelihood_incubation_exp/model_fitting.r", local=nicolas)

# Bayesian model
sim1 = simulate_bayesian(steps=50, phi=2, I0=1)
par(mfrow=c(1, 2))
plot(sim1$I, type="l", col="red", ylab="Count")
lines(sim1$expected, col="blue")
lines(sim1$expected_conditional)
legend("topleft", legend=c("I", "E[I]", "E[I | previous]"), col=c("red", "blue", "black"), pch=20)
R = bayesian$fit2(I=sim1$I, window=7, prior_shape=1, prior_rate=2, omega=sim1$omega)
plot(R$mode, type="l", ylab="Posterior R Mode")

# Almost Bayesian model with time-dependent serial interval
sim2 = simulate_inseparable_beta(steps=50, phi=2, tau0=5, I0=1)
par(mfrow=c(1, 2))
plot(sim2$I, type="l", col="red", ylab="Count")
lines(sim2$expected, col="blue")
lines(sim2$expected_conditional)
legend("topleft", legend=c("I", "E[I]", "E[I | previous]"), col=c("red", "blue", "black"), pch=20)
R = bayesian$fit3(I=sim2$I, window=7, prior_shape=1, prior_rate=2, omega=sim2$omega)
plot(R$mode, type="l", ylab="Posterior R Mode")
# fit without knowledge of the change in omega
R = bayesian$fit2(I=sim2$I, window=7, prior_shape=1, prior_rate=2, omega=sim1$omega)
lines(R$mode, type="l", ylab="Posterior R Mode", col="blue")
legend("topleft", legend=c("Timedep omega", "Independent"), col=c("black", "blue"), pch=20)