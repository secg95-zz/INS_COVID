source("simulations/generators.r")
source("simulations/fit_all.r")

# Scenario 1: Constant R
beta = c(rep(0.4, 50))
sim1 = simulate_impartial(steps=50, beta, tau1=3, tau2=7, I0=5)
# verify the simulation was close to expected values
plot(sim1$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(sim1$I, sim1$expected_I))))
lines(sim1$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
# fit models and store results
fit_all(sim1, "simulations/scenario1", ignore_beta_diff=NULL, lambda=c(0, 2^26, 2^27 - 2^26))

# Scenario 2: Linearly increasing R
beta = seq(0.1, 0.9, by = 0.015)
sim2 = simulate_impartial(steps=length(beta), beta, tau1=3, tau2=7, I0=1)
# verify the simulation was close to expected values
plot(sim2$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(sim2$I, sim2$expected_I))))
lines(sim2$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
# fit models and store results
fit_all(sim2, "simulations/scenario2", lambda=c(0, 2^9, 2^15))

# Scenario 3: Linearly decreasing R
beta = seq(0.9, 0.1, by = -0.015)
sim3 = simulate_impartial(steps=length(beta), beta, tau1=3, tau2=7, I0=1)
# verify the simulation was close to expected values
plot(sim3$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(sim3$I, sim3$expected_I))))
lines(sim3$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
# fit models and store results
fit_all(sim3, "simulations/scenario3", ignore_beta_diff=NULL, lambda=c(0, 2^8, 2^15))

# Scenario 4: Step-wise increasing beta
beta = c(rep(0.2, 25), rep(0.5, 25))
sim4 = simulate_impartial(steps=50, beta, tau1=3, tau2=7, I0=5)
# verify the simulation was close to expected values
plot(sim4$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(sim4$I, sim4$expected_I))))
lines(sim4$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
# fit models and store results
fit_all(sim4, "simulations/scenario4", ignore_beta_diff=25, lambda=c(0, 2^18 + 2^17, 2^30))

# Scenario 5: Step-wise decreasing beta
beta = c(rep(0.5, 25), rep(0.2, 25))
sim5 = simulate_impartial(steps=50, beta, tau1=3, tau2=7, I0=5)
# verify the simulation was close to expected values
plot(sim5$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(sim5$I, sim5$expected_I))))
lines(sim5$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
# fit models and store results
fit_all(sim5, "simulations/scenario5", ignore_beta_diff=25, lambda=c(0, 2^9, 2^20))

# Scenario 6: Large one-day pulses
beta = c(rep(0.3, 19), 1, rep(0.3, 9), 1, rep(0.3, 20))
sim6 = simulate_impartial(steps=50, beta, tau1=3, tau2=7, I0=5)
# verify the simulation was close to expected values
plot(sim6$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(sim6$I, sim6$expected_I))))
lines(sim6$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
# fit models and store results
fit_all(sim6, "simulations/scenario6", ignore_beta_diff=c(19, 20, 29, 30), lambda=c(0, 2^5, 2^15))

# Scenario 7: End of the epidemic
beta = c(rep(0.4, 30), rep(0.01, length=20))
sim7 = simulate_impartial(steps=50, beta, tau1=3, tau2=7, I0=5)
# verify the simulation was close to expected values
plot(sim7$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(sim7$I, sim7$expected_I))))
lines(sim7$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
# fit models and store results
fit_all(sim7, "simulations/scenario7", ignore_beta_diff=30, lambda=c(0, 2^10, 2^15))

# Scenario 8: Small number of cases
beta = c(rep(0.3, 50))
sim8 = simulate_impartial(steps=50, beta, tau1=3, tau2=7, I0=10)
# verify the simulation was close to expected values
plot(sim8$I, type="l", ylab="Count", xlab="t", ylim=c(0, max(c(sim8$I, sim8$expected_I))))
lines(sim8$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
# fit models and store results
fit_all(sim8, "simulations/scenario8", lambda=c(0, 2^6 + 1, 2^10))