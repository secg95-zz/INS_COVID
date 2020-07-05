source("simulations/generators.r")

# import each model
bayesian = new.env(); source("bayesian/model_fitting.r", local=bayesian)
ss = new.env(); source("state_space/model_fitting.r", local=ss)
nicolas = new.env(); source("likelihood_incubation_exp/model_fitting.r", local=nicolas)

# scenario 1
beta = c(rep(0.3, 25), rep(0.5, 25))
out_dir = "simulations/scenario1"
dir.create(out_dir, recursive=TRUE)
# simulate and verify the simulation was close to expected values
sim1 = simulate_impartial(steps=50, beta, tau1=3, tau2=7, I0=1)
png(paste(out_dir, "I.png", sep="/"))
plot(sim1$I, type="l", ylab="Count", xlab="t")
lines(sim1$expected_I, col="blue")
legend("topleft", legend=c("I", "E[I]"), col=c("black", "blue"), pch=20)
dev.off()
# fit models
bayesian_R = bayesian$fit2(sim1$I, window=7, prior_shape=4, prior_rate=1,
                           omega=sim1$omega)$mode
ss_R = ss$fit(sim1$I, f_inc=sim1$f_inc, f_inf=sim1$f_inf, prior_R=3)
nicolas_R = nicolas$fit_robust(sim1$I, beta_min=0.05, beta_max=2, tau1_min=1 / 10,
                               tau1_max=100, tau2_min=1 / 15, tau2_max=1,
                               N0_min=0, N0_max=5, A0_min=0, A0_max=5,
                               lambda=2 ^ 10, ignore_beta_diff=25, n_iter=10)