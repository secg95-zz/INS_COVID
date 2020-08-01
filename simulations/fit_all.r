library(RJSONIO)
library(ggplot2)

# import each model
bayesian = new.env(); source("bayesian/model_fitting.r", local=bayesian)
ss = new.env(); source("state_space/model_fitting.r", local=ss)
poisson = new.env(); source("likelihood_incubation_exp/model_fitting.r", local=poisson)
# R series comparison method
mape = function(truth, estimator) {
  mean(abs((truth - estimator) / truth), na.rm=TRUE)
}

fit_all = function(simulation, name, ignore_beta_diff=NULL, lambda) {
  out_dir = paste("simulations", name, sep="/")
  dir.create(out_dir, recursive=TRUE)

  x = seq(1, length.out = 50)
  data <- data.frame(unlist(x), Percent.Change = unlist(simulation$I))
  cols = c("dates", "I")
  colnames(data) = cols
  data$Expected <- simulation$expected_I
  data$observed <- simulation$I
  # 
  p =
    ggplot(data) +
    geom_line(data = data, aes(x = dates, y = observed, color = "I")) +
    geom_line(data = data, aes(x = dates, y = Expected, color = "E[I]")) +
    scale_color_manual(values = c(
      'I' = 'black',
      'E[I]' = 'blue'
      )) +
    xlab('t') +
    ylab('Incidencias')+
    theme(legend.title = element_blank(),legend.text=element_text(size=20), 
          axis.title=element_text(size=20), axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15), legend.position="bottom")
  # 
  png(paste(out_dir, "I.png", sep="/"), pointsize=20, width=460, height=390)
  print(p)
  dev.off()

  # assign priors by adding a small amount of noise to theoretical values
  prior_R = mean(simulation$R) + rnorm(1, sd=0.5)
  prior_rate = 1
  prior_shape = (prior_R * prior_rate) + 1
  tau1 = simulation$tau1
  tau2 = simulation$tau2
  steps = simulation$steps
  lags = 9
  # fit models
  bayesian_fit = bayesian$fit2(
    simulation$I, window=7, prior_shape=prior_shape, prior_rate=prior_rate,
    omega=simulation$omega
  )$mode
  ss_fit = ss$fit(
    simulation$I, f_inc=c(1e-30, simulation$f_inc[1:(steps - 1)]), f_inf=c(0, simulation$f_inf[1:(steps - 1)]), prior_R=prior_R, lags=lags
  )# $data$R
  poisson_fit = list()
  for (i in 1:length(lambda)) {
    poisson_fit[[i]] = poisson$fit_robust(
      simulation$I, beta_min=0.05, beta_max=2, tau1_min=1/tau1, tau1_max=1/tau1,
      tau2_min=1/tau2, tau2_max=1/tau2, N0_min=1, N0_max=5, A0_min=0, A0_max=5,
      lambda=lambda[i], ignore_beta_diff=ignore_beta_diff, use_history=FALSE, n_iter=20
    )
  }
  # store R graphical comparison

  x = seq(1, length.out = 49)
  data <- data.frame(unlist(x), Percent.Change = unlist(simulation$I[1:49]))
  cols = c("dates", "I")
  colnames(data) = cols
  data$real <- simulation$R[1:49]
  data$bayesian <- bayesian_fit[1:49]
  data$ss <- ss_fit$data$R[1:49]
  data$poisson <- poisson_fit[[2]]$R[1:49]
  # 
  p =
    ggplot(data) +
    geom_line(data = data, aes(x = dates, y = real, color = "Teórico")) +
    geom_line(data = data, aes(x = dates, y = bayesian, color = "Bayesiano")) +
    geom_line(data = data, aes(x = dates, y = ss, color = "Estado-Espacio")) +
    geom_line(data = data, aes(x = dates, y = poisson, color = "Poisson")) +
    scale_color_manual(values = c(
      'Teórico' = 'black',
      'Bayesiano' = 'blue',
      'Estado-Espacio' = 'darkgreen',
      'Poisson' = 'deeppink1'
      )) +
    xlab('t') +
    ylab('R(t)') + ylim(1,7) +
    guides(col = guide_legend(ncol = 2)) +
    theme(legend.title = element_blank(),legend.text=element_text(size=20),
          axis.title=element_text(size=20), axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15), legend.position="bottom")
  #
  png(paste(out_dir, "R.png", sep="/"), pointsize=20, width=460, height=390)
  print(p)
  dev.off()

  data$poisson_2 <- poisson_fit[[1]]$R[1:49]
  data$poisson_3 <- poisson_fit[[3]]$R[1:49]
  # 
  p =
    ggplot(data) +
    geom_line(data = data, aes(x = dates, y = real, color = "Teórico")) +
    geom_line(data = data, aes(x = dates, y = poisson, color = "Poisson")) +
    geom_line(data = data, aes(x = dates, y = poisson_2, color = "Poisson 1")) +
    geom_line(data = data, aes(x = dates, y = poisson_3, color = "Poisson 2")) +
    scale_color_manual(values = c(
      'Teórico' = 'black',
      'Poisson 1' = 'orange',
      'Poisson 2' = 'darkblue',
      'Poisson' = 'deeppink1'
      )) +
    xlab('t') +
    ylab('R(t)') +
    guides(col = guide_legend(ncol = 2)) +
    theme(legend.title = element_blank(),legend.text=element_text(size=20),
          axis.title=element_text(size=20), axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15), legend.position="bottom")
  #
  png(paste(out_dir, "R_Poisson.png", sep="/"), pointsize=20, width=460, height=390)
  print(p)
  dev.off()
  
  # store R SMAPE with respect to theoretical
  relevant = round(simulation$steps * c(0.1, 0.9))
  relevant = relevant[1]:relevant[2]
  R_mape = list(
    "bayesian"=mape(simulation$R[relevant], bayesian_fit[relevant]),
    "ss"=mape(simulation$R[relevant], ss_fit$data$R[relevant])
  )
  for (i in 1:length(lambda)) {
    R_mape[[paste("Poisson", i)]] = mape(simulation$R[relevant], poisson_fit[[i]]$R[relevant])
  }
  write(toJSON(R_mape), paste(out_dir, "R_mape.json", sep="/"))
  # store fitted models
  for (i in 1:length(lambda)) {
    write(toJSON(poisson_fit[[i]]), paste(out_dir, paste0("poisson_", i, "_fit.json"), sep="/"))
  }
  # store priors
  priors = list("prior_R"=prior_R, "prior_rate"=prior_rate,
                "prior_shape"=prior_shape, "tau1"=tau1, "tau2"=tau2)
  write(toJSON(priors), paste(out_dir, "priors.json", sep="/"))
}