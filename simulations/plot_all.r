library(ggplot2)
library(RJSONIO)
library(reshape2)

plot_all = function(out_dir, simulation, bayesian_fit, ss_fit, poisson_fit) {
  # plot incidence series
  x = seq(1, length.out = length(simulation$I))
  data <- data.frame(unlist(x), Percent.Change = unlist(simulation$I))
  cols = c("dates", "I")
  colnames(data) = cols
  data$Expected <- simulation$expected_I
  data$observed <- simulation$I
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
    theme(legend.title = element_blank(),legend.text=element_text(size=15), 
          axis.title=element_text(size=15), axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12), legend.position="bottom")
  png(paste(out_dir, "I.png", sep="/"), pointsize=15, width=420, height=350)
  print(p)
  dev.off()
  
  # plot R
  aux_length =  length(simulation$I) - 1
  x = seq(1, length.out = aux_length)
  data <- data.frame(unlist(x), Percent.Change = unlist(x))
  cols = c("dates", "I")
  colnames(data) = cols
  data$real <- simulation$R[1:aux_length]
  data$bayesian <- bayesian_fit$R[1:aux_length]
  data$ss <- ss_fit$R[1:aux_length]
  data$poisson <- poisson_fit[[2]]$R[1:aux_length]
  y_upper_bound = if (max(simulation$R) < 6) 7 else 8.5
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
    ylab('R(t)') + coord_cartesian(ylim=c(0, y_upper_bound)) +
    guides(col = guide_legend(ncol = 2)) +
    theme(legend.title = element_blank(),legend.text=element_text(size=15),
          axis.title=element_text(size=15), axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12), legend.position="bottom")
  png(paste(out_dir, "R.png", sep="/"), pointsize=15, width=420, height=350)
  print(p)
  dev.off()
  
  # plot Poisson regularization variants
  data$poisson_2 <- poisson_fit[[1]]$R[1:aux_length]
  data$poisson_3 <- poisson_fit[[3]]$R[1:aux_length]
  p =
    ggplot(data) +
    geom_line(data = data, aes(x = dates, y = real, color = "Teórico")) +
    geom_line(data = data, aes(x = dates, y = poisson, color = "Poisson")) +
    geom_line(data = data, aes(x = dates, y = poisson_2, color = "Poisson 1")) +
    geom_line(data = data, aes(x = dates, y = poisson_3, color = "Poisson 2")) +
    scale_color_manual(values = c(
      'Teórico' = 'black',
      'Poisson 1' = 'darkgreen',
      'Poisson 2' = 'blue',
      'Poisson' = 'deeppink1'
    )) +
    xlab('t') +
    ylab('R(t)') +
    guides(col = guide_legend(ncol = 2)) +
    theme(legend.title = element_blank(),legend.text=element_text(size=15),
          axis.title=element_text(size=15), axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12), legend.position="bottom")
  png(paste(out_dir, "R_Poisson.png", sep="/"), pointsize=15, width=420, height=350)
  print(p)
  dev.off()
}

# read fit results and plot
mape_table = NULL
scenario_count = 1
scenario_dirs = dir("simulations", pattern="scenario[0-9]", full.names=TRUE)
for (scenario_dir in scenario_dirs) {
  # R and I line plots
  simulation = fromJSON(paste(scenario_dir, "simulation.json", sep="/"), simplify = TRUE)
  bayesian_fit = fromJSON(paste(scenario_dir, "bayesian_fit.json", sep="/"), simplify = TRUE)
  ss_fit = fromJSON(paste(scenario_dir, "ss_fit.json", sep="/"), simplify = TRUE)
  poisson_fit = fromJSON(paste(scenario_dir, "poisson_fit.json", sep="/"), simplify = TRUE)
  plot_all(scenario_dir, simulation, bayesian_fit, ss_fit, poisson_fit)
  
  # build MAPE table
  scenario_mape = fromJSON(paste(scenario_dir, "R_mape.json", sep="/"), simplify = TRUE)
  scenario_mape = scenario_mape[c("bayesian", "ss", "Poisson 2")]
  names(scenario_mape) = c("Bayesiano", "Estado-Espacio", "Poisson")
  mape_subtable = data.frame(t(scenario_mape))
  mape_subtable[["scenario"]] = toString(scenario_count)
  scenario_count = scenario_count + 1
  mape_table = if (is.null(mape_table)) mape_subtable else rbind(mape_table, mape_subtable)
}
# store raw MAPE table
write.csv(mape_table, "simulations/mape.csv", row.names=FALSE)
# plot MAPE table
mape_melt = melt(mape_table, value.name="mape")
colnames(mape_melt) = c("scenario", "model", "mape")
mape_melt$model = as.character(mape_melt$model)
mape_melt$model[mape_melt$model == "Estado.Espacio"] = "Estado-Espacio"
mape_melt$mape[mape_melt$mape > 2] = NaN
p = ggplot(mape_melt) +
    geom_col(data=mape_melt, position="dodge", aes(x=scenario, y=mape, fill=model)) +
    scale_fill_manual(values = c('blue', 'darkgreen', 'deeppink1')) +
    xlab('Escenario') +
    ylab('MAPE(R)') +
    theme(legend.title = element_blank(), legend.text=element_text(size=15),
          axis.title=element_text(size=15), axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12), legend.position="bottom")
png("simulations/mape.png", pointsize=15, width=630, height=350)
print(p)
dev.off()