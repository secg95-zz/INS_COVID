library(RJSONIO)
library(readxl)
library(ggplot2)
library(ggforce)
library(cowplot)

# Load fitted models and Poisson confidence intervals
fitted_models = fromJSON("case_study/fitted_models.json", nullValue = NA)
poisson_confidence = readRDS("case_study/df_estimaIC_Rt_procesoAnalitico_.rds")

# Format dataframes for plotting with ggplot
incidence_df = data.frame(list(
  date = as.Date(fitted_models$date_range),
  I = fitted_models$I
))
bayesian_R_df = data.frame(fitted_models$bayesian)
bayesian_R_df$model = "bayesian"
bayesian_R_df$date = incidence_df$date
statespace_R_df = data.frame(fitted_models$state_space)
statespace_R_df$model = "state_space"
statespace_R_df$date = incidence_df$date
poissonexp_R_df = poisson_confidence
colnames(poissonexp_R_df) = c("R_lb", "R_mode", "R_ub")
rownames(poissonexp_R_df) = NULL
poissonexp_R_df$model = "poisson_exp"
poissonexp_R_df$date = incidence_df$date
R_df = rbind(bayesian_R_df, statespace_R_df, poissonexp_R_df)

# Create base plots
incidences_plot = ggplot(data=incidence_df, aes(x=date, y=I)) +
  geom_line() +
  labs(x="", y = "Incidences")

max_y = 4
R_df$R_lb[R_df$R_lb > max_y] = NA
R_df$R_mode[R_df$R_mode > max_y] = NA
R_df$R_ub[R_df$R_ub > max_y] = NA
complete_R_plot = ggplot(data=R_df, aes(x=date, y=R_mode, fill=model, colour=model)) +
  geom_line(size=0.4) +
  geom_ribbon(aes(ymin=R_lb, ymax=R_ub), linetype=0, alpha=0.4) +
  labs(x="", y="Effective R\n") +
  scale_fill_manual(
    values = c(
      'bayesian' = 'blue',
      'state_space' = 'darkgreen',
      'poisson_exp' = 'deeppink1'
    ),
    labels = c(
      "bayesian" = "Bayesian",
      "state_space" = "State-Space",
      "poisson_exp" = "Exp. Poisson"
    )
  ) +
  scale_colour_manual(
    values = c(
      'bayesian' = 'blue',
      'state_space' = 'darkgreen',
      'poisson_exp' = 'deeppink1'
    ),
    labels = c(
      "bayesian" = "Bayesian",
      "state_space" = "State-Space",
      "poisson_exp" = "Exp. Poisson"
    )
  ) +
  guides(colour=FALSE) +
  coord_cartesian(ylim=c(0.6, max_y))
  
# Create and store zoom plots
R_zoom_1 = complete_R_plot + facet_zoom(xy = date > "2020-03-01" & date < "2020-04-15", horizontal = FALSE) + coord_fixed()
plot_grid(incidences_plot, R_zoom_1, ncol=1, axis="lr", align="v", rel_heights = c(1, 2))
ggsave("case_study/Rt_1_eng.png", width = 6, height=4)

R_zoom_2 = complete_R_plot + facet_zoom(xy = date > "2020-04-15" & date < "2020-06-15", horizontal = FALSE) + coord_fixed()
plot_grid(incidences_plot, R_zoom_2, ncol=1, axis="lr", align="v", rel_heights = c(1, 2))
ggsave("case_study/Rt_2_eng.png", width = 6, height=4)

R_zoom_3 = complete_R_plot + facet_zoom(xy = date > "2020-06-15" & date < "2020-08-15", horizontal = FALSE) + coord_fixed()
plot_grid(incidences_plot, R_zoom_3, ncol=1, axis="lr", align="v", rel_heights = c(1, 2))
ggsave("case_study/Rt_3_eng.png", width = 6, height=4)