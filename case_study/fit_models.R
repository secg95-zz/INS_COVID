library(RJSONIO)

# Construct I series
# Import incidence time series and assert that it doesn't skip any steps
colombia_data = read.csv("./case_study/I.csv", sep=",")
colombia_data$Fecha.Inicio.Sintomas = as.Date(colombia_data$Fecha.Inicio.Sintomas)
date_range = seq(
  min(colombia_data$Fecha.Inicio.Sintomas),
  max(colombia_data$Fecha.Inicio.Sintomas),
  by="days"
)
date_range = data.frame(date_range)
colombia_data = merge(
  colombia_data, date_range,
  by.x = "Fecha.Inicio.Sintomas", by.y = "date_range",
  all = TRUE,
  sort = TRUE
)
# Subset the data to the period considered in the paper
I = colombia_data$cases[
  (colombia_data$Fecha.Inicio.Sintomas >= as.Date("2020-03-06"))
  & (colombia_data$Fecha.Inicio.Sintomas <= as.Date("2020-08-15"))
]
date_range = date_range$date_range[
  (date_range$date_range >= as.Date("2020-03-06"))
  & (date_range$date_range <= as.Date("2020-08-15"))
]

# Fit each model
SI_SHAPE = 2.23
SI_SCALE = 5.42
F_INC_SHAPE = 3.16
F_INC_RATE = 5.16
F_INF_SCALE = 24.2
F_INF_SHAPE = 2.98
R_PRIOR_RATE = 2.28 / (0.11735^2)
R_PRIOR_MEAN = 1.2
R_PRIOR_SHAPE = (2.68 * R_PRIOR_RATE) / R_PRIOR_MEAN
BAYESIAN_WINDOW = 5
MEAN_INFECTIOUS_TIME = F_INF_SCALE * base::gamma(1 + 1/(F_INF_SHAPE))
# Fit Bayesian model
bayesian = new.env(); source("bayesian/model_fitting.r", local=bayesian)
bayesian_fit = bayesian$fit2(
  I,
  window = 5,
  prior_shape = R_PRIOR_SHAPE, prior_rate = R_PRIOR_RATE,
  omega = diff(pweibull(0:length(I), shape=SI_SHAPE, scale=SI_SCALE))
)
# Fit State-Space model
ss = new.env(); source("state_space/model_fitting.r", local=ss)
ss_fit = ss$fit(
  I,
  f_inc = c(1e-30, diff(pgamma(0:(length(I) - 1), shape=F_INC_SHAPE, rate=F_INC_RATE))),
  f_inf = c(1e-30, diff(pweibull(0:(length(I) - 1), shape=F_INF_SHAPE, scale=F_INF_SCALE))),
  prior_R = R_PRIOR_MEAN,
  lags=9
)$data
# Fit Exponential Poisson model
poisson = new.env(); source("likelihood_incubation_exp/model_fitting.r", local=poisson)
N_ROBUST_ITER = 10
N_NLOPTR_ITER = 250000
# Robust algorithm with 1 iteration, and nloptr with one iteration, takes 0.07s
LAMBDA = 2^12
poisson_fit = poisson$fit_robust(
  observed_I = I,
  beta_min = 0.7 / 8.111421, beta_max = 3 / 8.111421,
  tau1_min = 5, tau1_max = 5,# tau1_min = 0.1176652, tau1_max = 0.1176652,# tau1_min = 1 / (4.5 + 2),  tau1_max = 1 / (4.5 - 2), # tau1_min = 1 / (1.915495 + 1),  tau1_max = 1 / (1.915495 - 1),
  tau2_min = 1/8.111421, tau2_max = 1/8.111421, # tau2_min = 1/ (8.111421 + 8),  tau2_max = 1/ (8.111421 - 2),
  N0_min = 10, N0_max = 60,# N0_min = 10, N0_max = 10,# N0_min = 1, N0_max = 50,
  A0_min = 10 / pexp(1, rate=1/8.111421), A0_max = 10 / pexp(1, rate=1/8.111421),# A0_min=1, A0_max=50,
  lambda=LAMBDA,
  ignore_beta_diff=NULL,
  use_history=FALSE,
  n_robust_iter=N_ROBUST_ITER,
  n_nloptr_iter=N_NLOPTR_ITER
)

# Format and store results
fitted_models = list(
  I = I,
  date_range = as.character(date_range),
  bayesian = list(
    R_mode = bayesian_fit$mode,
    R_ub = bayesian_fit$ub,
    R_lb = bayesian_fit$lb
  ),
  state_space = list(
    R_mode = ss_fit$R,
    R_ub = ss_fit$Rub,
    R_lb = ss_fit$Rlb
  ),
  poisson_exp = list(
    R_mode = poisson_fit$R,
    R_ub = NULL,
    R_lb = NULL
  )
)
write(toJSON(fitted_models), "./case_study/fitted_models.json")