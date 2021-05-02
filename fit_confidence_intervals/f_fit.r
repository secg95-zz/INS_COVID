#Castigo = 25

fit = function(observed_I, beta0, beta_min, beta_max, tau10, tau1_min, tau1_max,
               tau20, tau2_min, tau2_max,
               N00, N0_min, N0_max, Castigo, ignore_beta_diff, max_eval = 1000 ) {
  "
  Fits the model to observed daily new case counts.
  
  Parameters
  ----------
  observed_I : numeric vector
    Daily observed number of new cases. (son los incubadores!! o son los infecciosos nuevos)
  beta0 : numeric vector
    Initial value for beta; the expected number of cases (incubadores) stemming from a single
    person in a single day.
  beta_min : numeric
    Beta lower bound.
  beta_max : numeric
    Beta upper bound.
  ignore_beta_diff : numeric vector
    List of beta indices for which differences should be ignored while
    calculating the loss function. This amounts to moments in time in which we
    allow the beta series to be discontinuous.
  
  Returns
  -------
  beta : numeric vector
    Expected number of cases stemming from a single person in a single day.
  "
  # concatenate initial parameters into a single vector for optimization
  x0 = c(beta0, tau10, tau20, N00)
  steps = length(beta0)
  lb = c(rep(beta_min, steps), tau1_min, tau2_min, N0_min)
  ub = c(rep(beta_max, steps), tau1_max, tau2_max,  N0_max)
  # set optimization parameters
  opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-7, "maxeval" = max_eval)
  print(opts)
  # define loss function
  loss = function(x) {
    # unpack values
    beta = x[1:steps]
    tau1 = x[steps + 1]
    tau2 = x[steps + 2]
    N0 = x[steps + 3]
    # calculate loss
    expected = get_expected(beta, N0, tau1, tau2) # Obtiene el número de casos infectados
    beta_diff = diff(beta)
    regularization1 =  (beta_diff/beta[1:length(beta_diff)]) 
    regularization =  (diff(regularization1)) ^ 2
    regularization = regularization[!1:length(regularization) %in% ignore_beta_diff]
    observed_I=round(observed_I,0)
    factSum=function(x){
      return(sum(log(1:round(x,0))))
    }
    loglikelihood=observed_I*log(expected$NN)-expected$NN-sapply(as.matrix(observed_I), FUN=factSum)
    loglikelihood[observed_I==0]=-expected$NN[observed_I==0]
    
    loss = -mean(loglikelihood) + Castigo * mean(regularization)
    #loss=sum((observed_I-expected$I)^2)+ sum(regularization)
    loss
  }
  result = nloptr(x0=x0, eval_f=loss, eval_grad_f=NULL, lb=lb, ub=ub, opts=opts)
  model = list()
  model$beta = result$solution[1:steps]
  model$tau1 = result$solution[steps + 1]
  model$tau2 = result$solution[steps + 2]
  model$N0 = result$solution[steps + 3]
  model$loss = result$objective
  beta=model$beta
  beta_diff = diff(beta)
  regularization1 =  (beta_diff/beta[1:length(beta_diff)]) 
  regularization =  (diff(regularization1)) ^ 2
  regularization = regularization[!1:length(regularization) %in% ignore_beta_diff]
  model$likelihood=(model$loss-Castigo*mean(regularization))*length(beta)
  return(model)
}
