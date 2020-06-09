library(nloptr)
library(rlist)
library(Rlab)
##Funciones
gen_simulated=function(beta, N0, lambda){
  casos_Serie=rep(0,length=length(beta)+40)
  casos_Serie[1:gen_casosDisc(lambda)]=N0
  casos_Nuevo=N0
  casos_SerieV=N0
  for(i in 1:length(beta)){
    nCasos=rpois(1,beta[i]*casos_Serie[i])
    casos_SerieV[i]=casos_Serie[i+1]
    for(j in 1: nCasos){
      diasInf=gen_casosDisc(lambda)
      casos_Serie[(i+1):min(i+1+diasInf,length(casos_Serie))]=casos_Serie[(i+1):min(i+1+diasInf,length(casos_Serie))]+1
    }
    casos_Nuevo=c(casos_Nuevo,nCasos)
  }
  return(casos_Nuevo)
}
gen_casosDisc=function(lambda){
  casos=0
  r=1
  while(rbern(1,pexp(1, rate=lambda, lower.tail=FALSE))==1){
    casos=casos+1
  }
  return(casos)
}  

get_expected = function(beta, N0, lambda) {
  "
  Returns
  -------
  expected_I : numeric vector
    Daily expected number of new cases.
  expected_N : numeric vector
    Daily expected number of infectious cases.
  "
  moments = length(beta)
  expected_I = rep(0, moments)
  expected_N = rep(0, moments)
  N_previous = N0
  for(t in 1:moments){
    expected_N[t] = N_previous * (beta[t] + pexp(1, rate=lambda, lower.tail=FALSE))
    expected_I[t] = N_previous * beta[t]
    N_previous = expected_N[t]
  }
  expected_N = c(N0, expected_N)
  expected = list("I"=expected_I, "N"=expected_N)
  return(expected)
}

fit = function(observed_I, beta0, beta_min, beta_max, lambda0, lambda_min,
               lambda_max, N00, N0_min, N0_max, Castigo, ignore_beta_diff) {
  "
  Fits the model to observed daily new case counts.
  
  Parameters
  ----------
  observed_I : numeric vector
    Daily observed number of new cases.
  beta0 : numeric vector
    Initial value for beta; the expected number of cases stemming from a single
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
  x0 = c(beta0, lambda0, N00)
  steps = length(observed_I)
  lb = c(rep(beta_min, steps), lambda_min, N0_min)
  ub = c(rep(beta_max, steps), lambda_max, N0_max)
  # set optimization parameters
  opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-7, "maxeval" = 10000000)
  # define loss function
  loss = function(x) {
    # unpack values
    beta = x[1:steps]
    lambda = x[steps + 1]
    N0 = x[steps + 2]
    # calculate loss
    expected = get_expected(beta, N0, lambda)
    beta_diff = diff(beta)
    regularization =  (beta_diff/beta[1:length(beta_diff)]) ^ 2
    regularization = regularization[!1:length(regularization) %in% ignore_beta_diff]
    observed_I=round(observed_I,0)
    factSum=function(x){
      return(sum(log(1:round(x,0))))
    }
    loglikelihood=observed_I*log(expected$I)-expected$I-sapply(as.matrix(observed_I), FUN=factSum)
    loglikelihood[observed_I==0]=-expected$I[observed_I==0]
    loss = -mean(loglikelihood)+ Castigo*mean(regularization)
    #loss=sum((observed_I-expected$I)^2)+ sum(regularization)
    
  }
  result = nloptr(x0=x0, eval_f=loss, eval_grad_f=NULL, lb=lb, ub=ub, opts=opts)
  model = list()
  model$beta = result$solution[1:steps]
  model$lambda = result$solution[steps + 1]
  model$N0 = result$solution[steps + 2]
  model$loss = result$objective
  beta=model$beta
  beta_diff = diff(beta)
  regularization =  (beta_diff/beta[1:length(beta_diff)]) ^ 2
  regularization = regularization[!1:length(regularization) %in% ignore_beta_diff]
  model$likelihood=(model$loss-Castigo*mean(regularization))*length(beta)
  return(model)
}


##Ejercicios


beta=c(seq(0.4,0.1,length=10), seq(0.6,0.3,length=10), seq(0.1,0.3,length=10), seq(0.5,0.1,length=10), seq(0.5,0.1,length=10))
N0=1
lambda=1/7


observad0s=gen_simulated(beta, N0, lambda)
  for(h in 1:2){
    observad0s=(observad0s*h+gen_simulated(beta, N0, lambda))/(h+1)
  }

get_expected(beta, N0, lambda)
observed_I=observad0s[-1]
beta0=beta
beta_min=0.01
beta_max=0.8
lambda0=lambda
lambda_min=lambda
lambda_max=lambda
N00=1
N0_min=1
N0_max=1
alpha=1
ignore_beta_diff=100

LossFunction=NULL
castigos=0:20
for(C in castigos){
sol=fit(observed_I, beta0, beta_min, beta_max, lambda0, lambda_min,
               lambda_max, N00, N0_min, N0_max, Castigo=C, ignore_beta_diff)
LossFunction=c(LossFunction,sol$likelihood)
beta0=sol$beta
}

#Calcular Likelihood Ratio Test

LRTest=2*(LossFunction-LossFunction[1])

prob=pchisq(LRTest,length(beta), lower.tail=TRUE)
plot(castigos,prob)
lines(castigos,castigos*0+0.95, type="l", col="red")
#Elegir Parametro de suavizamiento
C_opt=min(castigos[prob>0.95])


sol=fit(observed_I, beta0, beta_min, beta_max, lambda0, lambda_min,
        lambda_max, N00, N0_min, N0_max, Castigo=C, ignore_beta_diff)

serie=1:length(beta)
#Beta estimado vs. Beta teorico
plot(serie, beta, type="l")
lines(serie, sol$beta, type="l", col="red")
#Ajuste
plot(serie, observed_I, type="l")
lines(serie, get_expected(sol$beta, N0, lambda)$I, type="l", col="red")
