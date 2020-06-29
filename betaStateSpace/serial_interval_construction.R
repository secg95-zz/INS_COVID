weightConstructionSerialInterval=function(lags=7, parameter1=2.235567, parameter2=5.419495){
  #Distribution: Weibull
  #Parameter 1: Shape Weibull
  #Parameter 2: Scale Weibull
  #computes Serial intervals?
  #Por que resto las diferencias?
  omegaLag=diff(pweibull(1:(lags+1),parameter1, parameter2 ))
  #ni idea que es esto?
  omegaLag=omegaLag/sum(omegaLag)
  # omegaLag=rev(omegaLag)
  return(omegaLag)
}
weightConstructionIncubacionInfeccion=function(lags=7, parameterINC1=3.169434, parameterINC2=5.163921, parameterINF1=24.206087, parameterINF2=2.984198 ){
  #Distribution Incubation: Weibull
  #Parameter INC 1: Shape Weibull
  #Parameter INC 2: Scale Weibull
  #Distribution Infection: Gamma
  #Parameter INF 1: Shape Gamma
  #Parameter INF 2: rate Gamma  

  prob=diff(pweibull(1:(lags+1),parameterINC1, parameterINC2 ))
  omegaLag=NULL
  for(Om in 1:lags){
    t=sum(prob[1:Om]*rev(pgamma(1:Om, parameterINF1, parameterINF2, lower.tail = FALSE)))
    omegaLag=c(omegaLag,t)
  }
  omegaLag=(parameterINF1/parameterINF2)*omegaLag/sum(omegaLag)
  omegaLag=rev(omegaLag)
  return(omegaLag)
}
weightConstructionIncubacion=function(lags=7, parameterINC1=3.169434, parameterINC2=5.163921 ){
  #Distribution Incubation: Weibull
  #Parameter INC 1: Shape Weibull
  #Parameter INC 2: Scale Weibull
  #Distribution Infection: Gamma
  #Parameter INF 1: Shape Gamma
  #Parameter INF 2: rate Gamma  
  
  omegaLag=diff(pweibull(1:(lags+1),parameterINC1, parameterINC2 ))
  omegaLag=omegaLag/sum(omegaLag)
  return(omegaLag)
}


betaThompson_non_constant_SI=function(base=base, window=7, initialCases=initialCases, lags=7,  update=1, parameter1=2.235567, parameter2=5.419495, omega){
  
  #Contructs de R with thompson method... (Cori et al)
  priorRate=3
  priorShape=1
  cases=base
  initialVector=initialCases
  
  #Generate expectation Matrix
  #calculo el I?
  longCases=c(initialVector,cases)


  omega_matrix=matrix(0,length(longCases)-lags,lags)
  for(i in 1:length(longCases)-lags){
    for(j in 1:lags){
    omega_matrix[i,j] = omega(i, j)
    }
  }
  cases=base
  initialVector=initialCases
  
  #Generate expectation Matrix
  #calculo el I?
  longCases=c(initialVector,cases)

  infectPotential_T=matrix(0,length(longCases)-lags,1)
  for(i in 1:dim(infectPotential_T)[1]){
    infectPotential_T[i]=sum(longCases[i:(i+lags-1)]*omega_matrix[i,])
  }
  #Seriel interval times observed cases, there are lags values
  # computing denominatortofequation 4?
  casesMatrix=matrix(0,length(infectPotential_T)-(window-1),window)
  for(i in 1:dim(casesMatrix)[1]){
    casesMatrix[i,]=infectPotential_T[i:(i+window-1)]
  }
  #doble suma ecuacion 7
  infectPotential=casesMatrix%*%matrix(1,window,1)
  
  #Generate observation Matrix
  postcasesMatrix=matrix(0,length(cases)-(window-1),window)
  for(i in 1:dim(postcasesMatrix)[1]){
    postcasesMatrix[i,]=cases[i:(i+window-1)]
  }
  
  #suma iw(t) se supone esto es Poisson
  infectObservados=postcasesMatrix%*%matrix(1,window,1)
  
  betaParameteres=cbind(priorShape, priorRate)
  
  for(i in 1:dim(casesMatrix)[1]){
    shape1=priorShape+infectObservados[i]
    rate1=priorRate+infectPotential[i]
    betaParameteres=rbind(betaParameteres,cbind(shape1, rate1))
  }
  shape1=betaParameteres[-1,1]
  rate1=betaParameteres[-1,2]
  
  lb=qgamma(0.05, shape=shape1, rate=rate1)
  ub=qgamma(0.95, shape=shape1, rate=rate1)
  beta=((shape1-1)/rate1)
  for(i in 1:window){
    lb=c(lb[1],lb)
    ub=c(ub[1],ub)
    beta=c(beta[1],beta)
  }
  
  
  data=cbind(lb, beta,ub)

  data=as.data.frame(data)
  data$day=1:length(lb)
  print(data$beta)
  data$R=data$beta*sum(omega_matrix[1,])
  data$Rlb=data$lb*sum(omega_matrix[1,])
  data$Rub=data$ub*sum(omega_matrix[1,]) 

  #R_Caso
  lb=data$lb
  ub=data$ub
  beta=data$beta
  longCaseslb=c(rep(lb[1], length=lags-1),lb)
  longCasesub=c(rep(ub[1], length=lags-1),ub)
  longCasesbeta=c(rep(beta[1], length=lags-1),beta)
  
  data$R_c=data$beta*sum(omega_matrix[1,])
  data$R_clb=data$lb*sum(omega_matrix[1,])
  data$R_cub=data$ub*sum(omega_matrix[1,])    
  
  for(i in 1:length(lb)){
    data$R_c[i]=sum(longCasesbeta[i:(i+lags-1)]*omega_matrix[1,])
    data$R_clb[i]=sum(longCaseslb[i:(i+lags-1)]*omega_matrix[1,])      
    data$R_cub[i]=sum(longCasesub[i:(i+lags-1)]*omega_matrix[1,])      
    
  }
  
  p<-ggplot(data=data, aes(x=day, y=R_c)) + 
    geom_line()+
    geom_ribbon(aes(ymin=R_clb, ymax=R_cub), linetype=2, alpha=0.1)+
    labs(title="Estimacion de la tasa de caso", x="Fecha de inicio de s?ntomas", y="R(t)")

  
  return(list(data=data, beta=data$beta, p=p))
}

betaThompson=function(base=base, window=7, initialCases=initialCases, lags=7,  update=1, parameter1=2.235567, parameter2=5.419495){
  omega=weightConstructionSerialInterval(lags, parameter1, parameter2)
  #Contructs de R with thompson method... (Cori et al)
  priorRate=3
  priorShape=1
  
  cases=base$newCases
  initialVector=initialCases
  
  #Generate expectation Matrix
  #calculo el I?
  longCases=c(initialVector,cases)

  casesMatrix_T=matrix(0,length(longCases)-lags,lags)
  print(casesMatrix_T)
  for(i in 1:dim(casesMatrix_T)[1]){
    casesMatrix_T[i,]=longCases[i:(i+lags-1)]
  }
  
  #Seriel interval times observed cases, there are lags values
  # computing denominatortofequation 4?
  infectPotential_T=casesMatrix_T%*%omega
  print(omega)
  print(sum(omega))
  casesMatrix=matrix(0,length(infectPotential_T)-(window-1),window)
  for(i in 1:dim(casesMatrix)[1]){
    casesMatrix[i,]=infectPotential_T[i:(i+window-1)]
  }
  #doble suma ecuacion 7
  infectPotential=casesMatrix%*%matrix(1,window,1)
  
  #Generate observation Matrix
  postcasesMatrix=matrix(0,length(cases)-(window-1),window)
  for(i in 1:dim(postcasesMatrix)[1]){
    postcasesMatrix[i,]=cases[i:(i+window-1)]
  }
  
  #suma iw(t) se supone esto es Poisson
  infectObservados=postcasesMatrix%*%matrix(1,window,1)
  
  betaParameteres=cbind(priorShape, priorRate)
  
  for(i in 1:dim(casesMatrix)[1]){
    shape1=priorShape+infectObservados[i]
    rate1=priorRate+infectPotential[i]
    betaParameteres=rbind(betaParameteres,cbind(shape1, rate1))
  }
  shape1=betaParameteres[-1,1]
  rate1=betaParameteres[-1,2]
  
  lb=qgamma(0.05, shape=shape1, rate=rate1)
  ub=qgamma(0.95, shape=shape1, rate=rate1)
  beta=((shape1-1)/rate1)
  print(beta)
  for(i in 1:window){
    lb=c(lb[1],lb)
    ub=c(ub[1],ub)
    beta=c(beta[1],beta)
  }
  
  
  data=cbind(lb, beta,ub)

  data=as.data.frame(data)
  print(1)
  data$day=base$Sintomas[1:length(lb)]
  print(2)
  data$R=data$beta*sum(omega)
  data$Rlb=data$lb*sum(omega)
  data$Rub=data$ub*sum(omega) 

  #R_Caso
  lb=data$lb
  ub=data$ub
  beta=data$beta
  longCaseslb=c(rep(lb[1], length=lags-1),lb)
  longCasesub=c(rep(ub[1], length=lags-1),ub)
  longCasesbeta=c(rep(beta[1], length=lags-1),beta)
  
  data$R_c=data$beta*sum(omega)
  data$R_clb=data$lb*sum(omega)
  data$R_cub=data$ub*sum(omega)    
  
  for(i in 1:length(lb)){
    data$R_c[i]=sum(longCasesbeta[i:(i+lags-1)]*omega)
    data$R_clb[i]=sum(longCaseslb[i:(i+lags-1)]*omega)      
    data$R_cub[i]=sum(longCasesub[i:(i+lags-1)]*omega)      
    
  }
  
  p<-ggplot(data=data, aes(x=day, y=R_c)) + 
    geom_line()+
    geom_ribbon(aes(ymin=R_clb, ymax=R_cub), linetype=2, alpha=0.1)+
    labs(title="Estimacion de la tasa de caso", x="Fecha de inicio de s?ntomas", y="R(t)")

  
  return(list(data=data, beta=data$beta, p=p))
}

#Incubacion-Infeccion
weightConstructionIncubacionInfeccion=function(lags=7, parameterINC1=3.169434, parameterINC2=5.163921, parameterINF1=24.206087, parameterINF2=2.984198 ){
  #Distribution Incubation: Weibull
  #Parameter INC 1: Shape Weibull
  #Parameter INC 2: Scale Weibull
  #Distribution Infection: Gamma
  #Parameter INF 1: Shape Gamma
  #Parameter INF 2: rate Gamma  

  prob=diff(pweibull(1:(lags+1),parameterINC1, parameterINC2 ))
  omegaLag=NULL
  for(Om in 1:lags){
    t=sum(prob[1:Om]*rev(pgamma(1:Om, parameterINF1, parameterINF2, lower.tail = FALSE)))
    omegaLag=c(omegaLag,t)
  }
  omegaLag=(parameterINF1/parameterINF2)*omegaLag/sum(omegaLag)
  omegaLag=rev(omegaLag)
  return(omegaLag)
}
weightConstructionIncubacion=function(lags=7, parameterINC1=3.169434, parameterINC2=5.163921 ){
  #Distribution Incubation: Weibull
  #Parameter INC 1: Shape Weibull
  #Parameter INC 2: Scale Weibull
  #Distribution Infection: Gamma
  #Parameter INF 1: Shape Gamma
  #Parameter INF 2: rate Gamma  
  
  omegaLag=diff(pweibull(1:(lags+1),parameterINC1, parameterINC2 ))
  omegaLag=omegaLag/sum(omegaLag)
  return(omegaLag)
}

betaStateSpace=function(base=base, initialCases=initialCases, lags=7, lags2=5, parameterINC1=3.169434, parameterINC2=5.163921, parameterINF1=24.206087, parameterINF2=2.984198){
  
  omega=weightConstructionIncubacionInfeccion(lags, parameterINC1, parameterINC2, parameterINF1, parameterINF2)
  f=weightConstructionIncubacion(lags2, parameterINC1, parameterINC2)
  
   priorMean=3/sum(omega)
  
  cases=base$newCases
  initialVector=initialCases
  
  #Dependiente
  casesMatrix_T=matrix(0,length(cases)-lags2+1,lags2)
  for(i in 1:dim(casesMatrix_T)[1]){
    casesMatrix_T[i,]=cases[i:(i+lags2-1)]
  }
  incubPotential_T=casesMatrix_T%*%f
  
  
  #Independiente
  longCases=c(initialVector,cases)
  casesMatrix_T=matrix(0,length(longCases)-lags,lags)
  for(i in 1:dim(casesMatrix_T)[1]){
    casesMatrix_T[i,]=longCases[i:(i+lags-1)]
  }
  infectPotential_T=casesMatrix_T%*%omega
  infectPotential_T=infectPotential_T[1:length(incubPotential_T)]
  
  pars=c(0,0)
    obj1 <- function(pars) {
      tryCatch({
      pars2=pars
      pars2=exp(pars2)
      x=log(infectPotential_T) 
      y=as.matrix(incubPotential_T)
      Z=t(cbind(as.matrix(x),1,0))
      n=length(incubPotential_T)
      dim(Z)=c(1,3,n)
      Z=as.array(Z)
      Ti=rep(c(1,0,0,0,1,0,0,1,1),n)
      dim(Ti)=c(3,3,n)
      
      R=rep(c(0,pars2[1],0,0,0,pars2[2]),n)
      dim(R)=c(3,2,n)
      
      a1=c(1,priorMean,priorMean)
      P1=diag(c(0,10^9,10^9))
      distribution="poisson"
      u=0
      model=ssm_mng(y, Z, Ti, R, a1, P1, distribution, state_names=c("Base", "Beta", "Beta_P"))
      
      objectiveFunction=-logLik(model, nsim = 1000, nsim_states = 0)
      if(objectiveFunction==-Inf){
        objectiveFunction=-objectiveFunction
      }
      }, error = function(e) {
        objectiveFunction=Inf
      } )
      return(objectiveFunction)
    }
    
    lb=c(-10,-10)
    ub=c(10,10)
    x0 = c(0,0)
    
    local_opts = list( "algorithm" = "NLOPT_LN_BOBYQA",
                       "xtol_rel" = 1.0e-7 )
    opts = list( "algorithm" = "NLOPT_LN_BOBYQA",
                 "xtol_rel" = 1.0e-7,
                 "maxeval" = 100000,
                 "local_opts" = local_opts )
    
    res = nloptr( x0=x0, eval_f=obj1, eval_grad_f=NULL,lb=lb, ub=ub, opts=opts)
    pars=res$solution
   
    
    pars2=pars
    pars2=exp(pars2)
    x=log(infectPotential_T) 
    y=as.matrix(incubPotential_T)
    Z=t(cbind(as.matrix(x),1,0))
    n=length(incubPotential_T)
    dim(Z)=c(1,3,n)
    Z=as.array(Z)
    Ti=rep(c(1,0,0,0,1,0,0,1,1),n)
    dim(Ti)=c(3,3,n)
    
    R=rep(c(0,pars2[1],0,0,0,pars2[2]),n)
    dim(R)=c(3,2,n)
    
    a1=c(1,priorMean,priorMean)
    P1=diag(c(0,10^9,10^9))
    distribution="poisson"
    u=0
    model=ssm_mng(y, Z, Ti, R, a1, P1, distribution, state_names=c("Base", "Beta", "Beta_P"))
    betaT=sim_smoother(model, nsim = 1)[-1,2,1] 
    for(h in 2:1000){
      betaT=cbind(betaT,sim_smoother(model, nsim = 1)[-1,2,1])  
    }
    
    beta=exp(betaT)

    simuQL=function(simu){
      return(quantile(simu,c(0.05,0.95)))
    }
    
    data=cbind(t(apply(beta, 1, simuQL)),apply(beta, 1, mean))
    colnames(data)=c("lb","ub", "beta")
    
    # lb=qgamma(0.05, shape=shape1, rate=rate1)
    # ub=qgamma(0.95, shape=shape1, rate=rate1)
    # beta=((shape1-1)/rate1)
    # 
    # data=cbind(betaParameteres,lb, beta,ub)
    data=as.data.frame(data)
    data$day=base$Sintomas[1:dim(beta)[1]]
    
    data$R=data$beta*sum(omega)
    data$Rlb=data$lb*sum(omega)
    data$Rub=data$ub*sum(omega)
    
    
    #R_Caso
    lb=data$lb
    ub=data$ub
    beta=data$beta
    longCaseslb=c(rep(lb[1], length=lags-1),lb)
    longCasesub=c(rep(ub[1], length=lags-1),ub)
    longCasesbeta=c(rep(beta[1], length=lags-1),beta)
    
    data$R_c=data$beta*sum(omega)
    data$R_clb=data$lb*sum(omega)
    data$R_cub=data$ub*sum(omega)    
    
    for(i in 1:dim(betaT)[1]){
      data$R_c[i]=sum(longCasesbeta[i:(i+lags-1)]*omega)
      data$R_clb[i]=sum(longCaseslb[i:(i+lags-1)]*omega)      
      data$R_cub[i]=sum(longCasesub[i:(i+lags-1)]*omega)      
      
    }

    p<-ggplot(data=data, aes(x=day, y=R_c)) + 
      geom_line()+
      geom_ribbon(aes(ymin=R_clb, ymax=R_cub), linetype=2, alpha=0.1)+
      labs(title="Estimacion de la tasa de caso", x="Fecha de inicio de s?ntomas", y="R(t)")
    print(p)
    
     return(list(data=data, p=p))
} 