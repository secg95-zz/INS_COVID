#Parametros para Incubacion-Infeccion
parameterINC1=3.169434
parameterINC2=5.163921
parameterINF1=24.206087
parameterINF2=2.984198
lags=21
lags2=5
window=7


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
    
     return(list(data=data))
} 


betaThompson(base=base, window=7, initialCases=initialCases, lags=21,  update=1, parameter1=2.235567, parameter2=5.419495)
