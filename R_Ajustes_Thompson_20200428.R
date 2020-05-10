#Loading package
library(extraDistr)
library(nloptr)
library(ggplot2)
library(deSolve)
library(dplyr)
library(matrixStats)
library(readxl)
library(R0)
library(gridExtra)
library(ggpubr)
library(moments)
#Funciones
simulacionProcesoA=function(BetaList=rep(0.4,times=50), casosIniciales=1, InfeccionP=5.6){
  moments=length(BetaList)
  nCasosT=matrix(0,1,moments)
  casosT=matrix(0,1,moments)
  ponderador_t=BetaList+exp(-1/InfeccionP)
  ponderador_n=BetaList
  temp=casosIniciales
  for(i in 1:moments){
    casosT[i]=temp*ponderador_t[i]
    nCasosT[i]=temp*ponderador_n[i]
    temp=casosT[i]
  }
  return(nCasosT)
}
simulacionProcesoA_Infeccion=function(base=ColombiaB, BetaList=rep(0.4,times=50), casosIniciales=1, InfeccionP=5.6, InfeccionPMin=5, InfeccionPMax=7){
  ColombiaB=base
  x=InfeccionP
  eval_f <- function( x ) {
    reales=ColombiaB$newCases
    simulados=simulacionProcesoA(BetaList, casosIniciales, x)
    simulados=simulados[1:length(reales)]
    asintomaticos=sum(simulados*reales)/sum(simulados*simulados)
    asintomaticos=1
    nResultado=mean((reales-asintomaticos*simulados)^4)
    return( nResultado )
  }
  # initial values
  x0 <- x
  # lower and upper bounds of control
  lb =InfeccionPMin #Infeccion
  ub =InfeccionPMax #Infeccion
  
  #Optimization Methods>
  #NLOPT_LN_BOBYQA
  #NLOPT_LN_COBYLA
  
  local_opts <- list( "algorithm" = "NLOPT_LN_BOBYQA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LN_BOBYQA",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 100000,
                "local_opts" = local_opts )
  #Correr las siguientes 3 linea si se cuenta con una hora
  res <- nloptr( x0=x0, eval_f=eval_f, eval_grad_f=NULL,lb=lb, ub=ub, opts=opts)
  return(res$solution)
}
BetaUpdateOLSA=function(base=ColombiaB,BetaList=BetaList,lambda=1600000){
  ColombiaB=base
  Colombia=base
  casosTotales=length(ColombiaB$newCases)
  #for(i in 2:casosTotales){
  #Colombia=ColombiaB[1:i,]
  Colombia=ColombiaB
  moments=length(ColombiaB$newCases)
  x=BetaList
  
  
  eval_f <- function( x ) {
    
    BetaList=x
    
    casosIniciales=1
    InfeccionP=simulacionProcesoA_Infeccion(ColombiaB, BetaList, casosIniciales, InfeccionP=5.6, InfeccionPMin=5, InfeccionPMax=7)
    reales=ColombiaB$newCases
    simulados=simulacionProcesoA(BetaList, casosIniciales, InfeccionP)
    simulados=simulados[1:length(reales)]
    asintomaticos=sum(simulados*reales)/sum(simulados*simulados)
    asintomaticos=1
    resultado=mean((reales-asintomaticos*simulados)^4)+lambda*mean((BetaList[-length(BetaList)]-BetaList[-1])^2)
    
    for(i in 2:5){
      casosIniciales=i
      InfeccionP=simulacionProcesoA_Infeccion(ColombiaB, BetaList, casosIniciales, InfeccionP=5.6, InfeccionPMin=5, InfeccionPMax=7)
      reales=ColombiaB$newCases
      simulados=simulacionProcesoA(BetaList, casosIniciales, InfeccionP)
      simulados=simulados[1:length(reales)]
      asintomaticos=sum(simulados*reales)/sum(simulados*simulados)
      asintomaticos=1
      nResultado=mean((reales-asintomaticos*simulados)^4)+lambda*mean((BetaList[-length(BetaList)]-BetaList[-1])^2)
      if(nResultado<resultado){
        resultado=nResultado
      }
    }
    return( resultado )
  }
  # initial values
  x0 <- x
  # lower and upper bounds of control
  lb =rep(0.1,times=moments) #Beta
  ub =rep(0.7,times=moments) #Beta
  
  #Optimization Methods>
  #NLOPT_LN_BOBYQA
  #NLOPT_LN_COBYLA
  
  start.time <- Sys.time()
  local_opts <- list( "algorithm" = "NLOPT_LN_BOBYQA",
                      "xtol_rel" = 1.0e-7 )
  opts <- list( "algorithm" = "NLOPT_LN_BOBYQA",
                "xtol_rel" = 1.0e-7,
                "maxeval" = 100000,
                "local_opts" = local_opts )
  #Correr las siguientes 3 linea si se cuenta con una hora
  res <- nloptr( x0=x0, eval_f=eval_f, eval_grad_f=NULL,lb=lb, ub=ub, opts=opts)
  cbind(lb,res$solution,ub)
  BetaList=res$solution
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  return(BetaList)
}
BetaAR0A=function(base=ColombiaB,BetaList=BetaList,lambda=1600000){
  ColombiaB=base
  casosIniciales=1
  InfeccionP=simulacionProcesoA_Infeccion(ColombiaB, BetaList, casosIniciales, InfeccionP=5.6, InfeccionPMin=5, InfeccionPMax=7)
  infP_final=InfeccionP
  reales=ColombiaB$newCases
  simulados=simulacionProcesoA(BetaList, casosIniciales, InfeccionP)
  simulados=simulados[1:length(reales)]
  asintomaticos=sum(simulados*reales)/sum(simulados*simulados)
  asintomaticos=1
  resultado=mean((reales-asintomaticos*simulados)^4)+lambda*mean((BetaList[-length(BetaList)]-BetaList[-1])^2)
  ni=1
  serieE=reales-asintomaticos*simulados
  for(i in 2:5){
    casosIniciales=i
    InfeccionP=simulacionProcesoA_Infeccion(ColombiaB, BetaList, casosIniciales, InfeccionP=5.6, InfeccionPMin=5, InfeccionPMax=7)
    reales=ColombiaB$newCases
    simulados=simulacionProcesoA(BetaList, casosIniciales, InfeccionP)
    simulados=simulados[1:length(reales)]
    asintomaticos=sum(simulados*reales)/sum(simulados*simulados)
    asintomaticos=1
    nResultado=mean((reales-asintomaticos*simulados)^4)+lambda*mean((BetaList[-length(BetaList)]-BetaList[-1])^2)
    if(nResultado<resultado){
      ni=i
      resultado=nResultado
      serieE=reales-asintomaticos*simulados
      infP_final=InfeccionP
    }
  }
  R0=BetaList*infP_final
  
  return(list(R0=R0, ni=ni, InfeccionP=infP_final, escala=asintomaticos, serieE=serieE))
}
weightConstruction=function(lags=7, distribution="uniform", parameter1=5, parameter2=7){
  #Distribution: uniform or gamma
  #Parameter 1: uniform lower bound or gamma shape
  #Parameter 2: uniform upper bound or gamma rate
  if(distribution=="uniform"){
    omegaLag=1-punif(1:lags,min=parameter1, max=parameter2)
    if(length(omegaLag)<lags){
      omegaLag=cbind(omegaLag, rep(0,length=lags-length(omegaLag)))
    }else{
      omegaLag=omegaLag[1:lags]
    }
  }else{
    omegaLag=1-pgamma(1:lags,shape=parameter1, rate = parameter2)
    
  }
  #Change order
  omegaLag=omegaLag[lags+1-(1:lags)]
  omegaLag/sum(omegaLag)
  return(omegaLag)
}
R0Thompson=function(base=ColombiaB,window=2, lags=7, initialCases=1, update=1, distribution="gamma", parameter1=5.8, parameter2=1){
  omega=weightConstruction(lags=lags, distribution, parameter1, parameter2)
  priorRate=2.28/0.11735^2
  priorShape=2.68*priorRate/sum(omega)
  
  cases=base$newCases
  initialVector=rep(0, length=lags)
  initialVector[lags]=initialCases
  
  #Generate expectation Matrix
  longCases=c(initialVector,cases)
  casesMatrix_T=matrix(0,length(longCases)-lags,lags)
  for(i in 1:dim(casesMatrix_T)[1]){
    casesMatrix_T[i,]=longCases[i:(i+lags-1)]
  }
  infectPotential_T=casesMatrix_T%*%omega
  casesMatrix=matrix(0,length(infectPotential_T)-(window-1),window)
  for(i in 1:dim(casesMatrix)[1]){
    casesMatrix[i,]=infectPotential_T[i:(i+window-1)]
  }
  infectPotential=casesMatrix%*%matrix(1,window,1)
  
  #Generate observation Matrix
  postcasesMatrix=matrix(0,length(cases)-(window-1),window)
  for(i in 1:dim(postcasesMatrix)[1]){
    postcasesMatrix[i,]=cases[i:(i+window-1)]
  }
  
  infectObservados=postcasesMatrix%*%matrix(1,window,1)
  
  R0Parameteres=cbind(priorShape, priorRate)
  
  for(i in 1:dim(casesMatrix)[1]){
    shape1=priorShape+infectObservados[i]
    rate1=priorRate+infectPotential[i]
    R0Parameteres=rbind(R0Parameteres,cbind(shape1, rate1))
    if(update==1){
      priorShape=shape1
      priorRate=rate1
    }
  }
  shape1=R0Parameteres[,1]
  rate1=R0Parameteres[,2]
  
  lb=qgamma(0.05, shape=shape1, rate=rate1)
  ub=qgamma(0.95, shape=shape1, rate=rate1)
  R0=((shape1-1)/rate1)
  
  R0Matrix_T=matrix(0,length(R0)-lags,lags)
  for(i in 1:dim(R0Matrix_T)[1]){
    R0Matrix_T[i,]=R0[i:(i+lags-1)]
  }
  R0_Adj=R0Matrix_T%*%omega
  
  lbMatrix_T=matrix(0,length(lb)-lags,lags)
  for(i in 1:dim(lbMatrix_T)[1]){
    lbMatrix_T[i,]=lb[i:(i+lags-1)]
  }
  lb_Adj=lbMatrix_T%*%omega
  
  ubMatrix_T=matrix(0,length(ub)-lags,lags)
  for(i in 1:dim(ubMatrix_T)[1]){
    ubMatrix_T[i,]=ub[i:(i+lags-1)]
  }
  ub_Adj=ubMatrix_T%*%omega
  
  
  
  
  data=cbind(R0Parameteres,lb,ub)
  data=as.data.frame(data)
  data$day=ColombiaB$Sintomas[1:length(lb)]
  data=data[(lags+1):length(lb),]
  data$lb=lb_Adj
  data$mode=R0_Adj
  data$ub=ub_Adj
  
  
  p<-ggplot(data=data, aes(x=day, y=mode)) + 
    geom_point() + geom_line()+geom_line(aes(x=day, y=1))+
    geom_ribbon(aes(ymin=lb, ymax=ub), linetype=2, alpha=0.1)+
    ylim(0,3)+labs(title="Estimación de la tasa efectiva de reproducción", x="Fecha de inicio de síntomas", y="R(t)")
  print(p)
  
  return(list(data=data, R0=data$mode))
}


#Cargar Bases

setwd("D:/Colombia")
load("CasosAjustados_Diagnostico.Rda")
data$Sintomas=as.Date(data$Sintomas, format="%Y/%m/%d")
COLstartDay <- min(data$Sintomas)
day <- COLstartDay
data$day <- as.numeric(difftime(data$Sintomas,day, units="days"))
tempCol <- data.frame(day = 0:max(data$day))
ColombiaB <- merge(tempCol, data, all.x = TRUE)
for(i in 1:length(ColombiaB$day)){
  if(is.na(ColombiaB$Sintomas[i])){
    ColombiaB$Sintomas[i]=ColombiaB$Sintomas[i+1]-1
  }
}
ColombiaB[is.na(ColombiaB)] <- 0
ColombiaB$newCases=ColombiaB$`Estimados Totales`
ColombiaB$cumCases <- cumsum(ColombiaB$newCases)


##Identificacion del R0 via Poisson Exponencial

result=NULL
resultC=NULL

BetaList = read.table("D:/Colombia/BetaInitials.txt", quote="\"", comment.char="")
BetaList = as.numeric(BetaList$V1)
if(length(BetaList)<length(ColombiaB$newCases)){
  while(length(BetaList)<length(ColombiaB$newCases)){
    BetaList=c(BetaList,0.4)
  }
}else{
  BetaList=BetaList[1:length(ColombiaB$newCases)]
}
#Grafica de residuos
lambda=6.4*10^9
BetaList=BetaUpdateOLSA(base=ColombiaB,BetaList=BetaList, lambda)
Ret=BetaAR0A(base=ColombiaB,BetaList, lambda)
data=as.data.frame(Ret$R0)
colnames(data)[1]=c(paste("RE(lambda=", lambda,")", sep=""))
daysR=ColombiaB$Sintomas[1:length(Ret$R0)]
data$Dias=daysR

lags=7
window=5
R0_Th=R0Thompson(base=ColombiaB,window=window, lags=lags, initialCases=5, update=0, distribution="uniform", parameter1=5, parameter2=7)
data$R0T=NA
data=data[data$Dias%in%R0_Th$data$day,]
data$R0T[data$Dias%in%R0_Th$data$day]=R0_Th$R0
colnames(data)[3]=paste("Th(w=",window,", lag=",lags,")",sep="")
  
  
ajusteR <- reshape2::melt(data, id.var="Dias")
colnames(ajusteR)[2]="Escenario"
p <- ggplot(ajusteR, aes(x=Dias, y=value, col=Escenario)) + 
  geom_line()+geom_line(aes(color=Escenario))+
  geom_point(aes(color=Escenario))+
  ggtitle("Periodo de Ajuste") +ylim(0,3)+
  labs(x="Dias", y="Tasa Efectiva de Reproduccion")+ scale_y_continuous(labels = scales::comma)
p





