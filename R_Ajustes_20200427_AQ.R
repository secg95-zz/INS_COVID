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


#Cargar Bases

# setwd("D:/Colombia")
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


##Identificacion del R0

result=NULL
resultC=NULL

moments=length(ColombiaB$newCases)
exponentials=seq(2,12, length=2)
lambdaValues2=64*10^exponentials
BetaList = read.table("BetaInitials.txt", quote="\"", comment.char="")
BetaList = as.numeric(BetaList$V1)
if(length(BetaList)<length(ColombiaB$newCases)){
  while(length(BetaList)<length(ColombiaB$newCases)){
    BetaList=c(BetaList,0.4)
  }
}else{
  BetaList=BetaList[1:length(ColombiaB$newCases)]
}
for(lambda in lambdaValues2){
  BetaList=BetaUpdateOLSA(base=ColombiaB,BetaList=BetaList, lambda)
  Ret=BetaAR0A(base=ColombiaB,BetaList, lambda)
  R0=Ret$R0
  casos=Ret$escala*simulacionProcesoA(BetaList, casosIniciales =Ret$ni,InfeccionP = Ret$InfeccionP )
  casos=t(casos)
  result=cbind( result,R0)
  resultC=cbind( resultC,casos)
  print(Ret$InfeccionP)
  print(lambda)
  print(Ret$ni)
}
write.table(BetaList, "BetaInitials.txt", row.names = FALSE, col.names = FALSE)

result=cbind(result)
resultC=cbind(resultC, ColombiaB$newCases[1:length(resultC[,1])])

colnames(result)=c(paste("PE(Lambda=",lambdaValues2,")",sep=""))

colnames(resultC)=c(paste("PE(Lambda=",lambdaValues2,")",sep=""), "Casos Observados")

daysR=ColombiaB$Sintomas[1:length(resultC[,1])]
data=as.data.frame(result)
data$Dias=daysR
data_R0=data
ajusteR <- reshape2::melt(data, id.var="Dias")
colnames(ajusteR)[2]="Escenario"
p <- ggplot(ajusteR, aes(x=Dias, y=value, col=Escenario)) + 
  geom_line()+geom_line(aes(color=Escenario))+
  geom_point(aes(color=Escenario))+
  ggtitle("Periodo de Ajuste") +
  labs(x="Dias", y="Tasa Efectiva de Reproduccion")+ scale_y_continuous(labels = scales::comma)
p

daysC=ColombiaB$Sintomas[1:length(resultC[,1])]
data=as.data.frame(resultC)
data$Dias=daysC
data_Ajuste=data
ajusteC <- reshape2::melt(data, id.var="Dias")
colnames(ajusteC)[2]="Escenario"
p <- ggplot(ajusteC, aes(x=Dias, y=value, col=Escenario)) + 
  geom_line()+geom_line(aes(color=Escenario))+
  geom_point(aes(color=Escenario))+
  ggtitle("Periodo de Ajuste") +
  labs(x="Dias", y="Infectados")+ scale_y_continuous(labels = scales::comma)
p
save(data_R0,file="D:/Colombia/dataR0_C.Rda")
save(data_Ajuste,file="D:/Colombia/dataAjuste_C.Rda")


#Grafica de residuos
lambda=6.4*10^10
# BetaList=BetaUpdateOLSA(base=ColombiaB,BetaList=BetaList, lambda)
BetaList=BetaUpdateOLSA(base=ColombiaB,BetaList=rep(0.3, 49), lambda)
Ret=BetaAR0A(base=ColombiaB,BetaList, lambda)
errores=as.data.frame(Ret$serieE)
colnames(errores)[1]=c("Residuos")
daysR=ColombiaB$Sintomas[1:length(resultC[,1])]
errores$Dias=daysR
p1 <- ggplot(errores, aes(x=Dias, y=Residuos)) + 
  geom_point()+
  ggtitle("Dispersion temporal") +
  labs(x="Dias", y="Residuos")+ scale_y_continuous(labels = scales::comma)
p2=ggplot(errores, aes(x=Residuos)) +
  geom_histogram(aes(y=..density..), fill="lightblue", color="blue", position="identity", alpha=0.5, bins = 20)+
  geom_vline(aes(xintercept=mean(Residuos, na.rm = TRUE)), color="blue", linetype="dashed")+
  geom_density(alpha=0.6, color="black")+
  labs(title="Histograma de residuos", x="Residuo", y="Densidad")

fig=ggarrange(p1, p2, ncol = 1, nrow = 2)
annotate_figure(fig,
                top = text_grob("Distribucion de residuos", color = "black", face = "bold", size = 14, hjust = 0.5),
)






mean(Ret$serieE)
##################################################################################################

