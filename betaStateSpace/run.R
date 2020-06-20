library(readxl)
library(ggplot2)
library(bssm)
library(nloptr)


source("betaStateSpace/serial_interval_construction.R")

#Parametros para Incubacion-Infeccion
params = list(
"parameterINC1"=3.169434,
"parameterINC2"=5.163921,
"parameterINF1"=24.206087,
"parameterINF2"=2.984198,
"lags"=21,
"lags2"=5,
"window"=7,
"parameter1"=2.235567, 
"parameter2"=5.419495
)

#read file
baseCasos <- read_excel("IRA_Departamentos_2013_2019.xlsx", sheet = "Hoja2")

#Base
municipioReferencia="BOGOTA"
initialCases=as.matrix(baseCasos[1:params$lags,municipioReferencia])
base=baseCasos[(params$lags+1):length(as.matrix(baseCasos[,1])),c("Sintomas",municipioReferencia)]
colnames(base)[2]="newCases"

initialCases=as.matrix(rep(10,params$lags))
base = list("newCases"=cos(seq(10,length.out= length(as.matrix(baseCasos[,1])) - (params$lags+1))/15)*4 + 4, "Sintomas"=baseCasos$Sintomas)


result_1 = betaThompson(base=base, window=params$window, initialCases=initialCases, lags=params$lags,
            update=1, parameter1=params$parameter1 , parameter2=params$parameter2)

result_2 = betaStateSpace(base=base, initialCases=initialCases, lags=params$lags, 
              lags2=params$lags2, parameterINC1=params$parameterINC1, 
              parameterINC2=params$parameterINC2, parameterINF1=params$parameterINF1, 
              parameterINF2=params$parameterINF2)


timestamp = format(Sys.time(), "%Y%m%d%H%M")
out_dir = paste("betaStateSpace", "tuned", timestamp, sep="/")
dir.create(out_dir, recursive=TRUE)
ggsave(filename = paste(out_dir, "State-space.png", sep="/"), plot = result_1$p)
ggsave(filename = paste(out_dir, "Thompson.png", sep="/"), plot = result_2$p)
dev.off()
