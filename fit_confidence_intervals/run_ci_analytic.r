library(rlist)
library(nloptr)
library(tictoc)
library(future.apply)
library(doParallel)
library(rjson)

#setwd("scripts/funciones")
source("f_get_expected.r", encoding = "UTF-8")
source("f_fit.r", encoding = "UTF-8")
source("f_Rcaso.r", encoding = "UTF-8")
source("f_bootstrap_IC.r", encoding = "UTF-8")

## Load Data ##
#currently working on dummy data.
result <- fromJSON(file = "../case_study/fitted_models.json")
ajuste_procesoAnalitico = list("beta"=result$poisson_exp$beta_mode,
                                "N0"= 81.11421,
                                "tau1"=1.632911,
                                "tau2"=0.123283)

#incidences
ColombiaB = list("newCases"=result$I)
#### Ajustes Ygorro y tiempos ---- 
Ygorro <- get_expected(beta = ajuste_procesoAnalitico$beta,
                       N0   = ajuste_procesoAnalitico$N0,
                       tau1 = ajuste_procesoAnalitico$tau1, 
                       tau2 = ajuste_procesoAnalitico$tau2)$NN
tiempos <- length(ajuste_procesoAnalitico$beta)

#### Corroboración por suma de cuadrados ----

prueba <- data.frame(Ygorro, newCases = ColombiaB$newCases)

scm <-  sum( ( prueba$Ygorro - mean(prueba$Ygorro) )^2 )
sct <-  sum( ( prueba$newCases - mean(prueba$newCases) )^2 ) 

print("sanity check R^2:", 1 - scm/sct)
### Número de iteracacciones ---- 

num_iteracciones <- 1000
#This has to be adjusted according to JF experiment
ejecucionCastigo <- 2^12

bootstrap_distribucion_NumCasos <- bootstrap_samples( expected_I = Ygorro,
                              observed_I = ColombiaB$newCases[1:tiempos],
                              number_samples = num_iteracciones , # Aumentar a minimo 1000
                              window_size=4 , use_wa = FALSE, 
                              beta = 0.9)

bootstrap_distribucion_NumCasos[ bootstrap_distribucion_NumCasos < 0 ] <- NA

for( j in 1:ncol(bootstrap_distribucion_NumCasos) ){
  temp_j <- bootstrap_distribucion_NumCasos[ , j ]
  temp_j[is.na(temp_j)] <- median(temp_j, na.rm = T )
  bootstrap_distribucion_NumCasos[ , j ] <- temp_j 
}

# #### Do Parallel ------
#beta_ultimo <- rep(ajuste_procesoAnalitico$beta[length( ajuste_procesoAnalitico$beta)], length( ajuste_procesoAnalitico$beta)) 
beta_ultimo <- ajuste_procesoAnalitico$beta
X   <- 1:nrow(bootstrap_distribucion_NumCasos)
tmp <- vector( length = nrow(bootstrap_distribucion_NumCasos) , mode = "list" )

F_INC_SHAPE = 3.16
F_INC_RATE = 5.16
F_INF_SCALE = 24.2
tau1 = 1 / (F_INC_SHAPE / F_INC_RATE)
tau2 = 1/8.111421 # 1 / MEAN_INFECTIOUS_TIME
A0 = ColombiaB$newCases[1] / pexp(1, rate=tau1)
N0_max = A0*pexp(1, rate=tau1)/(1 * tau2) # A0 = A0*(1 - pexp(1, rate=tau1)) + N0*beta0
N0_min = A0*pexp(1, rate=tau1)/(2 * tau2)

cl   <- parallel::makeCluster( detectCores()-1 )
doParallel::registerDoParallel(cl)

tic()
lista_OptimizacionBetas <- foreach( i = X ) %dopar% {
  library(nloptr)
  # setwd(ruta_proyecto)
  # source( "scripts/funciones/f_fit.R", encoding = "UTF-8" )
  
  cat(i, " ")
  tmp2 <- fit(observed_I = as.numeric(as.matrix(bootstrap_distribucion_NumCasos[i,1:(tiempos)])),
             beta0 = beta_ultimo,
             beta_min = 0.6 * tau2, beta_max = 2.75 * tau2,
             tau10 = tau1, tau1_min = tau1,
             tau1_max = tau1,
             tau20 = tau2, tau2_min = tau2,
             tau2_max = tau2,
             N00 = N0_min, N0_min = N0_min, N0_max = N0_max,
             Castigo = ejecucionCastigo, ignore_beta_diff = NULL,
             max_eval = 1000)
}
toc()

parallel::stopCluster(cl)


####  
m_beta <- matrix(NA, nrow = nrow(bootstrap_distribucion_NumCasos),
                 ncol =  ncol(bootstrap_distribucion_NumCasos))

m_Rt   <- matrix(NA, nrow = nrow(bootstrap_distribucion_NumCasos),
                 ncol =  ncol(bootstrap_distribucion_NumCasos))

m_Rcaso <- matrix(NA, nrow = nrow(bootstrap_distribucion_NumCasos),
               ncol =  ncol(bootstrap_distribucion_NumCasos))

tic()
lags1 <- 7
par1_inf=1
distribucion_inf="gamma"
for(i in 1:nrow(bootstrap_distribucion_NumCasos)){
  par2_inf <- lista_OptimizacionBetas[[i]]$tau2
  m_beta[i,] <- lista_OptimizacionBetas[[i]]$beta
  m_Rt[i,] <- m_beta[i,] * (1 / lista_OptimizacionBetas[[i]]$tau2) 
  omegas=weightConstructionIncubacionInfeccion(lags1, par1_inf, par2_inf, distribucion_inf )  
  m_Rcaso[i,] <- R_caso(lista_OptimizacionBetas[[i]]$beta, omegas)
  rm(omegas)
 }

# Rt = beta ¨(1/ tau2) 

toc(); 

colnames(m_beta) <- colnames(bootstrap_distribucion_NumCasos)
colnames(m_Rt) <- colnames(bootstrap_distribucion_NumCasos)
colnames(m_Rcaso) <- colnames(bootstrap_distribucion_NumCasos)



######################## Intervalos de confianza ##############
confianza <- 0.95
alpha <- 1 - confianza

# Beta
Limite_inferior_Beta <- apply(m_beta, 2, FUN = quantile, probs = alpha/2)
Beta_puntual <- colMeans(m_beta)
Limite_superior_Beta <- apply(m_beta, 2, FUN = quantile, probs = 1-alpha/2)

df_estimaIC_Beta <- data.frame(LI_beta = Limite_inferior_Beta, 
                               Beta_estima = Beta_puntual,
                               LS_beta = Limite_superior_Beta)

# Rt
Limite_inferior_Rt <- apply(m_Rt, 2, FUN = quantile, probs = alpha/2)
Rt_puntual <- colMeans(m_Rt)
Limite_superior_Rt <- apply(m_Rt, 2, FUN = quantile, probs = 1-alpha/2)

df_estimaIC_Rt <- data.frame(LI_Rt = Limite_inferior_Rt, 
                             Rt_estima = Rt_puntual,
                             LS_Rt = Limite_superior_Rt)
    
print(df_estimaIC_Rt)

# # Rt Caso
Limite_inferior_RtCaso <- apply(m_Rcaso, 2, FUN = quantile, probs = alpha/2)
RtCaso_puntual <- colMeans(m_Rcaso)
Limite_superior_RtCaso <- apply(m_Rcaso, 2, FUN = quantile, probs = 1-alpha/2)

df_estimaIC_RtCaso <- data.frame(LI_RtCaso = Limite_inferior_RtCaso, 
                                 RtCaso_estima = RtCaso_puntual,
                                 LS_RtCaso = Limite_superior_RtCaso)

print(df_estimaIC_RtCaso)

saveRDS(df_estimaIC_Rt, "df_estimaIC_Rt_procesoAnalitico_.rds")
saveRDS(df_estimaIC_RtCaso, "df_estimaIC_RtCaso_procesoAnalitico_.rds")