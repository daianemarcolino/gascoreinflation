# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
library(MARSS)
library(zoo)
source("functions/dcs_fk_estimation.R")
source("functions/dcs_estabilidade.R")
source("functions/dcs_fcst.R")
source("functions/ERRO.R")
source("functions/core.fcst.R")

# leitura -----------
ipc <- window(readRDS("data/ipc.rds"), start = c(2001,1), freq = 12)
nucleo_tf <- readRDS("data/nucleo_tf.rds")[,3]
nucleo_dcs <- window(readRDS("data/nucleo_dcs.rds"), end = c(2017,12), freq = 12)

# previsoes dcs --------------------------------------
initial_gamma <- c(0.418235294,-0.215294118,-0.011764706,-0.048235294,-0.183529412,-0.445294118,-0.320000000,-0.401764706,-0.374705882,-0.247647059, 0.008823529)

parametros3_d <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","df","beta","mu[1|0]",paste0("gamma",1:11)          ,"psi","phi","k3","d1"),
    value = c(0.3 ,0   ,0.5 ,5   ,8   ,0     ,0.5        ,as.vector(initial_gamma)[1:11],-1    ,0.6  ,0.5 ,0   ),
    lower = c(0.0 ,0   ,0.05,-Inf,4   ,0     ,-Inf     ,rep(-Inf,11)                  ,-Inf ,-1   ,0   ,-Inf),
    upper = c(Inf ,0   ,Inf ,Inf ,Inf ,0     ,Inf      ,rep(Inf,11)                   ,Inf  ,1    ,Inf ,Inf )
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros3_d

datas <- list(ano = as.numeric(substr(as.Date(ipc),1,4))[(length(ipc)-1):length(ipc)],
              mes = as.numeric(substr(as.Date(ipc),6,7))[(length(ipc)-1):length(ipc)])

iter <- list()
for(i in 1:length(datas$ano)){
  iter[[i]] <- dcs_fcst(y = ipc, start = c(datas$ano[i], datas$mes[i]), initial = parametros3_d, type = "BSM3", outlier = T, m = 1)
  message("\n \n \n \n \n \n \n ================== date:", datas$ano[i],"/",datas$mes[i])
}

saveRDS(iter, "data/fcst_dcs_mu.rds")