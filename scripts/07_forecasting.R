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
source("functions/bsm.R")
source("functions/diag.dcs.R", encoding = "utf8")

# leitura -----------
ipc <- window(readRDS("data/ipc.rds"), start = c(2001,1), freq = 12)
nucleo_tf <- readRDS("data/nucleo_tf.rds")

initial_gamma <- c(0.418235294,-0.215294118,-0.011764706,-0.048235294,-0.183529412,-0.445294118,-0.320000000,-0.401764706,-0.374705882,-0.247647059, 0.008823529)

parametros3_normald <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","df","beta","mu[1|0]",paste0("gamma",1:11)          ,"psi","phi","k3","d1"),
    value = c(0.1 ,0   ,0.5 ,5   ,0   ,0     ,0        ,as.vector(initial_gamma)[1:11],1    ,0.1  ,0.5 ,0   ),
    lower = c(0.0 ,0   ,0.0 ,-Inf,0   ,0     ,-Inf     ,rep(-Inf,11)                  ,-Inf ,-1   ,0   ,-Inf),
    upper = c(Inf ,0   ,Inf ,Inf ,0   ,0     ,Inf      ,rep(Inf,11)                   ,Inf  ,1    ,Inf ,Inf )
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros3_normald

# estabilidade
estabilidade <- dcs_estabilidade(y = ipc, start = c(2013, 1), initial = parametros3_normald, type = "BSM3_normal", outlier = T)
saveRDS(estabilidade, "data/estabilidade.rds")
ts.plot(round(estabilidade$out,2), col = c(1,2), lty = c(6,3))

# previsão



datas <- list(ano = as.numeric(substr(as.Date(ipc),1,4))[(length(ipc)-3):length(ipc)],
              mes = as.numeric(substr(as.Date(ipc),6,7))[(length(ipc)-3):length(ipc)])


iter <- list()
for(i in 1:length(datas$ano)){
  iter[[i]] <- dcs_fcst(y = ipc, start = c(datas$ano[i], datas$mes[i]), initial = parametros3_normald, type = "BSM3_normal", outlier = T, m = 1)
  message("\n \n \n \n \n \n \n ================== date:", datas$ano[i],"/",datas$mes[i])
}

ts.plot(round(na.omit(cbind(ipc, iter[[1]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
ts.plot(round(na.omit(cbind(ipc, iter[[2]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
ts.plot(round(na.omit(cbind(ipc, iter[[3]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))






dcs3_normald <- dcs_fk_estimation(ipc, initial = parametros3_normald, type = "BSM3_normal", outlier = T, otimo = T)

data.frame(name = parametros3_normald$par$name, 
           lower = parametros3_normald$par$lower,
           upper = parametros3_normald$par$upper,
           initial = round(parametros3_normald$par$value,4), 
           otimo = round(dcs3_normald$otimizados$par,4))
ts.plot(ipc,dcs3_normald$out[,"mu"], col = 1:2)

ts.plot(dcs3_normald$out[,"epsilon"], col = 1)
round(dcs3_normald$out[,"epsilon"],2)
dcs3_normal$out[,"beta"]

diag_dcs3_normald <- diag.dcs(out = dcs3_normald, type = "norm")
diag_dcs3_normald$stats

# leitura
# ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
# library(BETS)
# library(zoo)
# library(plotrix)
# source('functions/dcs_fcst.R', encoding = 'UTF-8')
# source('functions/dcs_fk_estimation.R', encoding = 'UTF-8')
#
# initial_gamma <- c(0.48070729,-0.05859547, 0.15370156, 0.10032539,-0.03874261,-0.26932386,-0.12546543,-0.18050091,-0.20419881,-0.07149072, 0.06560687)
# 
# parametros <- list(
#   par = data.frame(
#     name =  c("k1","k2","ks","f2","df","beta","mu[1|0]",paste0("gamma",1:11)          ,"psi","phi","k3"),
#     value = c(0.1 ,0   ,0.5 ,5   ,6   ,0     ,0        ,as.vector(initial_gamma)[1:11],1    ,0.1  ,0.5),
#     lower = c(0.05,0   ,0.15 ,-Inf,4   ,0     ,-Inf     ,rep(-Inf,11)                  ,-Inf ,-1   ,0),
#     upper = c(Inf ,0   ,Inf ,Inf ,Inf ,0     ,Inf      ,rep(Inf,11)                   ,Inf  ,1    ,Inf)
#   )
# )
# parametros
# 
# datas <- list(ano = as.numeric(substr(as.Date(ipc),1,4))[(length(ipc)-47):length(ipc)],
#               mes = as.numeric(substr(as.Date(ipc),6,7))[(length(ipc)-47):length(ipc)])
# 
# # iter <- list()
# # for(i in 1:length(datas$ano)){
# #   iter[[i]] <- dcs_fcst(y = ipc, start = c(datas$ano[i], datas$mes[i]), initial = parametros, type = "BSM3", outlier = F, m = 3)
# #   message("\n \n \n \n \n \n \n ================== date:", datas$ano[i],"/",datas$mes[i])
# # }
# # saveRDS(iter, "iter_BSM3_2.rds")
# ts.plot(round(na.omit(cbind(ipc, iter[[1]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
# ts.plot(round(na.omit(cbind(ipc, iter[[2]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
# ts.plot(round(na.omit(cbind(ipc, iter[[3]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
# 
# iter <- readRDS("dados/iter_BSM3.rds") 
# 
# # gráfico
# for(i in 1:49){
#   ts.plot((cbind(ipc,iter[[i]]$fcst[,1])), col = 1:2)
#   Sys.sleep(0.5)
# }
# 
# fcst_mean <- sapply(iter, FUN = function(x)  x$fcst[,1])
# fcst_lower <- sapply(iter, FUN = function(x)  x$fcst[,2])
# fcst_upper <- sapply(iter, FUN = function(x)  x$fcst[,3])
# 
# fcsts_mean  <- do.call(cbind, fcst_mean)
# fcsts_lower <- do.call(cbind, fcst_lower)
# fcsts_upper <- do.call(cbind, fcst_upper)
# 
# colnames(fcsts_mean) = colnames(fcsts_lower) = colnames(fcsts_upper) <- paste0("fcst", 1:ncol(fcsts_mean))
# 
# ts.plot(ipc, fcsts_lower[,1], fcsts_mean[,1], fcsts_upper[,1], col = c(1,4,2,4))
# ts.plot(ipc, fcsts_lower[,36], fcsts_mean[,36], fcsts_upper[,36], col = c(1,4,2,4))
# ts.plot(ipc, fcsts_lower[,40], fcsts_mean[,40], fcsts_upper[,40], col = c(1,4,2,4))
# 
# # previsoes 1 passo a frente
# 
# View(fcsts_mean)
# 
# k <- data.frame(fcsts_mean)
# fcsts_df <- apply(k, MARGIN = 2, FUN = function(x) c(na.omit(x),rep(NA,sum(is.na(x)))))
# 
# rmse <- matrix(NA,nrow = 49, ncol = 1)
# 
# for(i in 1:nrow(fcsts_df)){
#   ychapeu <- ts(c(fcsts_df[i,]), start = c(datas$ano[i], datas$mes[i]), freq = 12)
#   d <- cbind(ipc,ychapeu)
#   rmse[i,1] <- sqrt(mean((d[,1] - d[,2])^2, na.rm = T))
# }
# 
# # erro linha i passos a frente
# barp(rmse[c(1,3,6,9,12),], names.arg = c(1,3,6,9,12))
# barp(rmse)


