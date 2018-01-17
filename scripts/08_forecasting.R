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

# previsão

datas <- list(ano = as.numeric(substr(as.Date(ipc),1,4))[(length(ipc)-35):length(ipc)],
              mes = as.numeric(substr(as.Date(ipc),6,7))[(length(ipc)-35):length(ipc)])

# iter <- list()
# for(i in 1:length(datas$ano)){
#   iter[[i]] <- dcs_fcst(y = ipc, start = c(datas$ano[i], datas$mes[i]), initial = parametros3_normald, type = "BSM3_normal", outlier = T, m = 2000)
#   message("\n \n \n \n \n \n \n ================== date:", datas$ano[i],"/",datas$mes[i])
# }
# saveRDS(iter, "data/fcst_dcs.rds")
iter <- readRDS("data/fcst_dcs_1ano.rds")
iter <- readRDS("data/fcst_dcs.rds")

par(mfrow = c(1,3))
ts.plot(round(na.omit(cbind(ipc, iter[[1]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
ts.plot(round(na.omit(cbind(ipc, iter[[2]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
ts.plot(round(na.omit(cbind(ipc, iter[[3]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))

# gráfico
for(i in 1:12){
  ts.plot((cbind(ipc,iter[[i]]$fcst[,1])), col = 1:2)
  Sys.sleep(0.5)
}

fcst_mean <- sapply(iter, FUN = function(x)  x$fcst[,1])
fcst_lower <- sapply(iter, FUN = function(x)  x$fcst[,2])
fcst_upper <- sapply(iter, FUN = function(x)  x$fcst[,3])

fcsts_mean  <- do.call(cbind, fcst_mean)
fcsts_lower <- do.call(cbind, fcst_lower)
fcsts_upper <- do.call(cbind, fcst_upper)

colnames(fcsts_mean) = colnames(fcsts_lower) = colnames(fcsts_upper) <- paste0("fcst", 1:ncol(fcsts_mean))

ts.plot(ipc, fcsts_lower[,1], fcsts_mean[,1], fcsts_upper[,1], col = c(1,4,2,4))
ts.plot(ipc, fcsts_lower[,2], fcsts_mean[,2], fcsts_upper[,2], col = c(1,4,2,4))
ts.plot(ipc, fcsts_lower[,3], fcsts_mean[,3], fcsts_upper[,3], col = c(1,4,2,4))


# previsoes 1 passo a frente

View(fcsts_mean)

k <- data.frame(fcsts_mean)
fcsts_df <- apply(k, MARGIN = 2, FUN = function(x) c(na.omit(x),rep(NA,sum(is.na(x)))))

prev_1frente <- ts(t(fcsts_df)[,1], end = end(fcsts_mean), freq = 12)
prev_3frente <- ts(na.omit(t(fcsts_df)[,3]), end = end(fcsts_mean), freq = 12)
prev_6frente <- ts(na.omit(t(fcsts_df)[,6]), end = end(fcsts_mean), freq = 12)
prev_12frente <- ts(na.omit(t(fcsts_df)[,12]), end = end(fcsts_mean), freq = 12)

ts.plot(na.omit(cbind(ipc,prev_1frente)), col = 1:2)
ts.plot(na.omit(cbind(ipc,prev_3frente)), col = 1:2)
ts.plot(na.omit(cbind(ipc,prev_6frente)), col = 1:2)
ts.plot(na.omit(cbind(ipc,prev_12frente)), col = 1:2)

par(mar = c(2,4,1,2), mfrow = c(1,1))

prevs <- window(cbind(ipc,prev_1frente,prev_3frente,prev_6frente,prev_12frente), start= c(2016,1), freq = 12)

plot(prevs[,1], col = 1, lwd = c(2,1,2,2), lty = c(1,1,5,4), ylim = c(-0.3,2),
     xaxt = "n", ylab = "variação mensal percentual (%)", xlab = "", type = "o")
lines(prevs[,2], col = "#1874CD", lwd = 1, lty = 1, type = "o")
lines(prevs[,3], col = "#1874CD", lwd = 1, lty = 3, type = "o")
lines(prevs[,4], col = "#1874CD", lwd = 1, lty = 5, type = "o")
lines(prevs[,5], col = "#1874CD", lwd = 1, lty = 5, type = "o")

abline(v = seq(2016,2018,0.0833), h = seq(-0.25,2,0.25),lty = 3, col = "#C9C9C9")

axis(1, at = seq(2016,2018,0.0833*3) , labels = substr(as.Date(estabilidade_tf$out),1,7)[seq(1,25,3)])
legend(2017,0.8, legend = c("mês a mês","completo"), lwd = c(1,2), lty = c(5,1),# y.intersp = 1.5,
       col = c("#1874CD","#1874CD"), cex = 1.1,bg = "white", box.col = "white",box.lwd = 0)

sqrt(mean((prevs[,1] - prevs[,2])^2))
sqrt(mean((prevs[,1] - prevs[,3])^2))
sqrt(mean((prevs[,1] - prevs[,4])^2))
sqrt(mean((prevs[,1] - prevs[,5])^2))

rmse <- matrix(NA,nrow = 49, ncol = 1)

for(i in 1:nrow(fcsts_df)){
  ychapeu <- window(ts(c(fcsts_df[i,]), start = c(datas$ano[i], datas$mes[i]), freq = 12), start = c(2016,1), freq = 12)
  d <- cbind(ipc,ychapeu)
  rmse[i,1] <- sqrt(mean((d[,1] - d[,2])^2, na.rm = T))
}

# erro linha i passos a frente
barp(rmse[c(1,3,6,12),], names.arg = c(1,3,6,12))
barp(rmse[1:12], xlab = "passos à frente")



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



