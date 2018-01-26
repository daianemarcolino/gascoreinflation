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
# initial_gamma <- c(0.418235294,-0.215294118,-0.011764706,-0.048235294,-0.183529412,-0.445294118,-0.320000000,-0.401764706,-0.374705882,-0.247647059, 0.008823529)
# 
# parametros3_normald <- list(
#   par = data.frame(
#     name =  c("k1","k2","ks","f2","df","beta","mu[1|0]",paste0("gamma",1:11)          ,"psi","phi","k3","d1"),
#     value = c(0.1 ,0   ,0.5 ,5   ,0   ,0     ,0        ,as.vector(initial_gamma)[1:11],1    ,0.1  ,0.5 ,0   ),
#     lower = c(0.0 ,0   ,0.0 ,-Inf,0   ,0     ,-Inf     ,rep(-Inf,11)                  ,-Inf ,-1   ,0   ,-Inf),
#     upper = c(Inf ,0   ,Inf ,Inf ,0   ,0     ,Inf      ,rep(Inf,11)                   ,Inf  ,1    ,Inf ,Inf )
#   ),
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# parametros3_normald
# 
# # previsão
# 
# datas <- list(ano = as.numeric(substr(as.Date(ipc),1,4))[(length(ipc)-35):length(ipc)],
#               mes = as.numeric(substr(as.Date(ipc),6,7))[(length(ipc)-35):length(ipc)])
#
# iter <- list()
# for(i in 1:length(datas$ano)){
#   iter[[i]] <- dcs_fcst(y = ipc, start = c(datas$ano[i], datas$mes[i]), initial = parametros3_normald, type = "BSM3_normal", outlier = T, m = 2000)
#   message("\n \n \n \n \n \n \n ================== date:", datas$ano[i],"/",datas$mes[i])
# }
# saveRDS(iter, "data/fcst_dcs.rds")
iter <- readRDS("data/fcst_dcs_1ano.rds")
iter <- readRDS("data/fcst_dcs.rds")
iter <- readRDS("data/fcst_dcs_t_3anos.rds")
iter <- readRDS("data/fcst_dcs_t_3anos_exogen.rds")

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
k <- data.frame(fcsts_upper)
fcsts_df_up <- apply(k, MARGIN = 2, FUN = function(x) c(na.omit(x),rep(NA,sum(is.na(x)))))
k <- data.frame(fcsts_lower)
fcsts_df_low <- apply(k, MARGIN = 2, FUN = function(x) c(na.omit(x),rep(NA,sum(is.na(x)))))


prev_1frente <- window(ts(t(fcsts_df)[,1], end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_3frente <- window(ts(na.omit(t(fcsts_df)[,3]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_6frente <- window(ts(na.omit(t(fcsts_df)[,6]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_12frente <- window(ts(na.omit(t(fcsts_df)[,12]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)

prev_1frente_low <- window(ts(t(fcsts_df_low)[,1], end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_3frente_low <- window(ts(na.omit(t(fcsts_df_low)[,3]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_6frente_low <- window(ts(na.omit(t(fcsts_df_low)[,6]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_12frente_low <- window(ts(na.omit(t(fcsts_df_low)[,12]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)

prev_1frente_up <- window(ts(t(fcsts_df_up)[,1], end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_3frente_up <- window(ts(na.omit(t(fcsts_df_up)[,3]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_6frente_up <- window(ts(na.omit(t(fcsts_df_up)[,6]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)
prev_12frente_up <- window(ts(na.omit(t(fcsts_df_up)[,12]), end = end(fcsts_mean), freq = 12), start = c(2016,1), freq = 12)

sqrt(mean((ipc - prev_1frente)^2))
sqrt(mean((ipc - prev_3frente)^2))
sqrt(mean((ipc - prev_6frente)^2))
sqrt(mean((ipc - prev_12frente)^2))


par(mar = c(2,2,2,1), mfrow = c(2,2))

k1 <- na.omit(cbind(ipc, mean = prev_1frente, low = prev_1frente_low, up = prev_1frente_up))
k3 <- na.omit(cbind(ipc, mean = prev_3frente, low = prev_3frente_low, up = prev_3frente_up))
k6 <- na.omit(cbind(ipc, mean = prev_6frente, low = prev_6frente_low, up = prev_6frente_up))
k12 <- na.omit(cbind(ipc, mean = prev_12frente, low = prev_12frente_low, up = prev_12frente_up))

plot(1:nrow(k1), k[,2], type = "l", ylim = c(-0.5,2), xaxt = "n", main = "1", ylab = "")
axis(1, at = 1:nrow(k1) , labels = substr(as.Date(k1),1,7))
polygon(c(1:nrow(k1), rev(1:nrow(k1))),c(k1[,3],rev(k1[,4])),col = "#C6E2FF", border = FALSE)
abline(v = 1:24, lty = 3, col = "#C9C9C9")
lines(c(k1[,1]), lty = 4, ylab = "")
lines(c(k1[,2]), col = 2, lwd = 2, ylab = "")
lines(c(k1[,3]), col = "#75A1D0", lwd = 1, ylab = "")
lines(c(k1[,4]), col = "#75A1D0", lwd = 1, ylab = "")
legend(17,2,legend = c("IPC-Br","Previsão"), col = c(1,2), lwd = c(1,2), lty = c(4,1), 
       cex = 1.2, bg = "white", box.col = "white", box.lwd = 0)

plot(1:nrow(k3), k3[,2], type = "l", ylim = c(-0.5,2), xaxt = "n", main = "3", ylab = "")
axis(1, at = 1:nrow(k3) , labels = substr(as.Date(k3),1,7))
polygon(c(1:nrow(k3), rev(1:nrow(k3))),c(k3[,3],rev(k3[,4])),col = "#C6E2FF", border = FALSE)
abline(v = 1:24, lty = 3, col = "#C9C9C9")
lines(c(k3[,1]), lty = 4, ylab = "")
lines(c(k3[,2]), col = 2, lwd = 2, ylab = "")
lines(c(k3[,3]), col = "#75A1D0", lwd = 1, ylab = "")
lines(c(k3[,4]), col = "#75A1D0", lwd = 1, ylab = "")
legend(17,2,legend = c("IPC-Br","Previsão"), col = c(1,2), lwd = c(1,2), lty = c(4,1), 
       cex = 1.2, bg = "white", box.col = "white", box.lwd = 0)

plot(1:nrow(k6), k6[,2], type = "l", ylim = c(-0.5,2), xaxt = "n", main = "6", ylab = "")
axis(1, at = 1:nrow(k6) , labels = substr(as.Date(k6),1,7))
polygon(c(1:nrow(k6), rev(1:nrow(k6))),c(k6[,3],rev(k6[,4])),col = "#C6E2FF", border = FALSE)
abline(v = 1:24, lty = 3, col = "#C9C9C9")
lines(c(k6[,1]), lty = 4, ylab = "")
lines(c(k6[,2]), col = 2, lwd = 2, ylab = "")
lines(c(k6[,3]), col = "#75A1D0", lwd = 1, ylab = "")
lines(c(k6[,4]), col = "#75A1D0", lwd = 1, ylab = "")
legend(17,2,legend = c("IPC-Br","Previsão"), col = c(1,2), lwd = c(1,2), lty = c(4,1), 
       cex = 1.2, bg = "white", box.col = "white", box.lwd = 0)

plot(1:nrow(k12), k12[,2], type = "l", ylim = c(-0.5,2), xaxt = "n", main = "12", ylab = "")
axis(1, at = 1:nrow(k12) , labels = substr(as.Date(k12),1,7))
polygon(c(1:nrow(k12), rev(1:nrow(k12))),c(k12[,3],rev(k12[,4])),col = "#C6E2FF", border = FALSE)
abline(v = 1:24, lty = 3, col = "#C9C9C9")
lines(c(k12[,1]), lty = 4, ylab = "")
lines(c(k12[,2]), col = 2, lwd = 2, ylab = "")
lines(c(k12[,3]), col = "#75A1D0", lwd = 1, ylab = "")
lines(c(k12[,4]), col = "#75A1D0", lwd = 1, ylab = "")
legend(17,2,legend = c("IPC-Br","Previsão"), col = c(1,2), lwd = c(1,2), lty = c(4,1), 
       cex = 1.2, bg = "white", box.col = "white", box.lwd = 0)


# RMSE

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


# PREVISÕES DE DIFERENÇA -----------------------

# nucleo-dcs
fcst1_dcs <- core.fcst(ipc, nucleo_dcs, lag = 1, h = 12)
fcst3_dcs <- core.fcst(ipc, nucleo_dcs, lag = 3, h = 12)
fcst6_dcs <- core.fcst(ipc, nucleo_dcs, lag = 6, h = 12)
fcst12_dcs <- core.fcst(ipc, nucleo_dcs, lag = 12, h = 12)

round(sqrt(mean((ipc - lag(ipc, -1) - fcst1_dcs$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -3) - fcst3_dcs$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -6) - fcst6_dcs$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -12) - fcst12_dcs$fcst)^2)),2)

# nucleo-s
fcst1_tf <- core.fcst(ipc, nucleo_tf, lag = 1, h = 12)
fcst3_tf <- core.fcst(ipc, nucleo_tf, lag = 3, h = 12)
fcst6_tf <- core.fcst(ipc, nucleo_tf, lag = 6, h = 12)
fcst12_tf <- core.fcst(ipc, nucleo_tf, lag = 12, h = 12)

round(sqrt(mean((ipc - lag(ipc, -1)  - fcst1_tf$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -3)  - fcst3_tf$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -6)  - fcst6_tf$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -12) - fcst12_tf$fcst)^2)),2)

ts.plot(na.omit(cbind(fcst1_dcs$fcst, fcst1_tf$fcst, ipc - lag(ipc, -1))))
ts.plot(na.omit(cbind(fcst3_dcs$fcst, fcst3_tf$fcst, ipc - lag(ipc, -3))))
ts.plot(na.omit(cbind(fcst6_dcs$fcst, fcst6_tf$fcst, ipc - lag(ipc, -6))))
ts.plot(na.omit(cbind(fcst12_dcs$fcst, fcst12_tf$fcst, ipc - lag(ipc, -12))))

ipc0 <- window(ipc, start = c(2009,7), freq = 12)
nucleo_tf0 <- window(nucleo_tf, start = c(2009,7), freq = 12)
nucleo_dcs0 <- window(nucleo_dcs, start = c(2009,7), freq = 12)

# nucleo-dcs
fcst1_dcs <- core.fcst(ipc0, nucleo_dcs0, lag = 1, h = 12)
fcst3_dcs <- core.fcst(ipc0, nucleo_dcs0, lag = 3, h = 12)
fcst6_dcs <- core.fcst(ipc0, nucleo_dcs0, lag = 6, h = 12)
fcst12_dcs <- core.fcst(ipc0, nucleo_dcs0, lag = 12, h = 12)

round(sqrt(mean((ipc0 - lag(ipc0, -1) - fcst1_dcs$fcst)^2)),2)
round(sqrt(mean((ipc0 - lag(ipc0, -3) - fcst3_dcs$fcst)^2)),2)
round(sqrt(mean((ipc0 - lag(ipc0, -6) - fcst6_dcs$fcst)^2)),2)
round(sqrt(mean((ipc0 - lag(ipc0, -12) - fcst12_dcs$fcst)^2)),2)

# nucleo-s
fcst1_tf <- core.fcst(ipc0, nucleo_tf0, lag = 1, h = 12)
fcst3_tf <- core.fcst(ipc0, nucleo_tf0, lag = 3, h = 12)
fcst6_tf <- core.fcst(ipc0, nucleo_tf0, lag = 6, h = 12)
fcst12_tf <- core.fcst(ipc0, nucleo_tf0, lag = 12, h = 12)

round(sqrt(mean((ipc0 - lag(ipc0, -1)  - fcst1_tf$fcst)^2)),2)
round(sqrt(mean((ipc0 - lag(ipc0, -3)  - fcst3_tf$fcst)^2)),2)
round(sqrt(mean((ipc0 - lag(ipc0, -6)  - fcst6_tf$fcst)^2)),2)
round(sqrt(mean((ipc0 - lag(ipc0, -12) - fcst12_tf$fcst)^2)),2)

ts.plot(na.omit(cbind(fcst1_dcs$fcst, fcst1_tf$fcst, ipc - lag(ipc, -1))))
ts.plot(na.omit(cbind(fcst3_dcs$fcst, fcst3_tf$fcst, ipc - lag(ipc, -3))))
ts.plot(na.omit(cbind(fcst6_dcs$fcst, fcst6_tf$fcst, ipc - lag(ipc, -6))))
ts.plot(na.omit(cbind(fcst12_dcs$fcst, fcst12_tf$fcst, ipc - lag(ipc, -12))))


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



