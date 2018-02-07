# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
library(zoo)
source("functions/dcs_fk_estimation.R")
source("functions/dcs_estabilidade.R")
source("functions/dcs_fcst.R")
source("functions/ERRO.R")
source("functions/core.fcst.R")
source("functions/ifa_functions.R")

# leitura -----------
ipc <- window(readRDS("data/ipc.rds"), start = c(2001,1), freq = 12)
ipca <- window(readRDS("data/ipca.rds"), start = c(2001,1), freq = 12)
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

new_par <- c(0.119463244875724,0,0.05,-1.3745462986007,23.0409659173607,0,0.692015745581748,0.599825084746369,-0.0365868971327083,0.113220352342134,0.134774809026807,-0.0234199705179038,-0.298816328489241,-0.166058473073769,-0.2254689858412,-0.222905124231499,-0.101314231111584,0.0470387659877161,-0.662238154975964,0.544679342032902,0.538650742581951,1.75903961338264)
parametros3_d$par$value <- new_par
# fcst <- dcs_fcst(y = ipc, start = c(2018,1), initial = parametros3_d, type = "BSM3", outlier = T, m = 2000, h = 12, out = T, dist = "t")
# saveRDS(fcst, "data/fcst_ipc_2018_v2.rds")
fcst <- readRDS("data/fcst_ipc_2018.rds")
#fcst <- readRDS("data/fcst_ipc_2018.rds_v2")
par(mar = c(2,2,2,1))

k1 <- cbind(mean = fcst$fcst_y[,1], low = fcst$fcst_y[,2], up = fcst$fcst_y[,3])
#k1[nrow(k1)-11,1] <-tail(k1[,2],12)[1]

(prod(k1[,1]/100+1)-1)*100


plot(c(k1[,1]), type = "l", ylim = c(-0.6,1.5), xaxt = "n", ylab = "")
axis(1, at = 1:nrow(k1) , labels = 1:12)
polygon(c(1:nrow(k1), rev(1:nrow(k1))),c(k1[,2],rev(k1[,3])),col = "#C6E2FF", border = FALSE)
abline(v = 1:12, h = seq(-0.5,2,0.5), lty = 3, col = "#C9C9C9")
lines(c(k1[,1]), lty = 1, col = 2, lwd = 2)
lines(c(k1[,2]), col = "#75A1D0", lwd = 1, ylab = "")
lines(c(k1[,3]), col = "#75A1D0", lwd = 1, ylab = "")
legend(7,1.5,legend = c("Previsão", "IC de 95%"), col = c(2,"#C6E2FF"), 
       pch = c(NA,15), lty = c(1,0), lwd = c(2,0),
       cex = 1.2, pt.cex = 2.5, bg = "white", box.col = "white", box.lwd = 0)

round(k1,2)

k1 <- round(anual(cbind(mean = fcst$fcst_mu[,1], low = fcst$fcst_mu[,2], up = fcst$fcst_mu[,3])),2)

plot(c(k1[,1]), type = "l", ylim = c(2.5,7.5), xaxt = "n", ylab = "")
axis(1, at = 1:nrow(k1) , labels = 1:12)
polygon(c(1:nrow(k1), rev(1:nrow(k1))),c(k1[,2],rev(k1[,3])),col = "#C6E2FF", border = FALSE)
abline(v = 1:12, h = seq(-0.5,2,0.5), lty = 3, col = "#C9C9C9")
lines(c(k1[,1]), lty = 1, col = 2, lwd = 2)
lines(c(k1[,2]), col = "#75A1D0", lwd = 1, ylab = "")
lines(c(k1[,3]), col = "#75A1D0", lwd = 1, ylab = "")
legend(1,7.5,legend = c("Previsão", "IC de 95%"), col = c(2,"#C6E2FF"), 
       pch = c(NA,15), lty = c(1,0), lwd = c(2,0),
       cex = 1.2, pt.cex = 2.5, bg = "white", box.col = "white", box.lwd = 0)

(prod((k1[,3]/100+1)^(1/12))-1)*100


# PREVISAO IPCA ------------------------------------
k1 <- cbind(mean = fcst$fcst_y[,1], low = fcst$fcst_y[,2], up = fcst$fcst_y[,3])
x <- cbind(IPCA = ipca, IPC = ipc, IPC1 = lag(ipc,-1), IPC12 = lag(ipc,-12), IPCA12 = lag(ipca,-12))
x[(nrow(x)-11):nrow(x),"IPC"] <- k1[,1]
x[(nrow(x)-10):nrow(x),"IPC1"] <- k1[-12,1]
x

x_in <- window(x, start = c(2002,1), end = c(2017,12), freq = 12)
x_out <- window(x, start = c(2018,1), end = c(2018,12), freq = 12)[,-1]

xlow <- cbind(IPCA = ipca, IPC = ipc, IPC1 = lag(ipc,-1), IPC12 = lag(ipc,-12), IPCA12 = lag(ipca,-12))
xlow[(nrow(xlow)-11):nrow(xlow),"IPC"] <- k1[,2]
xlow[(nrow(xlow)-10):nrow(xlow),"IPC1"] <- k1[-12,2]
xlow

xlow_in <- window(xlow, start = c(2002,1), end = c(2017,12), freq = 12)
xlow_out <- window(xlow, start = c(2018,1), end = c(2018,12), freq = 12)[,-1]

xup <- cbind(IPCA = ipca, IPC = ipc, IPC1 = lag(ipc,-1), IPC12 = lag(ipc,-12), IPCA12 = lag(ipca,-12))
xup[(nrow(xup)-11):nrow(xup),"IPC"] <- k1[,3]
xup[(nrow(xup)-10):nrow(xup),"IPC1"] <- k1[-12,3]
xup

xup_in <- window(xup, start = c(2002,1), end = c(2017,12), freq = 12)
xup_out <- window(xup, start = c(2018,1), end = c(2018,12), freq = 12)[,-1]

m <- lm(IPCA ~ ., data = x_in)
summary(m)
acf(resid(m), 48)
plot(resid(m))
ipca_fcst <- ts(predict(m, newdata = x_out), start = c(2018,1), freq = 12)

m <- lm(IPCA ~ ., data = xlow_in)
ipca_fcst_low <- ts(predict(m, newdata = xlow_out), start = c(2018,1), freq = 12)

m <- lm(IPCA ~ ., data = xup_in)
ipca_fcst_up <- ts(predict(m, newdata = xup_out), start = c(2018,1), freq = 12)


ts.plot(ipc, ipca, ipca_fcst, k1[,1], col = c(1,2,2,1))

(prod(ipca_fcst/100+1)-1)*100


plot(c(ipca_fcst), type = "l", ylim = c(-0.7,1.2),xaxt = "n", ylab = "")
axis(1, at = 1:length(ipca_fcst) , labels = 1:12)
polygon(c(1:length(ipca_fcst), rev(1:length(ipca_fcst))),c(ipca_fcst_low,rev(ipca_fcst_up)),col = "#C6E2FF", border = FALSE)
abline(v = 1:12, h = seq(-0.5,2,0.5), lty = 3, col = "#C9C9C9")
lines(c(ipca_fcst), lty = 1, col = 2, lwd = 2)
lines(c(ipca_fcst_low), col = "#75A1D0", lwd = 1, ylab = "")
lines(c(ipca_fcst_up), col = "#75A1D0", lwd = 1, ylab = "")
legend(1,-0.3,legend = c("Previsão", "IC de 95%"), col = c(2,"#C6E2FF"), 
       pch = c(NA,15), lty = c(1,0), lwd = c(2,0),
       cex = 1.2, pt.cex = 2.5, bg = "white", box.col = "white", box.lwd = 0)

