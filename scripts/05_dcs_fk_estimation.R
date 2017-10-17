# pacotes, funções e leitura
library(dygraphs)
library(BETS)
source("functions/dcs_fk_estimation.R")
source("functions/ERRO.R")


# USGDP

usgdp <- window(ts(read.csv("dados/GDPC1.csv")[,2], start = c(1947,1), freq = 4), end = c(2012,1), freq = 4)
y <- diff(log(usgdp))
m <- dcs_fk_estimation(y, initial = c(0.2,0.2,0.2,4,5), type = "mean")
round(m$otimizados$par,4)
ts.plot(y,m$out[,c(1)], col = 1:2, lwd = 1:2)


# AWHMAN
aw <- window(ts(read.csv("dados/AWHMAN.csv")[,2], start = c(1939,1), freq = 12), start = c(1992,2), end = c(2010,5), freq = 12)

m <- dcs_fk_estimation(aw, initial = c(0.2,0.2,4,5,1), type = "cn")
round(m$otimizados$par,4)

k <- window(cbind(aw,lag(m$out[,c(1)],1)), end = c(1997,1), freq = 12)
ts.plot(k, col = 1:2, lwd = 1:2)

# IPC
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
m1 <- dcs_fk_estimation(ipc, initial = c(0.2,0.2,0.2,4,5), type = "mean")
m2 <- dcs_fk_estimation(ipc, initial = c(0.2,0.2,4,5,1), type = "cn")
m3 <- dcs_fk_estimation(ipc, type = "cn2")

round(m1$otimizados$par,4)
round(m2$otimizados$par,4)
ts.plot(ipc,m1$out[,c(1)], col = 1:2, lwd = 1:2)
ts.plot(ipc,m2$out[,c(1)], col = 1:2, lwd = 1:2)
ts.plot(ipc,m3$out[,"mu"] + m3$out[,"beta"] + m3$out[,"gamma"], col = 1:2, lwd = 1:2)
ts.plot(ipc,m3$out[,"mu"], col = 1:2, lwd = 1:2)
ts.plot(ipc,m3$out[,"mu"] + m3$out[,"beta"], col = 1:2, lwd = 1:2)
ts.plot(ipc,m3$out[,"beta"], col = 1:2, lwd = 1:2)

plot(m3$out)

m <- dcs_fk_estimation(ipc, type = "cn3")
y=ipc
ts.plot(y,m[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência + Sazonalidade")
ts.plot(y,m[,"f1"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
ts.plot(y,m[,"f2"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
ts.plot(y,m[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")
plot(m[,c("f1","f2")], col = 1:2, lwd = c(1,2), main = "Tendência e Sazonalidade")
