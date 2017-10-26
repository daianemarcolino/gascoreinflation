# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
source("functions/dcs_fk_estimation.R")
source("functions/ERRO.R")

# leitura -----------
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
plot(ipc)

turistas <- log(ts(read.csv2("dados/turistas.csv")[,2], start = c(2000,1), freq = 12))
plot(turistas)
saveRDS(turistas, "shiny_simulacao/data/turistas.rds")
# modelos -----------

# replicar 
bsm <- dcs_fk_estimation(turistas, type = "BSM_artigo")
saveRDS(bsm, "shiny_simulacao/data/out_turistas.rds")

par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(turistas,bsm$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
ts.plot(bsm$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
ts.plot(bsm$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
ts.plot(bsm$out[,c("scorestd","epsilonstd")], col = 1, lty = c(1,3), main = "Score e Epsilon padr.", ylab = "")

ts.plot(bsm$out[,c("epsilonstd")], col = 1, lty = c(1,3), main = "Epsilon padr.", ylab = "", ylim = c(-3,3))
ts.plot(bsm$out[,c("scorestd")], col = 1, lty = c(1,3), main = "Score padr.", ylab = "", ylim = c(-3,3))
acf(bsm$out[,c("epsilonstd")], 20, drop.lag.0 = T, main = "Epsilon padr.")
acf(bsm$out[,c("scorestd")], 20, drop.lag.0 = T, main = "Score padr.")


# modelo BSM
bsm <- dcs_fk_estimation(ipc, initial = c(0.2,0.2,0.2,5), type = "BSM")

ts.plot(ipc,bsm$out[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
ts.plot(ipc,bsm$out[,"beta"], col = 1:2, lwd = c(1,2), main = "IPC e Beta")
ts.plot(ipc,bsm$out[,"gamma"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
ts.plot(ipc,bsm$out[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")





















# # USGDP
# 
# usgdp <- window(ts(read.csv("dados/GDPC1.csv")[,2], start = c(1947,1), freq = 4), end = c(2012,1), freq = 4)
# y <- diff(log(usgdp))
# m <- dcs_fk_estimation(y, initial = c(0.2,0.2,0.2,4,5), type = "mean")
# round(m$otimizados$par,4)
# ts.plot(y,m$out[,c(1)], col = 1:2, lwd = 1:2)
# 
# 
# # AWHMAN
# aw <- window(ts(read.csv("dados/AWHMAN.csv")[,2], start = c(1939,1), freq = 12), start = c(1992,2), end = c(2010,5), freq = 12)
# 
# m <- dcs_fk_estimation(aw, initial = c(0.2,0.2,4,5,1), type = "cn")
# round(m$otimizados$par,4)
# 
# k <- window(cbind(aw,lag(m$out[,c(1)],1)), end = c(1997,1), freq = 12)
# ts.plot(k, col = 1:2, lwd = 1:2)
# 
# # IPC
# ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
# m1 <- dcs_fk_estimation(ipc, initial = c(0.2,0.2,0.2,4,5), type = "mean")
# m2 <- dcs_fk_estimation(ipc, initial = c(0.2,0.2,4,5,1), type = "cn")
# m3 <- dcs_fk_estimation(ipc, type = "cn2")
# 
# round(m1$otimizados$par,4)
# round(m2$otimizados$par,4)
# ts.plot(ipc,m1$out[,c(1)], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m2$out[,c(1)], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m3$out[,"mu"] + m3$out[,"beta"] + m3$out[,"gamma"], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m3$out[,"mu"], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m3$out[,"mu"] + m3$out[,"beta"], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m3$out[,"beta"], col = 1:2, lwd = 1:2)
# 
# plot(m3$out)
# 
# # initial dummies
# 
# dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 11),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 12)
# )
# 
# y_filter <- y[1:13]
# dummy_filter <- dummy[1:13,]
# t <- 1:13
# dados <- data.frame(cbind(y_filter,t,dummy_filter))
# dados <- data.frame(cbind(y,1:length(ipc),dummy))
# colnames(dados) <- c("y","t",paste0("D",1:12))
# md <- lm(y ~ ., data = dados)
# 
# data <- data.frame(param = c(paste0("D",1:11),            "df","w1","w2","beta","a11","a12","a13","a2","b2","f1","f2"),
#                    value = c(c(md$coefficients[-c(1,2,14)]), 3,   1,   1,  0.10,  0.5,  0.5,  0.5, 0.5, 0.9,   2, 0.5),
#                    lower = c(rep(-Inf,11),                   2,-Inf,-Inf, -Inf,-Inf,-Inf,-Inf,-Inf,-1,  -Inf, -Inf),
#                    upper = c(rep( Inf,11),                 Inf, Inf, Inf,  Inf, Inf, Inf, Inf, Inf, 1,   Inf,  Inf))   
# 
# data$value <- round(c(1.395931e+00, 7.717124e-01, 1.143014e+00, 1.017791e+00,-9.203240e+00, 7.213934e-01, 8.484142e-01, 7.339786e-01, 7.438413e-01, 
#  7.914208e-01, 9.271763e-01, 3.920452e+00,-8.925764e-01,-8.142339e-02,-3.349298e-03, 7.269610e-03, -7.303699e-05, 6.873406e-03,
#  2.295186e-01, 9.511262e-01, 6.880529e-01, 1.036890e+00),4)
# 
# data$value <- m$otimizados$par
# 
# round(m$otimizados$par[1:11], 3)
# round(md$coefficients[-c(1,2,14)], 3)
# 
# data_new <- data.frame(param = c(paste0("D",1:11),            "df","w1","w2","beta","a11","a12","a13","a2","b2","f1","f2"),
#                    value = c(m$otimizados$par))   
# 
# 
# data <- data.frame(param = c(paste0("D",1:11),"df","w1","w2","beta","a11","a12","a13","a2","b2","f1","f2"),
#                    value = c(rep(1,11),          4,   0, 0.1,  0.10,  0.3,  0.3,  0.5, 0.5, 0.1, 0.5, 0.5),
#                    lower = c(rep(-Inf,11),       2,-Inf,-Inf, -Inf,-Inf,-Inf,-Inf,-Inf,-1,  -Inf, -Inf),
#                    upper = c(rep( Inf,11),     Inf, Inf, Inf,  Inf, Inf, Inf, Inf, Inf, 1,   Inf,  Inf))   
# 
# 
# 
# a <- Sys.time()
# m <- dcs_fk_estimation(ipc, type = "cn3", initial = data, dummy = dummy)
# b <- Sys.time()
# b-a
# y=ipc
# ts.plot(y,m$out[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência + Sazonalidade")
# ts.plot(y,m$out[,"f1"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
# ts.plot(y,m$out[,"f2"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
# ts.plot(y,m$out[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")
# ts.plot(m[,c("f1","f2")], col = 1:2, lwd = c(1,2), main = "Tendência e Sazonalidade")
