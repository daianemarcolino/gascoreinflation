# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
library(MARSS)
source("functions/dcs_fk_estimation.R")
source("functions/ERRO.R")
source("functions/bsm.R")
source("functions/diag.dcs.R", encoding = "utf8")

# leitura -----------
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
plot(ipc)

turistas <- log(ts(read.csv2("dados/turistas.csv")[,2], start = c(2000,1), freq = 12))
plot(turistas)
#saveRDS(turistas, "shiny_simulacao/data/turistas.rds")

# replicar artigo ------------------------------------------------- 
bsm0 <- dcs_fk_estimation(turistas, type = "BSM_artigo")

pseudo.y <- bsm0$out[,"mu"] + bsm0$out[,"gamma"] + bsm0$out[,"u"]
ts.plot(turistas, pseudo.y, col = 1:2)
# mu smooth
fk1 <- bsm(pseudo.y, beta = T, iter = 10000)

# fig 6
par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(turistas,bsm0$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
ts.plot(bsm0$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
ts.plot(bsm0$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
abline(h=0, lty = 3, col = 2)
ts.plot(bsm0$out[,c("u")],bsm0$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")

# fig 7
par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(bsm0$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-0.2,0.2))
ts.plot(bsm0$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-0.2,0.2))
acf(bsm0$out[,c("nu")], 20, drop.lag.0 = T, main = "acf nu")
acf(bsm0$out[,c("u")], 20, drop.lag.0 = T, main = "acf u")

# fig 8
ts.plot(fk1[,"mu"],bsm0$out[,"mu"], col = 1, lty = c(1,3))

# diagnóstico
diag <- diag.dcs(y = bsm0$out[,"epsilon"], df = bsm0$otimizados$par[4])
diag$resid.q
diag$stats
plot(rq,bsm0$out[,"epsilon"])

# saveRDS(bsm0, "shiny_simulacao/data/bsm_turistas.rds")
# saveRDS(pseudo.y, "shiny_simulacao/data/pseudoy_turistas.rds")
# saveRDS(fk1, "shiny_simulacao/data/smooth_turistas.rds")
# saveRDS(diag, "shiny_simulacao/data/diag_turistas.rds")

# BSM para IPC ----------------------------------------------------


ipc0 <- window(ipc,start = c(2005,1),freq = 12)
initial_mu <- mean(ipc0)
gammas <- data.frame(ipc0, cycle(ipc0))
initial_gamma <- tapply(gammas[,1],gammas[,2], FUN = mean)

parametros <- list(
  par = data.frame(
    name = c("k1","ks","f2","df","mu[0]",paste0("gamma",1:11)),
    value = c(0.4016,-1.1053,-1.6886,10,mean(ipc),initial_gamma[1:11]),
    lower = c(0,-Inf,-Inf,2, rep(-Inf,12)),
    upper = c(1,Inf,Inf,Inf, rep(Inf,12))
  )
)
parametros
bsm2 <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2", initial.optim = T)
#bsm2 <- readRDS("bsm2.rds")
saveRDS(bsm2, "dados/bsm_ipc.rds")

comparar <- cbind(parametros$par, bsm2$otimizados$par)
comparar

pseudo.y <- bsm2$out[,"mu"] + bsm2$out[,"gamma"] + bsm2$out[,"u"]
ts.plot(ipc, pseudo.y, col = 1:2)
# mu smooth
fk2 <- bsm(pseudo.y, beta = T, iter = 10000)

# fig 6
par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(ipc,bsm2$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
ts.plot(bsm2$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
ts.plot(bsm2$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
abline(h=0, lty = 3, col = 2)
ts.plot(bsm2$out[,c("u")],bsm2$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")

# fig 7
par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(bsm2$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
ts.plot(bsm2$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
acf(bsm2$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
acf(bsm2$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")

# fig 8
ts.plot(fk2[,"mu"],bsm2$out[,"mu"], col = 1, lty = c(1,3))
ts.plot(fk2[,"mu"],ipc, col = 1, lty = c(1,3))
ts.plot(bsm2$out[,"mu"],ipc, col = 1, lty = c(1,3))

# diagnóstico
diag <- diag.dcs(y = bsm2$out[,"epsilon"], df = bsm2$otimizados$par[4])
hist(diag$resid.q)
diag$stats
plot(rq,bsm2$out[,"epsilon"])


#par(mfrow = c(2,2))
ts.plot(ipc,bsm2$out[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
ts.plot(bsm2$out[,"mu"], col = 1:2, lwd = c(1,2), main = "Tendência")
ts.plot(ipc,bsm2$out[,"gamma"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
#ts.plot(ipc,bsm2$out[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")
#ts.plot(bsm2$out[,c("epsilon")], col = 1:2, lwd = c(1,2), main = "epsilon")
#ts.plot(bsm2$out[,c("score")], col = 1:2, lwd = c(1,2), main = "score")
ts.plot(bsm2$out[,c("score","epsilon")], col = 1, lty = c(1,3), main = "Score e Epsilon", ylab = "")
legend("top", legend = c("score","epsilon"), lty = c(1,3), col = 1, bty = "n", cex = 1)

#saveRDS(list(bsm1 = bsm1, bsm2 = bsm2), "dcsipc_novosresultados.rds")

bsm1 <- readRDS("dcsipc_novosresultados.rds")$bsm1
bsm2 <- readRDS("dcsipc_novosresultados.rds")$bsm2

hist(bsm2$out[,"score"], main = "score")

betasuposto <- rbeta(length(bsm2$out[,"score"]), 0.5, bsm2$otimizados$par[4]/2)
q1 <- quantile(bsm2$out[,"score"], probs = seq(0,1,length.out = length(ipc)))
q2 <- quantile(betasuposto, probs = seq(0,1,length.out = length(ipc)))
plot(q1,q2, xlim = c(-0.4,0.4), ylim = c(0,0.4), xlab = "quantis score", ylab = "quantis beta(1/2,v/2)")
qqline(data.frame(q1,q2))
qqplot(bsm2$out[,"score"],betasuposto, xlim = c(-0.4,0.4), ylim = c(0,0.5))
qqline(data.frame(bsm2$out[,"score"],betasuposto))

pseudo.y <- bsm2$out[,"mu"] + bsm2$out[,"gamma"] + bsm2$out[,"score"]
ts.plot(ipc, pseudo.y, col = 1:2)
fk <- bsm(pseudo.y)


dygraph(cbind(bsm2$out[,"mu"], fk$filter[,"mu"]))
dygraph(cbind(fk$filter[,"mu"],fk$smooth[,"mu"]))
dygraph(cbind(ipc,fk$filter[,"mu"],fk$smooth[,"mu"]))
dygraph(cbind(ipc,fk$smooth[,"mu"]))

dados <- cbind(ipc, bsm2$out, fk$smooth)
colnames(dados) <- c("ipc",colnames(bsm2$out),"mu_smooth","gamma_smooth")

dados_save <- data.frame(data = zoo::as.Date(dados), dados)
write.csv2(dados_save, "ultimos_resultados.csv", row.names = F)
saveRDS(dados_save, "ultimos_resultados.rds")
dygraph(((dados/100+1)^12-1)*100)

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
