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
fk1 <- bsm(pseudo.y, beta = T, iter = 9)

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

saveRDS(bsm0, "shiny_simulacao/data/bsm_turistas.rds")
saveRDS(pseudo.y, "shiny_simulacao/data/pseudoy_turistas.rds")
saveRDS(fk1, "shiny_simulacao/data/smooth_turistas.rds")
saveRDS(diag, "shiny_simulacao/data/diag_turistas.rds")

# BSM para IPC ----------------------------------------------------


ipc0 <- window(ipc,start = c(2004,1),freq = 12)
initial_mu <- median(ipc)
gammas <- data.frame((ipc0 - mean(ipc0))/sd(ipc0), cycle(ipc0))
gammas <- data.frame((ipc - mean(ipc))/sd(ipc), cycle(ipc))
gammas <- data.frame(ipc0, cycle(ipc0))
initial_gamma <- tapply(gammas[,1],gammas[,2], FUN = median)



parametros <- list(
  par = data.frame(
    name = c("k1","ks","f2","df","mu[0]"),
    value = c(0.5,-1,-0.5,12, 0.68278),
    lower = c(0,-Inf,-Inf,2,-Inf),
    upper = c(1,Inf,Inf,Inf,Inf)
  ),
  mu =  0.68278,#initial_mu,
  gamma = #c(0.50767857,-0.13357143, 0.07642857,-0.09607143,-0.19357143,-0.54357143,-0.54107143,-0.46982143,-0.29142857,-0.20142857,-0.04857143)
    #c(0.43032967,-0.20752747, 0.03032967,-0.06038462,-0.13109890,-0.41538462,-0.39324176,-0.39967033,-0.35461538,-0.25538462,-0.09615385)
    c(0.37195906,-0.20488304,-0.02277778,-0.09172515,-0.22804094,-0.45383041,-0.21172515,-0.37804094,-0.39611111,-0.25611111,0.01666667)
    #as.vector(initial_gamma[1:11])
)

# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]",paste0("gamma",1:11)),
#     value = c(0.5,-1,-1.22,11,mean(ipc0),c(0.63634560, -0.11379156, 0.08584090, 0.04461767, 0.27043329, -0.11460482, -0.09033392, -0.07185322, -0.51438998, -0.29311867, -0.02459362)),
#     lower = c(0,-Inf,-Inf,2, rep(-Inf,12)),
#     upper = c(1,Inf,Inf,Inf, rep(Inf,12))
#   ),
#   mu = 0.92148341, 
#   gamma = c(0.04130164, 0.10012847, 0.01462955,-0.20743980,-0.39845356,-0.15230706, 0.67272960,-0.02919001,-0.41970979, 0.06692283, 0.27037323)
# )

parametros
bsm2 <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2", initial.optim = F)
# bsm2 <- readRDS("dados/bsm_ipc.rds")
# saveRDS(bsm2, "dados/bsm_ipc.rds")

comparar <- cbind(parametros$par, bsm2$otimizados$par)
comparar
#saveRDS(comparar, "dados/bsm_ipc0_parametros.rds")
pseudo.y <- bsm2$out[,"mu"] + bsm2$out[,"gamma"] + bsm2$out[,"u"]
# mu smooth
fk2 <- bsm(pseudo.y, beta = T, iter = 10000)
psd <- list(pseudo.y = pseudo.y, fk2 = fk2)
#saveRDS(psd, "dados/bsm_ipc_pseudo.rds")

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
plot(fk2[,"mu"])
plot(bsm2$out[,"mu"])

# diagnóstico
diag <- diag.dcs(y = bsm2$out[,"epsilon"], df = bsm2$otimizados$par[4])
hist(diag$resid.q)
diag$stats
plot(diag$resid.q,bsm2$out[,"epsilon"])

saveRDS(bsm2, "shiny_simulacao/data/bsm_ipc.rds")
saveRDS(pseudo.y, "shiny_simulacao/data/pseudoy_ipc.rds")
saveRDS(fk2, "shiny_simulacao/data/smooth_ipc.rds")
saveRDS(diag, "shiny_simulacao/data/diag_ipc.rds")
saveRDS(psd, "shiny_simulacao/data/smooth_ipc_psd.rds")

# 
library(zoo)
exportar <- data.frame(data = as.Date(cbind(bsm2$out, fk2)),cbind(bsm2$out, fk2))
colnames(exportar) <- c("data", colnames(bsm2$out), paste0(colnames(fk2),"_smoother"))
write.csv2(exportar,"shiny_simulacao/www/ultimos_resultados.csv", row.names = F)

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
dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 11),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 12)
)
# 
# y_filter <- y[1:13]
# dummy_filter <- dummy[1:13,]
# t <- 1:13
# dados <- data.frame(cbind(y_filter,t,dummy_filter))
dados <- data.frame(window(cbind(ipc,dummy),start = c(1999,1), freq = 12))
colnames(dados) <- c("y",paste0("D",1:12))
md <- lm(y ~ ., data = dados)
data.frame(as.vector(md$coefficients)[2:12])
c(0.37195906,-0.20488304,-0.02277778,-0.09172515,-0.22804094,-0.45383041,-0.21172515,-0.37804094,-0.39611111,-0.25611111,0.01666667)
c(0.43032967,-0.20752747, 0.03032967,-0.06038462,-0.13109890,-0.41538462,-0.39324176,-0.39967033,-0.35461538,-0.25538462,-0.09615385)
c(0.50767857,-0.13357143, 0.07642857,-0.09607143,-0.19357143,-0.54357143,-0.54107143,-0.46982143,-0.29142857,-0.20142857,-0.04857143)


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
