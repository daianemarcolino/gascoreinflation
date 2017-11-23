rm(list = ls())
# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
library(MARSS)
library(zoo)
source("functions/dcs_fk_estimation.R")
source("functions/dcs_fk_estimation_exercise.R")
source("functions/ERRO.R")
source("functions/bsm.R")
source("functions/diag.dcs.R", encoding = "utf8")

# leitura -----------
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
plot(ipc)

# filtro ipc
ipc0 <- window(ipc,start = c(2004,1),freq = 12)
# initial_mu <- median(ipc)
# gammas <- data.frame((ipc0 - mean(ipc0))/sd(ipc0), cycle(ipc0))
# gammas <- data.frame((ipc - mean(ipc))/sd(ipc), cycle(ipc))
# gammas <- data.frame(ipc0, cycle(ipc0))
# initial_gamma <- tapply(gammas[,1],gammas[,2], FUN = mean)
# 
# y = ipc0
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
# dados <- cbind(y, dummy)
# colnames(dados) <- c("y",paste0("D",1:12))
# md <- lm(y ~ ., data = dados)
# data.frame(as.vector(md$coefficients)[2:12])
#initial_gamma <- c(0.37195906,-0.20488304,-0.02277778,-0.09172515,-0.22804094,-0.45383041,-0.21172515,-0.37804094,-0.39611111,-0.25611111, 0.01666667)
#initial_gamma <- c(0.50767857,-0.13357143, 0.07642857,-0.09607143,-0.19357143,-0.54357143,-0.54107143,-0.46982143,-0.29142857,-0.20142857,-0.04857143)
initial_gamma <- c(0.05516970,0.09563170,0.05699505,-0.18377041,-0.43230580,-0.17613783,0.69144612,-0.04360480,-0.44064156,0.07471822,0.27193868)
initial_gamma <- c(0.153357232, 0.017716173, 0.132833222,-0.139159857,-0.404060430,-0.145550980, 0.660517961,-0.048933837,-0.449613456, 0.022034873, 0.223473635)
#initial_gamma <- c(0.512692417,-0.041189734, 0.144115903, 0.088499102,-0.060261330,-0.266698219,-0.138449064,-0.195460536,-0.221822482,-0.069103134, 0.092489454)
initial_gamma <- c(0.5193,-0.0469, 0.1511, 0.0878,-0.0606,-0.2625,-0.1520,-0.1967,-0.2275,-0.0687, 0.0995)
initial_gamma <- c(0.48070729,-0.05859547, 0.15370156, 0.10032539,-0.03874261,-0.26932386,-0.12546543,-0.18050091,-0.20419881,-0.07149072, 0.06560687)


# BSM NOVO --------------------------------------------------------

# RESTRITO
parametros <- list(
  par = data.frame(
    name = c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11)),
    value = c(0.1,0.1,-1,12, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11]),
    lower = c(0,0,-Inf,4,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11)),
    upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11))
  ),
  gamma = NA
)
parametros
dcs_res <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2_beta_psi")
ts.plot(ipc,dcs_res$out[,"mu"], col = 1:2)
ts.plot(dcs_res$out[,"mu"], col = 1:2)
ts.plot(dcs_res$out[,"gamma"], col = 1:2)
pseudo.y <- (1 - dcs_res$out[,"b"])*ipc + dcs_res$out[,"b"]*(dcs_res$out[,"mu"] + dcs_res$out[,"gamma"])

x <- data.frame(parametros$par, otimo = round(dcs_res$otimizados$par,4))
#x[2,3:5] <- NA 
x
# mu smooth
bsm_res <- bsm(pseudo.y, type = "BSM2_beta", iter = 10000)
ts.plot(bsm_res[,"mu"],ipc, col = 2:1)

ts.plot(bsm_res[,"mu"], dcs_res$out[,"mu"], col = 2:1)

saveRDS(list(core_res = dcs_res$out[,"mu"], core_ress = bsm_res[,"mu"]), "./dados/nucleo.rds")


# fig 6
par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(ipc,dcs_res$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
ts.plot(dcs_res$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
ts.plot(dcs_res$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
abline(h=0, lty = 3, col = 2)
ts.plot(dcs_res$out[,c("u")],dcs_res$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")

# fig 7
par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(dcs_res$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
ts.plot(dcs_res$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
acf(dcs_res$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
acf(dcs_res$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")

# fig 8
# ts.plot(fk2[,"mu"],bsm2$out[,"mu"], col = 1, lty = c(1,3))
# ts.plot(fk2_beta[,"mu"],bsm2_beta$out[,"mu"], col = 1, lty = c(1,3))
# ts.plot(fk2[,"mu"],fk2_beta[,"mu"], col = 1, lty = c(1,3))
# ts.plot(fk2[,"mu"],ipc, col = 1, lty = c(1,3))
# ts.plot(fk2_beta[,"mu"],ipc, col = 1, lty = c(1,3), lwd = c(2,1))

write.csv2(data.frame(data = as.Date(bsm2_beta$out[,"mu"]), round(((bsm2_beta$out[,"mu"]/100+1)^12-1)*100,2)), "tendencia5.csv", row.names = F)

# diagnóstico
diag <- diag.dcs(y = dcs_res$out[,"epsilon"], df = dcs_res$otimizados$par[4])

hist(diag$resid.q, main = "residuo quantilico")
hist(dcs_res$out[,"epsilon"], main = "residuo de pearson")

# BSM para IPC ----------------------------------------------------
# valores iniciais
parametros <- list(
  par = data.frame(
    name = c("k1","ks","f2","df","mu[0]","beta"),
    value = c(0.02,0.5,0.5,10, 0.49,-0.01),
    lower = c(0,0.001,-Inf,4,-Inf,-Inf),
    upper = c(1,Inf,Inf,Inf,Inf,Inf)
  ),
  gamma = as.vector(initial_gamma)[1:11]
)

# IRRESTRITO
parametros <- list(
  par = data.frame(
    name = c("k1","ks","f2","df","mu[0]","beta", paste0("gamma",1:11)),
    value = c(0.1,0.1,-1,12, 0.5,-0.1, as.vector(initial_gamma)[1:11]),
    lower = c(0,0,-Inf,4,-Inf,-Inf, rep(-Inf,11)),
    upper = c(1,Inf,Inf,Inf,Inf,Inf,rep(Inf,11))
  ),
  gamma = NA
)
parametros
dcs_full <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2_beta")
ts.plot(ipc,dcs_full$out[,"mu"], col = 1:2)

pseudo.y <- (1 - dcs_full$out[,"b"])*ipc + dcs_full$out[,"b"]*(dcs_full$out[,"mu"] + dcs_full$out[,"gamma"])

# mu smooth
bsm_full <- bsm(pseudo.y, type = "BSM2_beta", iter = 10000)
ts.plot(bsm_full[,"mu"],ipc, col = 2:1)

ts.plot(bsm_full[,"mu"], dcs_full$out[,"mu"], col = 2:1)

# RESTRITO
parametros <- list(
  par = data.frame(
    name = c("k1","ks","f2","df","mu[0]","beta", paste0("gamma",1:11)),
    value = c(0.1,0.1,-1,12, 0.5,-0.1, as.vector(initial_gamma)[1:11]),
    lower = c(0,0,-Inf,4,-Inf,-Inf, rep(-Inf,11)),
    upper = c(0.1,Inf,Inf,Inf,Inf,Inf,rep(Inf,11))
  ),
  gamma = NA
)
parametros
dcs_res <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2_beta")
ts.plot(ipc,dcs_res$out[,"mu"], col = 1:2)
ts.plot(dcs_res$out[,"gamma"], col = 1:2)
pseudo.y <- (1 - dcs_res$out[,"b"])*ipc + dcs_res$out[,"b"]*(dcs_res$out[,"mu"] + dcs_res$out[,"gamma"])

x <- data.frame(parametros$par, otimo = round(dcs_res$otimizados$par,4))
x[2,3:5] <- NA 
x
# mu smooth
bsm_res <- bsm(pseudo.y, type = "BSM2_beta", iter = 10000)
ts.plot(bsm_res[,"mu"],ipc, col = 2:1)

ts.plot(bsm_res[,"mu"], dcs_res$out[,"mu"], col = 2:1)

saveRDS(list(core_res = dcs_res$out[,"mu"], core_ress = bsm_res[,"mu"]), "./dados/nucleo.rds")


# fig 6
par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(ipc,dcs_res$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
ts.plot(dcs_res$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
ts.plot(dcs_res$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
abline(h=0, lty = 3, col = 2)
ts.plot(dcs_res$out[,c("u")],dcs_res$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")

# fig 7
par(mfrow = c(2,2), mar = c(3,3,2,3))
ts.plot(dcs_res$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
ts.plot(dcs_res$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
acf(dcs_res$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
acf(dcs_res$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")

# fig 8
# ts.plot(fk2[,"mu"],bsm2$out[,"mu"], col = 1, lty = c(1,3))
# ts.plot(fk2_beta[,"mu"],bsm2_beta$out[,"mu"], col = 1, lty = c(1,3))
# ts.plot(fk2[,"mu"],fk2_beta[,"mu"], col = 1, lty = c(1,3))
# ts.plot(fk2[,"mu"],ipc, col = 1, lty = c(1,3))
# ts.plot(fk2_beta[,"mu"],ipc, col = 1, lty = c(1,3), lwd = c(2,1))

write.csv2(data.frame(data = as.Date(bsm2_beta$out[,"mu"]), round(((bsm2_beta$out[,"mu"]/100+1)^12-1)*100,2)), "tendencia5.csv", row.names = F)

# diagnóstico
diag <- diag.dcs(y = dcs_res$out[,"epsilon"], df = dcs_res$otimizados$par[4])

hist(diag$resid.q, main = "modelo sem beta")

diag$stats

saveRDS(bsm2, "shiny_simulacao/data/bsm_ipc.rds")
saveRDS(pseudo.y, "shiny_simulacao/data/pseudoy_ipc.rds")
saveRDS(fk2, "shiny_simulacao/data/smooth_ipc.rds")
saveRDS(diag, "shiny_simulacao/data/diag_ipc.rds")
#saveRDS(psd, "shiny_simulacao/data/smooth_ipc_psd.rds")

# nucleo final
nucleo <- bsm2$out[,"mu"]
nucleo <- fk2[,"mu"]
# exercício

out <- NULL
vero <- NULL
k1 <- NULL
ks <- seq(0,1,by = 0.01)
mu <- NULL
beta <- NULL
dfs <- 4:30
for(i in 1:length(ks)){
  #out[[i]] <- dcs_fk_estimation_exercise(ipc, initial = parametros, type = "BSM2_beta", k1 = k1[i])
  #out[[i]] <- dcs_fk_estimation_exercise(ipc, initial = parametros, type = "BSM2_beta", df = dfs[i])
  out[[i]] <- dcs_fk_estimation_exercise(ipc, initial = parametros, type = "BSM2_beta", ks = ks[i])
  vero[i] <- out[[i]]$loglik
  
}
k <- cbind(ks = ks, loglik = -vero)
plot(k[,1:2], xlab = "graus de liberdade", ylab = "loglik")
round(k,3)

library(zoo)
exportar <- data.frame(data = as.Date(cbind(bsm2$out, fk2)),cbind(bsm2$out, fk2))
colnames(exportar) <- c("data", colnames(bsm2$out), paste0(colnames(fk2),"_smoother"))
write.csv2(exportar,"shiny_simulacao/www/ultimos_resultados.csv", row.names = F)
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

# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]",paste0("gamma",1:11)),
#     value = c(0.8,-1.1,-1.6,6,mean(ipc),initial_gamma[1:11]),
#     lower = c(0,-Inf,-Inf,2, rep(-Inf,12)),
#     upper = c(1,Inf,Inf,Inf, rep(Inf,12))
#   ),
#   mu = initial_mu,
#   gamma = as.vector(initial_gamma[1:11])
# )
# 
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
y = ipc0
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
dados <- data.frame(window(cbind(ipc0,dummy),start = start(ipc0), freq = 12))
colnames(dados) <- c("y",paste0("D",1:12))
md <- lm(y ~ ., data = dados)
data.frame(as.vector(md$coefficients)[2:12])
c(0.37195906,-0.20488304,-0.02277778,-0.09172515,-0.22804094,-0.45383041,-0.21172515,-0.37804094,-0.39611111,-0.25611111,0.01666667)
c(0.43032967,-0.20752747, 0.03032967,-0.06038462,-0.13109890,-0.41538462,-0.39324176,-0.39967033,-0.35461538,-0.25538462,-0.09615385)
c(0.50767857,-0.13357143, 0.07642857,-0.09607143,-0.19357143,-0.54357143,-0.54107143,-0.46982143,-0.29142857,-0.20142857,-0.04857143)
c(0.42891026,-0.19647436, 0.04583333,-0.04032051,-0.14724359,-0.45878205,-0.42032051,-0.44262821,-0.33250000,-0.23250000,-0.08250000)
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
