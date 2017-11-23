# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
library(MARSS)
source("functions/dcs_fk_estimation.R")
source("functions/dcs_fk_estimation_exercise.R")
source("functions/ERRO.R")
source("functions/bsm.R")
source("functions/diag.dcs.R", encoding = "utf8")

# leitura -----------

turistas <- log(ts(read.csv2("dados/turistas.csv")[,2], start = c(2000,1), freq = 12))
plot(turistas)
#saveRDS(turistas, "shiny_simulacao/data/turistas.rds")

# replicar artigo ------------------------------------------------- 
bsm0 <- dcs_fk_estimation(turistas, type = "BSM_artigo")

pseudo.y <- (1 - bsm0$out[,"b"])*turistas + bsm0$out[,"b"]*bsm0$out[,"mu"]
ts.plot(turistas, pseudo.y, col = 1:2)

# mu smooth
fk1 <- bsm(pseudo.y,type = "BSM2", iter = 9)

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