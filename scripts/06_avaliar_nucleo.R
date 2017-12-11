library(urca)
library(dynlm)

# ler núcleos
nucleos <- readRDS("dados/tendencias.rds")
ipc_ma <- readRDS("dados/ipc_medias_aparadas.rds")
ipc_ma <- readRDS("dados/ipc_medias_aparadas_sa.rds")
ipc_ma <- readRDS("dados/ipc_medias_aparadas_sa_mm3.rds")
nucleo_restricao <- readRDS("dados/nucleo_restricaok1df.rds")
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)

y <- ipc
core <- nucleo_restricao$out[,"mu"]#nucleos[,"dcs_padrao_psi_dummy"]
conf <- 0.95
# 
# nucleo <- ((nucleo/100+1)^12-1)*100
# ipc <- ((ipc/100+1)^12-1)*100
# y=ipc
# core=nucleos$core_res
# 
# y=window(ipc, start = c(2008,5), freq = 12)
# core=window(nucleo, start = c(2008,5), freq = 12)
# 
# ts.plot(y,core)
# write.csv2(data.frame(data = as.Date(core), round(((core/100+1)^12-1)*100,2)), "tendencia_final.csv", row.names = F)

# NOVO ---------------------------------------
source("functions/dcs_fk_estimation.R")

initial_gamma <- c(0.48070729,-0.05859547, 0.15370156, 0.10032539,-0.03874261,-0.26932386,-0.12546543,-0.18050091,-0.20419881,-0.07149072, 0.06560687)

parametros_psi_dummy1 <- list(
  par = data.frame(
    name = c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11), "d1", "d2"),
    value = c(0.1,0.1,-1,12, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11], 1,1),
    lower = c(0,0,-Inf,4,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11), -Inf, -Inf),
    upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11), Inf, Inf)
  ),
  gamma = NA,
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2000,7)),
                d2 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)

out <- dcs_padrao_psi_dummy <- dcs_fk_estimation(ipc, initial = parametros_psi_dummy1, type = "BSM2_beta_psi", outlier = T, otimo = T)

y = ipc
core = out$out[,"mu"]
conf = 0.95

# não viesado
# não estacionário
# cointegrado com a inflação

# estabilidade do núcleo -----------------------------
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
source('functions/dcs_estabilidade.R', encoding = 'UTF-8')
source("functions/dcs_fk_estimation.R")
initial_gamma <- c(0.48070729,-0.05859547, 0.15370156, 0.10032539,-0.03874261,-0.26932386,-0.12546543,-0.18050091,-0.20419881,-0.07149072, 0.06560687)

parametros_psi_dummy1 <- list(
  par = data.frame(
    name = c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11), "d1", "d2"),
    value = c(0.1,0.1,-1,12, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11], 1,1),
    lower = c(0,0,-Inf,4,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11), -Inf, -Inf),
    upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11), Inf, Inf)
  ),
  gamma = NA,
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2000,7)),
                d2 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)

out <- dcs_estabilidade(ipc, start = c(2013,9), initial = parametros_psi_dummy1, type = "BSM2_beta_psi", outlier = T)
saveRDS(out, "estabilidade.rds")
ts.plot(round(out$out,2), col = 1:2, lty = c(1,2), type = "o")
abline(h = seq(0.3,0.8, 0.1), v = 2014:2017, col = "lightgrey", lty = 3)
legend("topleft", legend = c("mês a mês","completo"), col = 1:2, lty = c(1,3), bty = "n",pch = 1)

ts.plot(ipc, core, lwd = 1:2)
y = ipc
core = out$out[,"mu"]
conf = 0.95

