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

