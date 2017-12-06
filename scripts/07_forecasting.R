# leitura
setwd("C:/Users/master/Desktop/DAIANE")
ipc <- ts(read.csv2("VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
library(BETS)
library(zoo)
source('dcs_fcst.R', encoding = 'UTF-8')
source('dcs_fk_estimation.R', encoding = 'UTF-8')

initial_gamma <- c(0.48070729,-0.05859547, 0.15370156, 0.10032539,-0.03874261,-0.26932386,-0.12546543,-0.18050091,-0.20419881,-0.07149072, 0.06560687)

parametros_iniciais <- list(
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

datas <- list(ano = as.numeric(substr(as.Date(ipc),1,4))[(length(ipc)-48):length(ipc)],
              mes = as.numeric(substr(as.Date(ipc),6,7))[(length(ipc)-48):length(ipc)])

iter <- list()
for(i in 1:length(datas$ano)){
  iter[[i]] <- dcs_fcst(y = ipc, start = c(datas$ano[i], datas$mes[i]), initial = parametros_iniciais, type = "BSM2_beta_psi", outlier = T, m = 1000)
}
saveRDS(iter, "iter.rds")
ts.plot(round(na.omit(cbind(ipc, iter1$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
  