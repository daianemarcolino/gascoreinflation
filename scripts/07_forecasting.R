# leitura
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
library(BETS)
library(zoo)
library(plotrix)
source('functions/dcs_fcst.R', encoding = 'UTF-8')
source('functions/dcs_fk_estimation.R', encoding = 'UTF-8')

initial_gamma <- c(0.48070729,-0.05859547, 0.15370156, 0.10032539,-0.03874261,-0.26932386,-0.12546543,-0.18050091,-0.20419881,-0.07149072, 0.06560687)

parametros <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","df","beta","mu[1|0]",paste0("gamma",1:11)          ,"psi","phi","k3"),
    value = c(0.1 ,0   ,0.5 ,5   ,6   ,0     ,0        ,as.vector(initial_gamma)[1:11],1    ,0.1  ,0.5),
    lower = c(0.05,0   ,0.15 ,-Inf,4   ,0     ,-Inf     ,rep(-Inf,11)                  ,-Inf ,-1   ,0),
    upper = c(Inf ,0   ,Inf ,Inf ,Inf ,0     ,Inf      ,rep(Inf,11)                   ,Inf  ,1    ,Inf)
  )
)
parametros

datas <- list(ano = as.numeric(substr(as.Date(ipc),1,4))[(length(ipc)-47):length(ipc)],
              mes = as.numeric(substr(as.Date(ipc),6,7))[(length(ipc)-47):length(ipc)])

# iter <- list()
# for(i in 1:length(datas$ano)){
#   iter[[i]] <- dcs_fcst(y = ipc, start = c(datas$ano[i], datas$mes[i]), initial = parametros, type = "BSM3", outlier = F, m = 3)
#   message("\n \n \n \n \n \n \n ================== date:", datas$ano[i],"/",datas$mes[i])
# }
# saveRDS(iter, "iter_BSM3_2.rds")
ts.plot(round(na.omit(cbind(ipc, iter[[1]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
ts.plot(round(na.omit(cbind(ipc, iter[[2]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))
ts.plot(round(na.omit(cbind(ipc, iter[[3]]$fcst)),2), col = c(1,2,4,4), lty = c(3,2,1,1))

iter <- readRDS("dados/iter_BSM3.rds") 

# gráfico
for(i in 1:49){
  ts.plot((cbind(ipc,iter[[i]]$fcst[,1])), col = 1:2)
  Sys.sleep(0.5)
}

fcst_mean <- sapply(iter, FUN = function(x)  x$fcst[,1])
fcst_lower <- sapply(iter, FUN = function(x)  x$fcst[,2])
fcst_upper <- sapply(iter, FUN = function(x)  x$fcst[,3])

fcsts_mean  <- do.call(cbind, fcst_mean)
fcsts_lower <- do.call(cbind, fcst_lower)
fcsts_upper <- do.call(cbind, fcst_upper)

colnames(fcsts_mean) = colnames(fcsts_lower) = colnames(fcsts_upper) <- paste0("fcst", 1:ncol(fcsts_mean))

ts.plot(ipc, fcsts_lower[,1], fcsts_mean[,1], fcsts_upper[,1], col = c(1,4,2,4))
ts.plot(ipc, fcsts_lower[,36], fcsts_mean[,36], fcsts_upper[,36], col = c(1,4,2,4))
ts.plot(ipc, fcsts_lower[,40], fcsts_mean[,40], fcsts_upper[,40], col = c(1,4,2,4))

# previsoes 1 passo a frente

View(fcsts_mean)

k <- data.frame(fcsts_mean)
fcsts_df <- apply(k, MARGIN = 2, FUN = function(x) c(na.omit(x),rep(NA,sum(is.na(x)))))

rmse <- matrix(NA,nrow = 49, ncol = 1)

for(i in 1:nrow(fcsts_df)){
  ychapeu <- ts(c(fcsts_df[i,]), start = c(datas$ano[i], datas$mes[i]), freq = 12)
  d <- cbind(ipc,ychapeu)
  rmse[i,1] <- sqrt(mean((d[,1] - d[,2])^2, na.rm = T))
}

# erro linha i passos a frente
barp(rmse[c(1,3,6,9,12),], names.arg = c(1,3,6,9,12))
barp(rmse)
