# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
library(MARSS)
library(zoo)
library(seasonal)
source("functions/dcs_fk_estimation.R")
source("functions/dcs_estabilidade.R")
source("functions/dcs_fcst.R")
source("functions/ERRO.R")
source("functions/bsm.R")
source("functions/diag.dcs.R", encoding = "utf8")
source("functions/ifa_functions.R", encoding = "utf8")

# leitura -----------
ipc <- window(readRDS("data/ipc.rds"), start = c(2001,1), freq = 12)
nucleo_tf <- readRDS("data/nucleo_tf.rds")

# NÚCLEO-DCS -----------------------------------------------
initial_gamma <- c(0.418235294,-0.215294118,-0.011764706,-0.048235294,-0.183529412,-0.445294118,-0.320000000,-0.401764706,-0.374705882,-0.247647059, 0.008823529)

parametros3_normald <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","df","beta","mu[1|0]",paste0("gamma",1:11)          ,"psi","phi","k3","d1"),
    value = c(0.1 ,0   ,0.5 ,5   ,0   ,0     ,0        ,as.vector(initial_gamma)[1:11],1    ,0.1  ,0.5 ,0   ),
    lower = c(0.0 ,0   ,0.0 ,-Inf,0   ,0     ,-Inf     ,rep(-Inf,11)                  ,-Inf ,-1   ,0   ,-Inf),
    upper = c(Inf ,0   ,Inf ,Inf ,0   ,0     ,Inf      ,rep(Inf,11)                   ,Inf  ,1    ,Inf ,Inf )
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros3_normald

# estabilidade
# estabilidade <- dcs_estabilidade(y = ipc, start = c(2016,1), initial = parametros3_normald, type = "BSM3_normal", outlier = T)
# saveRDS(estabilidade, "data/estabilidade.rds")
estabilidade <- readRDS("data/estabilidade.rds")
ts.plot(round(estabilidade$out,2), col = 1:2, lty = c(6,3), ylim = c(0.3,0.9))

# NÚCLEO-S -------------------------------------------------

# filtro1 <- nucleo_tf[,1]
# 
# ajustar_filtro1 <- function(x, end = NULL){
#   
#   y <- window(x, end = end, freq = 12)
#   mes <- casefold(month.abb,upper = F)[as.numeric(substr(as.Date(tail(y,1)),6,7))]
#   data_fim <- paste0("2009.jan,",substr(as.Date(tail(y,1)),1,4),".",mes)
#   
#   ajuste <- seas(y, regression.aictest = NULL,
#                  transform.function = "none", 
#                  arima.model = "(1 1 0)(1 0 0)",
#                  series.modelspan = data_fim)
#   
#   ajuste$series$s11
#   
# }
# 
# filtro2_inicio <- matrix(NA, ncol = 24, nrow = length(na.omit(filtro1)))
# 
# ajuste <- NULL
# k <- 0
# for(ano in 2016:2017){
#   for(mes in 1:12){
#     ajuste <- ajustar_filtro1(filtro1, end = c(ano,mes))
#     filtro2_inicio[,mes + k] <- c(ajuste,rep(NA,nrow(filtro2_inicio) - length(ajuste)))
#     message(paste0("- data: ", mes,"/",ano))
#   }
#   k <- 12
# }
# 
# filtro2 <- ts(filtro2_inicio, start = start(filtro1), freq = 12)
# 
# filtro3 <- filtro2*NA
# 
# # guardar os resultados com o terceiro filtro mm3
# for(i in 1:ncol(filtro2)){
#   filtro3[,i] <- round(geom3(filtro2[,i]),2)
# }
# 
# filtro3
# 
# # criar série do núcleo guardando o último ajuste de cada mês
# nucleo_tf_full <- round(window(nucleo_tf[,3], start = c(2016,1), freq = 12),2)
# nucleo_tf_month <- nucleo_tf_full*NA 
# 
# for(i in 1:ncol(filtro3)){
#   nucleo_tf_month[i] <- tail(na.omit(filtro3[,i]),1)
# }
#
# estabilidade_tf <- list(out = cbind(core_month = nucleo_tf_month, core_full = nucleo_tf_full))
# saveRDS(estabilidade_tf, "data/estabilidade_tf.rds")
estabilidade_tf <- readRDS("data/estabilidade_tf.rds")
ts.plot(round(estabilidade_tf$out,2), col = 1:2, lty = c(6,3), ylim = c(0.3,0.9))

# GRÁFICO COMPARAÇÕES -------------------------------------------------
par(mar = c(2,4,2,1), mfrow = c(1,2))
plot(estabilidade_tf$out[,1], col = "#1874CD", lty = 5, ylim = c(0.25,0.9), xaxt = "n", ylab = "variação mensal percentual (%)", xlab = "", main = "Núcleo-S")
lines(estabilidade_tf$out[,2], col = "#1874CD", lty = 1, lwd = 2)
abline(v = seq(2016,2018,0.0833), h = seq(0.3,0.9,0.1),lty = 3, col = "#C9C9C9")
axis(1, at = seq(2016,2018,0.0833*3) , labels = substr(as.Date(estabilidade_tf$out),1,7)[seq(1,25,3)])
legend(2017,0.8, legend = c("mês a mês","completo"), lwd = c(1,2), lty = c(5,1),# y.intersp = 1.5,
       col = c("#1874CD","#1874CD"), cex = 1.1,bg = "white", box.col = "white",box.lwd = 0)

plot(estabilidade$out[,1], col = "#CD0000", lty = 5, ylim = c(0.25,0.9), xaxt = "n", ylab = "variação mensal percentual (%)", xlab = "", main = "Núcleo-DCS")
lines(estabilidade$out[,2], col = "#CD0000", lty = 1, lwd = 2)
abline(v = seq(2016,2018,0.0833), h = seq(0.3,0.9,0.1),lty = 3, col = "#C9C9C9")
axis(1, at = seq(2016,2018,0.0833*3) , labels = substr(as.Date(estabilidade_tf$out),1,7)[seq(1,25,3)])
legend(2017,0.8, legend = c("mês a mês","completo"), lwd = c(1,2), lty = c(5,1),# y.intersp = 1.5,
       col = c("#CD0000","#CD0000"), cex = 1.1,bg = "white", box.col = "white",box.lwd = 0)



