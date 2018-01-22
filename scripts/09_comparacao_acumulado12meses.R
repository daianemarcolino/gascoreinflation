# pacotes, funções e leitura
library(BETS)
library(TSA)
library(MARSS)
library(zoo)
source("functions/dcs_fk_estimation.R")
source("functions/dcs_estabilidade.R")
source("functions/dcs_fcst.R")
source("functions/ERRO.R")
source("functions/bsm.R")
source("functions/diag.dcs.R", encoding = "utf8")
source("functions/ifa_functions.R", encoding = "utf8")

# leitura -----------
ipc <- readRDS("data/ipc.rds")
nucleo_tf <- readRDS("data/nucleo_tf.rds")[,3]
nucleo_dcs <- window(readRDS("data/nucleo_dcs.rds"), end = c(2017,12), freq = 12)
ipca <- readRDS("data/ipca.rds")

ipc12 <- window(acum12(ipc), start = c(2001,1), freq = 12)
ipca12 <- window(acum12(ipca), start = c(2001,1), freq = 12)

# GRÁFICO COMPARAÇÕES -------------------------------------------------
par(mar = c(2,4,1,1))#, mfrow = c(2,1))
plot(ipca12, col = "white", lty = 5,  xaxt = "n", ylab = "variação mensal percentual (%)", xlab = "", main = "", ylim = c(0,20))
axis(1, at = seq(2001,2018,0.0833*3)[-length(seq(2001,2018,0.0833*3))] , labels = substr(as.Date(ipca12),1,7)[seq(1,204,3)])
polygon(c(2000.999,2001.999,2001.999,2000.999),c(2,2,6,6), col = "#CAE1FF", border = NA)
polygon(c(2001.999,2002.999,2002.999,2001.999),c(1.5,1.5,5.5,5.5), col = "#CAE1FF", border = NA)
polygon(c(2002.999,2003.999,2003.999,2002.999),c(1.5,1.5,6.5,6.5), col = "#CAE1FF", border = NA)
polygon(c(2003.999,2004.999,2004.999,2003.999),c(3,3,8,8), col = "#CAE1FF", border = NA)
polygon(c(2004.999,2005.999,2005.999,2004.999),c(2,2,7,7), col = "#CAE1FF", border = NA)
polygon(c(2005.999,2016.999,2016.999,2005.999),c(2.5,2.5,6.5,6.5), col = "#CAE1FF", border = NA)
polygon(c(2016.999,2017.999,2017.999,2016.999),c(3,3,6,6), col = "#CAE1FF", border = NA)

lines(anual(nucleo_tf), col = "#1874CD", lty = 1, lwd = 2)
lines(anual(nucleo_dcs), col = "#CD0000", lty = 1, lwd = 3)
lines(ipca12, col = 1, lty = 4)
abline(v = seq(2001,2018,0.0833*3), h = seq(0,20,2.5),lty = 3, col = "#C9C9C9")
legend(2010,17, legend = c("IPCA 12", "Núcleo-S","Núcleo-DCS"), lwd = c(1,2,3), lty = c(4,1,1),# y.intersp = 1.5,
       col = c(1,"#1874CD","#CD0000"), cex = 1.1, bg = "white", box.col = "white", box.lwd = 0)

plot(window(ipc, start = c(2001,1),freq = 12), col = 1, lty = 5,  xaxt = "n", ylab = "variação mensal percentual (%)", xlab = "", main = "")
axis(1, at = seq(2001,2018,0.0833*4)[-52] , labels = substr(as.Date(ipca12),1,7)[seq(1,204,4)])
abline(v = seq(2001,2018,0.0833*4), h = seq(0,3,0.5),lty = 3, col = "#C9C9C9")
legend(2010,3, legend = c("IPC-Br"), lwd = c(1), lty = c(5),# y.intersp = 1.5,
       col = c(1), cex = 1.1, bg = "white", box.col = "white", box.lwd = 0)


dez_nucleo_tf <- anual(ts(nucleo_tf[cycle(nucleo_tf) == 12], start = c(2001,12), freq = 1))
dez_nucleo_dcs <- anual(ts(nucleo_dcs[cycle(nucleo_dcs) == 12], start = c(2001,12), freq = 1))
dez_ipc12 <- ts(ipc12[cycle(ipc12) == 12], start = c(2001,12), freq = 1)

ts.plot(dez_nucleo_tf, dez_nucleo_dcs, dez_ipc12, col = c(1,2,4), lwd = c(1,1,2))
