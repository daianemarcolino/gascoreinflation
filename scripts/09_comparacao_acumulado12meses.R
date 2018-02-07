# pacotes, funções e leitura
library(BETS)
library(TSA)
library(zoo)
library(lmtest)
library(dynlm)
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
nucleo_comb <- readRDS("data/nucleo_comb.rds")
ipca <- readRDS("data/ipca.rds")

ipc12 <- window(acum12(ipc), start = c(2001,1), freq = 12)
ipca12 <- window(acum12(ipca), start = c(2001,1), freq = 12)
selic <- window( readRDS("data/selic.rds"), start = c(2001,1), end = c(2017,12), freq = 12)
expec <- ts(read.csv2("data/expectativa_ipca_focus.csv")[,3], start = c(2002,1), freq = 12)

# GRÁFICO COMPARAÇÕES -------------------------------------------------
par(mar = c(2,4,1,4))#, mfrow = c(2,1))
plot(ipca12, col = "white", lty = 5,  xaxt = "n", ylab = "variação mensal percentual anualizada (%)", xlab = "", main = "", ylim = c(0,20))
axis(1, at = seq(2001,2018,0.0833*12)[-length(seq(2001,2018,0.0833*12))] , labels = substr(as.Date(ipca12),1,4)[seq(1,204,12)])
polygon(c(2000.999,2001.999,2001.999,2000.999),c(2,2,6,6), col = "#CAE1FF", border = NA)
polygon(c(2001.999,2002.999,2002.999,2001.999),c(1.5,1.5,5.5,5.5), col = "#CAE1FF", border = NA)
polygon(c(2002.999,2003.999,2003.999,2002.999),c(1.5,1.5,6.5,6.5), col = "#CAE1FF", border = NA)
polygon(c(2003.999,2004.999,2004.999,2003.999),c(3,3,8,8), col = "#CAE1FF", border = NA)
polygon(c(2004.999,2005.999,2005.999,2004.999),c(2,2,7,7), col = "#CAE1FF", border = NA)
polygon(c(2005.999,2016.999,2016.999,2005.999),c(2.5,2.5,6.5,6.5), col = "#CAE1FF", border = NA)
polygon(c(2016.999,2017.999,2017.999,2016.999),c(3,3,6,6), col = "#CAE1FF", border = NA)
#lines(anual(nucleo_tf), col = "#1874CD", lty = 1, lwd = 2)
lines(anual(nucleo_dcs), col = "#CD0000", lty = 1, lwd = 3)
lines(ipca12, col = 1, lty = 1)
abline(v = seq(2000,2019,0.0833*6), h = seq(0,30,2.5),lty = 3, col = "#C9C9C9")
par(new = TRUE)
plot(selic, col = "#78AB46", lty = 5, lwd = 2, yaxt = "n", xaxt = "n", ylab = "")
axis(4, col = "#78AB46", col.ticks = "#78AB46", font = 2, col.axis = "#78AB46" )#, at = seq(0,25,5) , labels = seq(0,25,5))
mtext("Selic", side = 4, line = 2.5, font = 2, col = "#78AB46")
mtext("IPCA 12", side = 2, line = 2, font = 1, col = "#000000", adj = 0.3)
mtext("&", side = 2, line = 2, font = 1, col = "#000000", adj = 0.47)
mtext("Núcleo-DCS", side = 2, line = 2, font = 2, col = "#CD0000", adj = 0.7)
legend(2010,25, legend = c("IPCA 12","Selic","Núcleo-DCS"), lwd = c(1,2,3), lty = c(1,5,1),# y.intersp = 1.5,
       col = c(1,"#78AB46","#CD0000"), cex = 1.1, bg = "white", box.col = "white", box.lwd = 0)

plot(window(ipc, start = c(2001,1),freq = 12), col = 1, lty = 5,  xaxt = "n", ylab = "variação mensal percentual (%)", xlab = "", main = "")
axis(1, at = seq(2001,2018,0.0833*4)[-52] , labels = substr(as.Date(ipca12),1,7)[seq(1,204,4)])
abline(v = seq(2001,2018,0.0833*4), h = seq(0,3,0.5),lty = 3, col = "#C9C9C9")
legend(2010,3, legend = c("IPC-Br"), lwd = c(1), lty = c(5),# y.intersp = 1.5,
       col = c(1), cex = 1.1, bg = "white", box.col = "white", box.lwd = 0)


dez_nucleo_tf <- anual(ts(nucleo_tf[cycle(nucleo_tf) == 12], start = c(2001,12), freq = 1))
dez_nucleo_dcs <- anual(ts(nucleo_dcs[cycle(nucleo_dcs) == 12], start = c(2001,12), freq = 1))
dez_ipc12 <- ts(ipc12[cycle(ipc12) == 12], start = c(2001,12), freq = 1)

ts.plot(dez_nucleo_tf, dez_nucleo_dcs, dez_ipc12, col = c(1,2,4), lwd = c(1,1,2))

# expectativa
ts.plot(expec,ipca12)
x <- window(cbind(ipca12,expec,nucleo_dcs = anual(nucleo_dcs)), start = c(2014,1), freq = 12)
par(mar = c(2,4,1,1))#, mfrow = c(2,1))
plot(x[,1], col = "white", lty = 5,  xaxt = "n", ylab = "variação mensal percentual (%)", xlab = "", main = "", ylim = c(2,12.5))
axis(1, at = seq(2001,2018,0.0833*12)[-length(seq(2001,2018,0.0833*12))] , labels = substr(as.Date(ipca12),1,4)[seq(1,204,12)])
polygon(c(2013.999,2016.999,2016.999,2013.999),c(2.5,2.5,6.5,6.5), col = "#CAE1FF", border = NA)
polygon(c(2016.999,2017.999,2017.999,2016.999),c(3,3,6,6), col = "#CAE1FF", border = NA)

#lines(anual(nucleo_tf), col = "#1874CD", lty = 1, lwd = 2)
lines(x[,"expec"], col = "#F0A804", lty = 6, lwd = 3)
lines(x[,"nucleo_dcs"], col = "#CD0000", lty = 1, lwd = 3)
lines(x[,"ipca12"], col = 1, lty = 1)
abline(v = seq(2000,2019,0.0833*6), h = seq(0,12,2),lty = 3, col = "#C9C9C9")
legend(2014,12, legend = c("IPCA 12","Expectativa","Núcleo-DCS"), lwd = c(1,3,3), lty = c(1,6,1),# y.intersp = 1.5,
       col = c(1,"#F0A804","#CD0000"), cex = 1.1, bg = "white", box.col = "white", box.lwd = 0)

# GRANGER CAUSA -----------------------------------------

grangertest(nucleo_comb ~ ipca12, order = 6)
grangertest(ipca12 ~ nucleo_comb, order = 6)

grangertest(nucleo_dcs ~ ipca12, order = 6)
grangertest(ipca12 ~ nucleo_dcs, order = 6)

grangertest(nucleo_tf ~ ipca12, order = 6)
grangertest(ipca12 ~ nucleo_tf, order = 6)

grangertest(nucleo_tf ~ nucleo_dcs, order = 6)
grangertest(nucleo_dcs ~ nucleo_tf, order = 6)

grangertest(nucleo_tf ~ nucleo_comb, order = 1)
grangertest(nucleo_dcs ~ nucleo_comb, order = 1)

grangertest(nucleo_comb ~ nucleo_dcs, order = 12)
grangertest(nucleo_comb ~ nucleo_tf, order = 12)

m <- dynlm(nucleo_dcs ~ L(nucleo_tf,c(1,2,4,5)) +  L(nucleo_dcs,c(1,6)))
summary(m)

m <- dynlm(nucleo_tf ~ L(nucleo_dcs,c(5,6)) +  L(nucleo_tf,c(1:7)))
summary(m)

m <- dynlm(nucleo_dcs ~ L(ipca12,c(0,1,12)) +  L(nucleo_dcs,c(1,12)))
summary(m)
m <- dynlm(ipca12 ~ L(ipca12,c(1,12)) +  L(nucleo_dcs,c(0,12)))
summary(m)


m <- dynlm(nucleo_tf ~ L(ipca12,c(0,1,12)) +  L(nucleo_tf,c(1,12)))
summary(m)
m <- dynlm(ipca12 ~ L(ipca12,c(1)) +  L(nucleo_tf,c(0,1,12)))
summary(m)

acf(m$residuals,48)


