library(urca)
library(dynlm)
library(BETS)
library(seasonal)
source("./functions/core.diag.R", encoding = "utf-8")


# leitura IPC-Br e núcleos
ipc  <- window(readRDS("data/ipc.rds"), start = c(2001,1), freq = 12)
nucleo_dcs <- readRDS("data/nucleo_dcs.rds")
nucleo_tf <- readRDS("data/nucleo_tf.rds")

nucleos <- window(cbind(nucleo_tf[,3],nucleo_dcs), start = c(2001,1), end = c(2017,12), freq = 12)
colnames(nucleos) <- c("S","DCS")

# gráfico dos núcleos - lado a lado -----------------------------------------
par(mar = c(2,4,1,2), mfrow = c(1,1))

# núcleos + ipc-br
plot(ipc, main = "", lwd = 1, lty = 4, ylim = c(-0.5,3.5),
        col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"DCS"], lwd = 2, lty = 1, col = "#CD0000")
lines(nucleos[,"S"], lwd = 2, lty = 5, col = "#1874CD")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
#abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3.5, legend = c("IPC-Br","Núcleo-S","Núcleo-DCS"), lwd = c(1,2,2), lty = c(4,5,1), y.intersp = 1.5,
       col = c(1,"#1874CD","#CD0000"), cex = 1.3,bg = "white", box.col = "white",box.lwd = 0)


# núcleos sem ipc-br
plot(nucleos[,"DCS"], lwd = 2, lty = 1, col = "#CD0000", ylim = c(0.1,1.2),
     ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"S"], lwd = 2, lty = 5, col = "#1874CD")
abline(h = seq(0.2,1.2,0.2), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
#abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,1.2, legend = c("Núcleo-S","Núcleo-DCS"), lwd = c(2,2), lty = c(5,1), y.intersp = 1.5,
       col = c("#1874CD","#CD0000"), cex = 1.3,bg = "white", box.col = "white",box.lwd = 0)


# núcleo-dcs + ipc-br
plot(ipc, main = "", lwd = 1, lty = 4, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"DCS"], lwd = 2, lty = 1, col = "#CD0000")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
#abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3.5, legend = c("IPC-Br","Núcleo-DCS"), lwd = c(1,2), lty = c(4,1), y.intersp = 1.5,
       col = c(1,"#CD0000"), cex = 1.3,bg = "white", box.col = "white",box.lwd = 0)

# núcleo-s + ipc-br
plot(ipc, main = "", lwd = 1, lty = 4, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"S"], lwd = 2, lty = 5, col = "#1874CD")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
#abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3.5, legend = c("IPC-Br","Núcleo-S"), lwd = c(1,2), lty = c(4,1), y.intersp = 1.5,
       col = c(1,"#1874CD"), cex = 1.3,bg = "white", box.col = "white",box.lwd = 0)


# diagnósticos -------------------------------------

# viés
bias_full <- rbind(core.diag(ipc, nucleos[,"S"], test = "bias")$out$full,
                   core.diag(ipc, nucleos[,"DCS"], test = "bias")$out$full[2,])
bias_part1 <- rbind(core.diag(ipc, nucleos[,"S"], test = "bias")$out$part1,
                    core.diag(ipc, nucleos[,"DCS"], test = "bias")$out$part1[2,])
bias_part2 <- rbind(core.diag(ipc, nucleos[,"S"], test = "bias")$out$part2,
                    core.diag(ipc, nucleos[,"DCS"], test = "bias")$out$part2[2,])

rownames(bias_full) = rownames(bias_part1) = rownames(bias_part2) <- c("IPC-Br", "Núcleo-S", "Núcleo-DCS")

bias_full
bias_part1
bias_part2

# teste de raiz unitária
adf <- core.diag(ipc, nucleos[,"S"], test = "root", lag_level = c(12,12,12,13,12,12), lag_diff = c(11,6,10,12,11,11))
adf <- core.diag(ipc, nucleos[,"DCS"], test = "root", lag_level = c(12,12,12,0,0,0), lag_diff = c(11,6,10,0,0,0))

summary(adf$adf$full$y$level)
summary(adf$adf$part1$y$level)
summary(adf$adf$part2$y$level)
acf(adf$adf$part2$y$level@res,36)
summary(adf$adf$full$core$level)
summary(adf$adf$part1$core$level)
summary(adf$adf$part2$core$level)
acf(adf$adf$part2$core$level@res,36)
summary(adf$adf$full$y$diff)
summary(adf$adf$part1$y$diff)
summary(adf$adf$part2$y$diff)
acf(adf$adf$full$y$diff@res,36)
summary(adf$adf$full$core$diff)
summary(adf$adf$part1$core$diff)
summary(adf$adf$part2$core$diff)
acf(adf$adf$full$core$diff@res,36)


# teste de cointegração
coint <- core.diag(ipc, nucleos[,"S"], test = "coint")
coint <- core.diag(ipc, nucleos[,"DCS"], test = "coint")

# atratividade
atratividade <- core.diag(ipc, nucleos[,"S"], test = "attract", 
                          lags_y = c(1:10,12), lags_core = c(1:6))
atratividade <- core.diag(ipc, nucleos[,"DCS"], test = "attract", 
                        lags_y = c(7,12), lags_core = c(1,2))


ipc0 <- window(ipc, start = c(2009,7), freq = 12)
nucleos0 <- window(nucleos, start = c(2009,7), freq = 12)

atratividade <- core.diag(ipc0, nucleos0[,"S"], test = "attract", 
                          lags_y = c(1,2,3,4,5,6,7,8,9,10,12), lags_core = c(1,2,3,4,6,7))
atratividade <- core.diag(ipc0, nucleos0[,"DCS"], test = "attract", 
                          lags_y = c(12), lags_core = c(1,2))

# sazonalidade

a <- seas(nucleo_dcs, regression.aictest = NULL, arima.model = "(1 1 0)")
summary(a)
qs(a)

core.diag(ipc, nucleos[,"S"], test = "seas")
core.diag(ipc, nucleos[,"DCS"], test = "seas")

# previsao

core.diag(ipc, nucleos[,"S"], test = "fcst")
core.diag(ipc, nucleos[,"DCS"], test = "fcst")
core.diag(ipc0, nucleos0[,"S"], test = "fcst")
core.diag(ipc0, nucleos0[,"DCS"], test = "fcst")


