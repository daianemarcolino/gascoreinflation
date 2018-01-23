library(urca)
library(dynlm)
library(BETS)
library(seasonal)
source("./functions/core.diag.R", encoding = "utf-8")
source("./functions/core.fcst.R", encoding = "utf-8")

# baixar núcleos
# nucleo_fgv <- BETS.get(4467)
# nucleo_ipca_mas <- BETS.get(4466) # Núcleo médias aparadas com suavização
# nucleo_ipca_ma <- BETS.get(11426) # Núcleo médias aparadas sem suavização
# nucleo_ipca_ex1 <- BETS.get(11427) # Núcleo por exclusão - Sem monitorados e alimentos no domicílio 
# nucleo_ipca_ex2 <- BETS.get(16121) # Núcleo por exclusão - ex2
# nucleo_ipca_dp <- BETS.get(16122) # Núcleo de dupla ponderação
# nucleos <- window(cbind(nucleo_fgv, nucleo_ipca_mas, nucleo_ipca_ma,nucleo_ipca_ex1,nucleo_ipca_ex2, nucleo_ipca_dp),
#                   start = c(2001,1), freq = 12)
# colnames(nucleos) <- c("FGV","IPCA_MAS","IPCA_MA","IPCA_EX1","IPCA_EX2","IPCA_DP")
# saveRDS(nucleos, "data/nucleos.rds")
nucleos <- readRDS("data/nucleos.rds")

# IPCA e IPC-Br
# ipca <- window(BETS.get(433), start = c(1999,1), freq = 12)
# ipc  <- window(BETS.get(191), start = c(1999,1), freq = 12)
# saveRDS(ipc, "data/ipc.rds")
# saveRDS(ipca, "data/ipca.rds")
ipca <- window(readRDS("data/ipca.rds"), start = c(2001,1), freq = 12)
ipc <- window(readRDS("data/ipc.rds"), start = c(2001,1), freq = 12)

# gráfico dos núcleos - lado a lado -----------------------------------------
par(mar = c(2,4,1,2), mfrow = c(2,3))

# FGV
plot(ipc, main = "", lwd = 1, lty = 1, ylim = c(-0.5,3.5),
        col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"FGV"], lwd = 3, lty = 1,
     col = "#33A1C9", ylab = "variação mensal percentual (%)", xlab = "")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3, legend = c("IPC-Br","IPC-Br-MAS"), lwd = c(1,3), lty = c(1,1), y.intersp = 1.5,
       col = c(1,"#33A1C9"), cex = 1.5,bg = "white", box.col = "white",box.lwd = 0)

# IPCA - MAS
plot(ipca, main = "", lwd = 1, lty = 1, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"IPCA_MAS"], lwd = 3, lty = 1,
      col = "#33A1C9", ylab = "variação mensal percentual (%)", xlab = "")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3, legend = c("IPCA","Núcleo IPCA-MAS"), lwd = c(1,3), lty = c(1,1), y.intersp = 1.5,
       col = c(1,"#33A1C9"), cex = 1.5,bg = "white", box.col = "white",box.lwd = 0)

# IPCA - MA
plot(ipca, main = "", lwd = 1, lty = 1, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"IPCA_MA"], lwd = 3, lty = 1,
      col = "#33A1C9", ylab = "variação mensal percentual (%)", xlab = "")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3, legend = c("IPCA","Núcleo IPCA-MA"), lwd = c(1,3), lty = c(1,1), y.intersp = 1.5,
       col = c(1,"#33A1C9"), cex = 1.5,bg = "white", box.col = "white",box.lwd = 0)


# IPCA - EX1
plot(ipca, main = "", lwd = 1, lty = 1, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"IPCA_EX1"], lwd = 3, lty = 1,
      col = "#33A1C9", ylab = "variação mensal percentual (%)", xlab = "")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3, legend = c("IPCA","Núcleo IPCA-EX1"), lwd = c(1,3), lty = c(1,1), y.intersp = 1.5,
       col = c(1,"#33A1C9"), cex = 1.5,bg = "white", box.col = "white",box.lwd = 0)

# IPCA - EX2
plot(ipca, main = "", lwd = 1, lty = 1, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"IPCA_EX2"], lwd = 3, lty = 1,
      col = "#33A1C9", ylab = "variação mensal percentual (%)", xlab = "")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3, legend = c("IPCA","Núcleo IPCA-EX2"), lwd = c(1,3), lty = c(1,1), y.intersp = 1.5,
       col = c(1,"#33A1C9"), cex = 1.5,bg = "white", box.col = "white",box.lwd = 0)


# IPCA - DP
plot(ipca, main = "", lwd = 1, lty = 1, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleos[,"IPCA_DP"], lwd = 3, lty = 1,
      col = "#33A1C9", ylab = "variação mensal percentual (%)", xlab = "")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3, legend = c("IPCA","Núcleo IPCA-DP"), lwd = c(1,3), lty = c(1,1), y.intersp = 1.5,
       col = c(1,"#33A1C9"), cex = 1.5,bg = "white", box.col = "white",box.lwd = 0)


# diagnósticos -------------------------------------

# viés
bias_full <- rbind(core.diag(ipc, nucleos[,"FGV"], test = "bias")$out$full,
                   core.diag(ipca, nucleos[,"IPCA_MAS"], test = "bias")$out$full,
                   core.diag(ipca, nucleos[,"IPCA_MA"], test = "bias")$out$full[2,],
                   core.diag(ipca, nucleos[,"IPCA_EX1"], test = "bias")$out$full[2,],
                   core.diag(ipca, nucleos[,"IPCA_EX2"], test = "bias")$out$full[2,],
                   core.diag(ipca, nucleos[,"IPCA_DP"], test = "bias")$out$full[2,])

bias_part1 <- rbind(core.diag(ipc, nucleos[,"FGV"], test = "bias")$out$part1,
                    core.diag(ipca, nucleos[,"IPCA_MAS"], test = "bias")$out$part1,
                    core.diag(ipca, nucleos[,"IPCA_MA"], test = "bias")$out$part1[2,],
                    core.diag(ipca, nucleos[,"IPCA_EX1"], test = "bias")$out$part1[2,],
                    core.diag(ipca, nucleos[,"IPCA_EX2"], test = "bias")$out$part1[2,],
                    core.diag(ipca, nucleos[,"IPCA_DP"], test = "bias")$out$part1[2,])

bias_part2 <- rbind(core.diag(ipc, nucleos[,"FGV"], test = "bias")$out$part2,
                    core.diag(ipca, nucleos[,"IPCA_MAS"], test = "bias")$out$part2,
                    core.diag(ipca, nucleos[,"IPCA_MA"], test = "bias")$out$part2[2,],
                    core.diag(ipca, nucleos[,"IPCA_EX1"], test = "bias")$out$part2[2,],
                    core.diag(ipca, nucleos[,"IPCA_EX2"], test = "bias")$out$part2[2,],
                    core.diag(ipca, nucleos[,"IPCA_DP"], test = "bias")$out$part2[2,])

rownames(bias_full) = rownames(bias_part1) = rownames(bias_part2) <- c("IPC-Br", "Núcleo IPC-Br", "IPCA", colnames(nucleos)[-1])
bias_full
bias_part1
bias_part2

# teste de raiz unitária
adf <- core.diag(ipc, nucleos[,"FGV"], test = "root", lag_level = c(12,12,12,12,12,12), lag_diff = c(11,6,10,11,9,11))
adf <- core.diag(ipca, nucleos[,"IPCA_MAS"], test = "root",lag_level = c(7,7,12,12,7,12), lag_diff = c(11,6,7,11,6,11))
adf <- core.diag(ipca, nucleos[,"IPCA_MA"], test = "root",lag_level = c(12,7,12,7,7,7), lag_diff = c(11,6,7,6,6,6))
adf <- core.diag(ipca, nucleos[,"IPCA_EX1"], test = "root",lag_level = c(12,7,12,11,8,11), lag_diff = c(11,6,7,10,8,10))
adf <- core.diag(ipca, nucleos[,"IPCA_EX2"], test = "root",lag_level = c(12,7,12,11,12,8), lag_diff = c(11,6,7,10,11,7))
adf <- core.diag(ipca, nucleos[,"IPCA_DP"], test = "root",lag_level = c(12,7,12,7,7,7), lag_diff = c(11,6,7,6,6,12))


summary(adf$adf$full$y$level)
summary(adf$adf$part1$y$level)
summary(adf$adf$part2$y$level)
summary(adf$adf$full$core$level)
summary(adf$adf$part1$core$level)
summary(adf$adf$part2$core$level)
acf(adf$adf$full$y$level@res,36)
acf(adf$adf$part2$core$level@res,36)
# summary(adf$adf$full$y$diff)
# summary(adf$adf$part1$y$diff)
# summary(adf$adf$part2$y$diff)
summary(adf$adf$full$core$diff)
summary(adf$adf$part1$core$diff)
summary(adf$adf$part2$core$diff)
#acf(adf$adf$full$y$diff@res,36)
acf(adf$adf$part2$core$diff@res,36)


# teste de cointegração
coint <- core.diag(ipc, nucleos[,"FGV"], test = "coint")
coint <- core.diag(ipca, nucleos[,"IPCA_MAS"], test = "coint")
coint <- core.diag(ipca, nucleos[,"IPCA_MA"], test = "coint")
coint <- core.diag(ipca, nucleos[,"IPCA_EX1"], test = "coint")
coint <- core.diag(ipca, nucleos[,"IPCA_EX2"], test = "coint")
coint <- core.diag(ipca, nucleos[,"IPCA_DP"], test = "coint")

# atratividade
atratividade <- core.diag(ipc, nucleos[,"FGV"], test = "attract", 
                          lags_y = c(7,12), lags_core = c(1,2,3,12))

atratividade <- core.diag(ipca, nucleos[,"IPCA_MA"], test = "attract", 
                          lags_y = c(5,7), lags_core = c(1,2,5,6,7))

atratividade <- core.diag(ipca, nucleos[,"IPCA_MAS"], test = "attract", 
                          lags_y = c(7), lags_core = c(1,2,3,12))

atratividade <- core.diag(ipca, nucleos[,"IPCA_EX1"], test = "attract", 
                          lags_y = c(7), lags_core = c(1,2,3,4,5,6,7,8,9,11))

atratividade <- core.diag(ipca, nucleos[,"IPCA_EX2"], test = "attract", 
                          lags_y = c(1,2,3,5,6,7,8), lags_core = c(1,2,3,4,5,6,7,8,9,10,11))

atratividade <- core.diag(ipca, nucleos[,"IPCA_DP"], test = "attract", 
                          lags_y = c(5,7), lags_core = c(1,2,5,7))


# atratividade - part2
ipc0 <- window(ipc, start = c(2009,7), freq = 12)
ipca0 <- window(ipca, start = c(2009,7), freq = 12)
nucleos0 <- window(nucleos, start = c(2009,7), freq = 12)

atratividade <- core.diag(ipc0, nucleos0[,"FGV"], test = "attract", 
                          lags_y = c(12), lags_core = c(1,2,12))

atratividade <- core.diag(ipca0, nucleos0[,"IPCA_MA"], test = "attract", 
                          lags_y = c(5,7), lags_core = c(1,2,7))

atratividade <- core.diag(ipca0, nucleos0[,"IPCA_MAS"], test = "attract", 
                          lags_y = c(12), lags_core = c(1,2,3,12))

atratividade <- core.diag(ipca0, nucleos0[,"IPCA_EX1"], test = "attract", 
                          lags_y = c(12), lags_core = c(1,2,3,4,5,6,7,8,9,11))

atratividade <- core.diag(ipca0, nucleos0[,"IPCA_EX2"], test = "attract", 
                          lags_y = c(12), lags_core = c(1,2))

atratividade <- core.diag(ipca0, nucleos0[,"IPCA_DP"], test = "attract", 
                          lags_y = c(12), lags_core = c(1,2,5,6,7))


# sazonalidade

core.diag(ipc, nucleos[,"FGV"], test = "seas")
core.diag(ipca, nucleos[,"IPCA_MAS"], test = "seas")
core.diag(ipca, nucleos[,"IPCA_MA"], test = "seas")
core.diag(ipca, nucleos[,"IPCA_EX1"], test = "seas")
core.diag(ipca, nucleos[,"IPCA_EX2"], test = "seas")
core.diag(ipca, nucleos[,"IPCA_DP"], test = "seas")


# previsao

core.diag(ipc, nucleos[,"FGV"], test = "fcst")
core.diag(ipca, nucleos[,"IPCA_MA"], test = "fcst")
core.diag(ipca, nucleos[,"IPCA_MAS"], test = "fcst")
core.diag(ipca, nucleos[,"IPCA_EX1"], test = "fcst")
core.diag(ipca, nucleos[,"IPCA_EX2"], test = "fcst")
core.diag(ipca, nucleos[,"IPCA_DP"], test = "fcst")

core.diag(ipc0, nucleos0[,"FGV"], test = "fcst")
core.diag(ipca0, nucleos0[,"IPCA_MA"], test = "fcst")
core.diag(ipca0, nucleos0[,"IPCA_MAS"], test = "fcst")
core.diag(ipca0, nucleos0[,"IPCA_EX1"], test = "fcst")
core.diag(ipca0, nucleos0[,"IPCA_EX2"], test = "fcst") # parei aqui
core.diag(ipca0, nucleos0[,"IPCA_DP"], test = "fcst")

# PREVISÕES DE DIFERENÇA -----------------------

# nucleo-s
fcst1 <- core.fcst(ipc,  nucleos[,"FGV"], lag = 1, h = 12)
fcst3 <- core.fcst(ipc,  nucleos[,"FGV"], lag = 3, h = 12)
fcst6 <- core.fcst(ipc,  nucleos[,"FGV"], lag = 6, h = 12)
fcst12 <- core.fcst(ipc, nucleos[,"FGV"], lag = 12, h = 12)

round(sqrt(mean((ipc - lag(ipc, -1)  - fcst1$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -3)  - fcst3$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -6)  - fcst6$fcst)^2)),2)
round(sqrt(mean((ipc - lag(ipc, -12) - fcst12$fcst)^2)),2)


fcst1 <- core.fcst(ipca,  nucleos[,"IPCA_DP"], lag = 1, h = 12)
fcst3 <- core.fcst(ipca,  nucleos[,"IPCA_DP"], lag = 3, h = 12)
fcst6 <- core.fcst(ipca,  nucleos[,"IPCA_DP"], lag = 6, h = 12)
fcst12 <- core.fcst(ipca, nucleos[,"IPCA_DP"], lag = 12, h = 12)

round(sqrt(mean((ipca - lag(ipca, -1)  - fcst1$fcst)^2)),2)
round(sqrt(mean((ipca - lag(ipca, -3)  - fcst3$fcst)^2)),2)
round(sqrt(mean((ipca - lag(ipca, -6)  - fcst6$fcst)^2)),2)
round(sqrt(mean((ipca - lag(ipca, -12) - fcst12$fcst)^2)),2)
