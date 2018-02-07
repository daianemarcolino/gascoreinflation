
# funções e leitura -----------------------------

# leitura
ipc <- readRDS("data/ipc.rds")
nucleo_tf <- readRDS("data/nucleo_tf.rds")[,3]
nucleo_dcs <- window(readRDS("data/nucleo_dcs.rds"), end = c(2017,12), freq = 12)

# cálculo da combinação --------------------------------

# parte 1: calcular desvios segundo janela móvel de 48 meses

nucleos <- cbind(nucleo_tf, nucleo_dcs)
desvios <- NA*nucleos
for(i in 2:nrow(nucleos)){
  ifelse(i < 48, 
         desvios[i,] <- apply(nucleos[1:i,], MARGIN = 2, FUN = sd),
         desvios[i,] <- apply(nucleos[(i-47):i,], MARGIN = 2, FUN = sd))
}

desvios
novos_pesos <- (1/desvios)/rowSums(1/desvios) 
nucleo_comb <- ts(rowSums(nucleos * novos_pesos), start = start(nucleos), freq = 12)

ts.plot(nucleos, nucleo_comb, col = c(1,2,4))

# gráficos --------------------

par(mar = c(2,4,1,2), mfrow = c(1,1))
# IPC-Br vs. núcleo
plot(nucleo_dcs, main = "", lwd = 1, lty = 4, ylim = c(0,1.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(nucleo_tf, lwd = 2, lty = 1, col = "#1874CD")
lines(nucleo_comb, lwd = 2, lty = 1, col = "orangered")

abline(h = seq(0,1.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.5, legend = c("Núcleo-DCS","Núcleo-S","Núcleo-COMB"), lwd = c(1,2,2), lty = c(4,1,1), y.intersp = 1.5,
       col = c(1,"#1874CD","orangered"), cex = 1,bg = "white", box.col = "white",box.lwd = 0)

saveRDS(nucleo_comb, "data/nucleo_comb.rds")
