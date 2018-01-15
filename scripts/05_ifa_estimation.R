
# funções e leitura -----------------------------
source("functions/ifa_functions.R")
library(seasonal)

# leitura variações e pesos do IPC
variacao <- read.csv2("data/arquivosIPC/2 VARIAÇÃO IPC-DI.csv")
pesos <- read.csv2("data/arquivosIPC/1 PONDERAÇÃO IPC-DI.csv")
variacao.ts <- ipcts(variacao)
pesos.ts <- ipcts(pesos)
ipc <- window(variacao.ts$ipc, start = c(2001,1), freq = 12)
variacao.ts$subitens <- window(variacao.ts$subitens, start = c(2000,11), freq = 12)
pesos.ts$subitens <- window(pesos.ts$subitens, start = c(2000,11), freq = 12)
# saveRDS(ipc, "data/ipc.rds")

# cálculo do ifa --------------------------------

# filtro 1: médias aparadas
filtro1 <- list()
filtro2 <- list()
filtro3 <- list()
medias <- integer(11)
for(i in 10:20){
  filtro1[[i-9]] <- core.ma(variacao.ts$subitens, pesos.ts$subitens, inf = 20, sup = i, suave = T)
  filtro2[[i-9]] <- seas(filtro1[[i-9]], regression.aictest = NULL,
                         transform.function = "none", 
                         arima.model = "(1 1 0)(1 0 0)",
                         series.modelspan = "2009.jan,2017.dec")$series$s11
  filtro3[[i-9]] <- window(geom3(filtro2[[i-9]]), start = c(2001,1), freq = 12)
  print(i)
}
medias <- sapply(filtro3, FUN = mean, na.rm = T)
medias
mean(ipc)
i <- (10:20)[6]

# # diag filtro 2: ajuste sazonal
# f <- seas(filtro1[[7]], regression.aictest = NULL,
#           transform.function = "none",
#           arima.model = "(1 1 0)(1 0 0)",
#           series.modelspan = "2009.jan,2017.dec")
# # diagnóstico
# summary(f)
# qs(f)
# plot(f)
# sliding <- series(f, "sfs")
# sum(sliding[,"Max_._DIFF"] > 3, na.rm = T) / length(na.omit(sliding[,"Max_._DIFF"])) * 100
# monthplot(f)


# filtros finais
filtro1_final <- filtro1[[6]]
filtro2_final <- filtro2[[6]]
filtro3_final <- filtro3[[6]]


# exportar --------------------
# ifa <- window(cbind(filtro1_final,filtro2_final,filtro3_final), start = c(2001,1),freq = 12)
# colnames(ifa) <- paste0("filtro",1:3)
# saveRDS(ifa, "data/nucleo_tf.rds")
ifa <- readRDS("data/nucleo_tf.rds")
# gráficos --------------------

par(mar = c(2,4,1,2), mfrow = c(1,1))
# IPC-Br vs. núcleo
plot(window(ipc, start = c(2001,1), freq = 12), main = "", lwd = 1, lty = 4, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(ifa[,"filtro3"], lwd = 2, lty = 1, col = "#1874CD")
abline(h = 0:3, col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
#abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3, legend = c("IPC-Br","Núcleo-S"), lwd = c(1,2), lty = c(4,1), y.intersp = 1.5,
       col = c(1,"#1874CD"), cex = 1.3,bg = "white", box.col = "white",box.lwd = 0)


# núcleo: os três filtros
plot(ifa[,1], main = "", lwd = 1, lty = 4, # ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(ifa[,2], lwd = 2, lty = 5, col = "orangered")
lines(ifa[,3], lwd = 2, lty = 1, col = "#1874CD")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,1.5, legend = c("passo 1: médias ap.","passo 2: aj. sazonal","passo 3: núcleo final"), lwd = c(1,2,2), lty = c(4,5,1), y.intersp = 1.5,
       col = c(1,"orangered","#1874CD"), cex = 1.3,bg = "white", box.col = "white",box.lwd = 0)

