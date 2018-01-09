
# funções e leitura -----------------------------
source("functions/ifa_functions.R")

# leitura variações e pesos do IPC
variacao <- read.csv2("data/arquivosIPC/2 VARIAÇÃO IPC-DI.csv")
pesos <- read.csv2("data/arquivosIPC/1 PONDERAÇÃO IPC-DI.csv")
variacao.ts <- ipcts(variacao)
pesos.ts <- ipcts(pesos)
ipc <- variacao.ts$ipc

# cálculo do ifa --------------------------------

# filtro 1: médias aparadas
filtro1 <- core.ma(variacao.ts$subitens, pesos.ts$subitens, inf = 20, sup = 15, suave = T)

# média do IPC e do filtro 1
mean(ipc)
mean(filtro1)

# filtro 2: ajuste sazonal
filtro2 <- seas(filtro1, regression.aictest = NULL,
                transform.function = "none", 
                arima.model = "(1 1 0)(1 0 0)",
                series.modelspan = "2008.jan,2017.dec")

# diagnóstico
# summary(filtro2)
# qs(filtro2)
# plot(filtro2)
# ts.plot(filtro1,filtro2$series$s11,filtro2$series$s12, col = c("darkgray", 1,2), lwd = c(1),
#         lty = c(2,1,1))
# sliding <- series(filtro2, "sfs")
# sum(sliding[,"Max_._DIFF"] > 3, na.rm = T) / length(na.omit(sliding[,"Max_._DIFF"])) * 100
# monthplot(filtro2)

# núcleo com ajuste sazonal
filtro2 <- filtro2$series$s11

# filtro 3: médias móveis
filtro3 <- geom3(filtro2)

# exportar --------------------
ifa <- cbind(filtro1,filtro2,filtro3)
saveRDS(ifa, "data/nucleo_tf.rds")
