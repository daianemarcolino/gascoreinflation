# pacotes, funções e leitura
library(betategarch)
library(dygraphs)
source("functions/betategarch_estimation.R")
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)

# estimar modelo: somente variância
m1 <- betategarch_estimation(ipc, initial = c(0.5,0.2,0.2,3), type = "var")
m1$otimizados
plot(m1$out)
ver <- cbind(ipc,m1$out[,"sigma"])
colnames(ver) <- c("ipc","sigma")
d1y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d1y

d1ep <- dygraph((m1$out[,"epsilon"] - mean(m1$out[,"epsilon"]))/sd(m1$out[,"epsilon"])) %>%
  dyOptions(colors = "black", strokeWidth = 1) %>% dySeries("V1", label = "epsilon")
d1ep

# estimar modelo: variância + média constante
m2 <- betategarch_estimation(ipc, initial = c(0.5,0.2,0.2,3,0.5), type = "mean-var")
plot(m2$out)
ver <- cbind(ipc,m2$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d2y <-dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d2y

d2ep <- dygraph((m2$out[,"epsilon"] - mean(m2$out[,"epsilon"]))/sd(m2$out[,"epsilon"])) %>%
  dyOptions(colors = "black", strokeWidth = 1) %>% dySeries("V1", label = "epsilon")
d2ep

# estimar modelo: variância + média 
m3 <- betategarch_estimation(ipc, initial = c(0.5,0.5,0.2,0.2,0.2,0.2,3), type = "mean-var2")
plot(m3$out)
ver <- cbind(ipc,m3$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d3y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d3y 

d3ep <- dygraph((m3$out[,"epsilon"] - mean(m3$out[,"epsilon"]))/sd(m3$out[,"epsilon"])) %>%
  dyOptions(colors = "black", strokeWidth = 1) %>% dySeries("V1", label = "epsilon")
d3ep

# teste de autocorrelação
Box.test(m1$out[,"epsilon"], lag = 12, type = "Ljung-Box")
Box.test(m2$out[,"epsilon"], lag = 12, type = "Ljung-Box")
Box.test(m3$out[,"epsilon"], lag = 12, type = "Ljung-Box")

# salvar resultados para shiny
resultados <- list(dygraphs_y = list(d1y = d1y, d2y = d2y, d3y = d3y),
                   dygraphs_ep = list(d1ep = d1ep, d2ep = d2ep, d3ep = d3ep),
                   estimation = list(m1 = m1, m2 = m2, m3 = m3))
saveRDS(resultados, "shiny_simulacao/data/resultados_betategarch_estimation.rds")
