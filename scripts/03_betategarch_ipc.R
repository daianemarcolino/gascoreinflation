# pacotes, funções e leitura
library(betategarch)
library(dygraphs)
library(BETS)
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

# ACF
acf(m3$out[,"epsilon"], lag.max = 48)
m3$otimizados

# salvar resultados para shiny
resultados <- list(dygraphs_y = list(d1y = d1y, d2y = d2y, d3y = d3y),
                   dygraphs_ep = list(d1ep = d1ep, d2ep = d2ep, d3ep = d3ep),
                   estimation = list(m1 = m1, m2 = m2, m3 = m3))
saveRDS(resultados, "shiny_simulacao/data/resultados_betategarch_estimation.rds")



# estimar modelo: variância + média, s(t-1)
d <- data.frame(param = c("w1", "w2",  "A1_00","A1_01","B1_00","B1_01","A2_00","B2_00","df"),
                valor = c(0.21, -0.03,  0.05  ,0.05  ,0.57   ,0.02   ,0.97   , 4.73))
m4 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var3")
plot(m4$out)
ver <- cbind(ipc,m4$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d4y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d4y 

plot(m4$out[,"epsilon"])
acf(m4$out[,"epsilon"], lag.max = 48)
Box.test(m4$out[,"epsilon"], lag = 12, type = "Ljung-Box")
m4$otimizados

# estimar modelo: variância + média, s(t-1) f(t-1)
d <- data.frame(param = c("w1", "w2",  "A1_00","A1_01","B1_00","B1_01","A2_00","B2_00","df"),
                valor = c(0.29, -0.002,  0.04  ,0.01  ,0.40   ,0.5,   0.03   ,0.99   , 4.47))
m5 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var4")
plot(m5$out)
ver <- cbind(ipc,m5$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d5y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d5y 

plot(m5$out[,"epsilon"])
acf(m5$out[,"epsilon"], lag.max = 48)
Box.test(m5$out[,"epsilon"], lag = 12, type = "Ljung-Box")
m5$otimizados

# estimar modelo: variância + média, s(t-1) s(t-11) f(t-1)
d <- data.frame(param = c("w1", "w2",  "A1_00","A1_01","A1_11","B1_00","B1_01","A2_00","B2_00","df"),
                valor = c(0.29, -0.002,  0.04  ,0.01  ,0.05   , 0.40  ,0.5    ,0.03   ,0.99   , 4.47))
m6 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var5")
plot(m6$out)
ver <- cbind(ipc,m6$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d6y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d6y 

plot(m6$out[,"epsilon"])
acf(m6$out[,"epsilon"], lag.max = 48)
Box.test(m6$out[,"epsilon"], lag = 12, type = "Ljung-Box")
m6$otimizados


resp <- (m6$out[,"epsilon"] - mean(m6$out[,"epsilon"]))/sd(m6$out[,"epsilon"])
dygraph(resp) %>%
  dySeries("V1", label = "resíduo") %>%
  dyShading(from = -3, to = 3, axis = "y")

# OUTLIERS: JULHO DE 2000 E NOVEMBRO DE 2002
d1 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2000, month = 7)
d2 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2002, month = 11)
d3 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2015, month = 1)
d4 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2016, month = 1)
dummies <- cbind(d1,d2,d3,d4)

# estimar modelo: variância + média, s(t-1) s(t-11) f(t-1)
d <- data.frame(param = c("w1", "w2",  "A1_00","A1_01","A1_11","B1_00","B1_01","A2_00","B2_00","df","D1","D2","D3","D4"),
                valor = c(0.29, -0.002,  0.04  ,0.01  ,0.05   , 0.40  ,0.5    ,0.03   ,0.99   , 4.47, 1,1,1,1 ))
m7 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var6", dummy = dummies)
plot(m7$out)
ver <- cbind(ipc,m7$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d7y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d7y 

plot(m7$out[,"epsilon"])
acf(m7$out[,"epsilon"], lag.max = 48)
Box.test(m7$out[,"epsilon"], lag = 24, type = "Ljung-Box")
m7$otimizados


resp <- (m7$out[,"epsilon"] - mean(m7$out[,"epsilon"]))/sd(m7$out[,"epsilon"])
dygraph(resp) %>%
  dySeries("V1", label = "resíduo") %>%
  dyShading(from = -3, to = 3, axis = "y")

# previsão
prev <- betategarch_forecasting(out = m7, y = ipc, type = "mean-var6", dummy = dummies, h = 4)
ts.plot(ipc,prev[,1], col = 1:2)


# EXERCÍCIO DENTRO DA AMOSTRA
y0 <- window(ipc, end = c(2016,12), freq = 12)

# OUTLIERS: JULHO DE 2000 E NOVEMBRO DE 2002
d1 <- BETS.dummy(start = start(y0), end = end(y0), year = 2000, month = 7)
d2 <- BETS.dummy(start = start(y0), end = end(y0), year = 2002, month = 11)
d3 <- BETS.dummy(start = start(y0), end = end(y0), year = 2015, month = 1)
d4 <- BETS.dummy(start = start(y0), end = end(y0), year = 2016, month = 1)
dummies <- cbind(d1,d2,d3,d4)

# estimar modelo: variância + média, s(t-1) s(t-11) f(t-1)
d <- data.frame(param = c("w1", "w2",  "A1_00","A1_01","A1_11","B1_00","B1_01","A2_00","B2_00","df","D1","D2","D3","D4"),
                valor = c(0.29, -0.002,  0.04  ,0.01  ,0.05   , 0.40  ,0.5    ,0.03   ,0.99   , 4.47, 1,1,1,1 ))
m <- betategarch_estimation(y0, initial = d[,2], type = "mean-var6", dummy = dummies)
prev <- betategarch_forecasting(out = m, y = y0, type = "mean-var6", dummy = dummies, h = 8)
ts.plot(ipc,prev[,1], col = 1:2)
