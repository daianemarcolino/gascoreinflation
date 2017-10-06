# pacotes, funções e leitura
library(betategarch)
library(dygraphs)
library(BETS)
source("functions/betategarch_estimation.R")
source("functions/betategarch_forecasting.R")
source("functions/ERRO.R")
ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)

# MODELO 1 ------------------------------------------------------------------------------------
# estimar modelo: variância + média | s(t) f(t)
m1 <- betategarch_estimation(ipc, initial = c(0.5,0.5,0.2,0.2,0.2,0.2,3), type = "mean-var2")
plot(m1$out)
ver <- cbind(ipc,m1$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d1y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d1y 

resp1 <- (m1$out[,"epsilon"] - mean(m1$out[,"epsilon"]))/sd(m1$out[,"epsilon"])
d1ep <- dygraph(resp1) %>%
  dySeries("V1", label = "resíduo", color = "#AD2D1F", strokeWidth = 2) %>%
  dyShading(from = -3, to = 3, axis = "y", color = "#F0DFDB")
d1ep

# ACF
acf(resp1, lag.max = 96, ci.col = "red", main = "", xaxt = "n")
axis(1, at = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9)*12)
title(main = list("Função de autocorrelação nos resíduos", font = 1, cex = 0.9 ))
grid(nx = NA, ny = NULL)

# # previsão fora da amostra simples 
# prev1 <- betategarch_forecasting(out = m1, y = ipc, type = "mean-var2", h = 4)
# prev1. <- ts(c(tail(ipc,1), prev1[,1]), end = end(prev1), freq = 12)
# ts.plot(ipc, prev1., col = 1:2)
# 
# # previsão dentro da amostra
# data <- data.frame(ano = c(rep(2016,12),rep(2017,7)), mes = c(1:12,1:7))
# initial <- m1$otimizados$par
# prev <- NULL
# for(i in 1:nrow(data)){
#   y0 <- window(ipc, end = c(data$ano[i],data$mes[i]), freq = 12)
#   m <- betategarch_estimation(y0, initial = initial, type = "mean-var2")
#   if(m$otimizados$message != "relative convergence (4)"){
#     m <- betategarch_estimation(y0, initial = c(0.2,0.2,0.2,0.2,0.2,0.2,3), type = "mean-var2")
#   }
#   prev <- rbind(prev, betategarch_forecasting(out = m, y = y0, type = "mean-var2", h = 1))
#   initial <- m$otimizados$par
# }
# prev_m1 <- na.omit(cbind(ipc, ts(prev[,1], start = c(data$ano[2], data$mes[2]), freq = 12)))
# colnames(prev_m1) <- c("ipc","prev")
# ts.plot(prev_m1, col = c(1,"orangered"), lty = c(3,1))
# ERRO(prev_m1[,"ipc"],prev_m1[,"prev"])

# MODELO 2 ------------------------------------------------------------------------------------
# estimar modelo: variância + média, s(t-1) f(t-1)
d <- data.frame(param = c("w1", "w2",  "A1_00","A1_01","B1_00","B1_01","A2_00","B2_00","df"),
                valor = c(0.29, -0.002,  0.04  ,0.01  ,0.40   ,0.5,   0.03   ,0.99   , 4.47))
m2 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var4")
plot(m2$out)
ver <- cbind(ipc,m2$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d2y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d2y 

resp2 <- (m2$out[,"epsilon"] - mean(m2$out[,"epsilon"]))/sd(m2$out[,"epsilon"])
d2ep <- dygraph(resp2) %>%
  dySeries("V1", label = "resíduo", color = "#AD2D1F", strokeWidth = 2) %>%
  dyShading(from = -3, to = 3, axis = "y", color = "#F0DFDB")
d2ep

# ACF
acf(resp2, lag.max = 96, ci.col = "red", main = "", xaxt = "n")
axis(1, at = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9)*12)
title(main = list("Função de autocorrelação nos resíduos", font = 1, cex = 0.9 ))
grid(nx = NA, ny = NULL)

# # previsão fora da amostra simples 
# prev2 <- betategarch_forecasting(out = m2, y = ipc, type = "mean-var4", h = 4)
# prev2. <- ts(c(tail(ipc,1), prev2[,1]), end = end(prev2), freq = 12)
# ts.plot(ipc, prev2., col = 1:2)
# 
# # previsão dentro da amostra
# data <- data.frame(ano = c(rep(2016,12),rep(2017,7)), mes = c(1:12,1:7))
# initial <- m2$otimizados$par
# prev <- NULL
# for(i in 1:nrow(data)){
#   y0 <- window(ipc, end = c(data$ano[i],data$mes[i]), freq = 12)
#   m <- betategarch_estimation(y0, initial = initial, type = "mean-var4")
#   if(!(m$otimizados$message %in% c("relative convergence (4)", "X-convergence (3)"))){
#     m <- betategarch_estimation(y0, initial = c(0.3,-0.03,0.04,0.02,0.06,0.21,0.03,0.9,5), type = "mean-var4")
#     message(m$otimizados$par)
#   }
#   prev <- rbind(prev, betategarch_forecasting(out = m, y = y0, type = "mean-var4", h = 1))
#   initial <- m$otimizados$par
# }
# prev_m2 <- na.omit(cbind(ipc, ts(prev[,1], start = c(data$ano[2], data$mes[2]), freq = 12)))
# colnames(prev_m2) <- c("ipc","prev")
# ts.plot(prev_m2, col = c(1,"orangered"), lty = c(3,1))
# ERRO(prev_m2[,"ipc"],prev_m2[,"prev"])

# MODELO 3 ------------------------------------------------------------------------------------
# estimar modelo: variância + média, s(t-1) s(t-11) f(t-1) f(t-11)
d <- data.frame(param = c("w1","w2"  ,"A1_00","A1_01","A1_11","B1_00","B1_01","B1_11","A2_00","B2_00","df"),
                valor = c(0.29,-0.002,0.04   ,0.01   ,0.05   ,0.40   ,0.5    ,0.05   ,0.03   ,0.99   , 4.47))
m3 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var7")
plot(m3$out)
ver <- cbind(ipc,m3$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d3y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d3y 

resp3 <- (m3$out[,"epsilon"] - mean(m3$out[,"epsilon"]))/sd(m3$out[,"epsilon"])
d3ep <- dygraph(resp3) %>%
  dySeries("V1", label = "resíduo", color = "#AD2D1F", strokeWidth = 2) %>%
  dyShading(from = -3, to = 3, axis = "y", color = "#F0DFDB")
d3ep

# ACF
acf(resp3, lag.max = 96, ci.col = "red", main = "", xaxt = "n")
axis(1, at = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9)*12)
title(main = list("Função de autocorrelação nos resíduos", font = 1, cex = 0.9 ))
grid(nx = NA, ny = NULL)

# # previsão fora da amostra simples 
# prev3 <- betategarch_forecasting(out = m3, y = ipc, type = "mean-var7", h = 4)
# prev3. <- ts(c(tail(ipc,1), prev3[,1]), end = end(prev3), freq = 12)
# ts.plot(ipc, prev3., col = 1:2)
# 
# # previsão dentro da amostra
# data <- data.frame(ano = c(rep(2016,12),rep(2017,7)), mes = c(1:12,1:7))
# initial <- m3$otimizados$par
# prev <- NULL
# for(i in 1:nrow(data)){
#   y0 <- window(ipc, end = c(data$ano[i],data$mes[i]), freq = 12)
#   m <- betategarch_estimation(y0, initial = initial, type = "mean-var7")
#   if(!(m$otimizados$message %in% c("relative convergence (4)", "X-convergence (3)"))){
#     m <- betategarch_estimation(y0, initial = m3$otimizados$par, type = "mean-var7")
#     message(m$otimizados$par)
#   }
#   prev <- rbind(prev, betategarch_forecasting(out = m, y = y0, type = "mean-var7", h = 1))
#   initial <- m$otimizados$par
# }
# prev_m3 <- na.omit(cbind(ipc, ts(prev[,1], start = c(data$ano[2], data$mes[2]), freq = 12)))
# colnames(prev_m3) <- c("ipc","prev")
# ts.plot(prev_m3, col = c(1,"orangered"), lty = c(3,1))
# ERRO(prev_m3[,"ipc"],prev_m3[,"prev"])

# MODELO 4 ------------------------------------------------------------------------------------
# estimar modelo: variância + média, s(t-1) s(t-11) f(t-1) f(t-11) c/ dummies

# OUTLIERS: JULHO DE 2000 E NOVEMBRO DE 2002
d1 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2002, month = 11)
d2 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2000, month = 7)
d3 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2015, month = 1)
d4 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2016, month = 1)
dummies <- cbind(d1,d2,d3,d4)

d <- data.frame(param = c("w1","w2"  ,"A1_00","A1_01","A1_11","B1_00","B1_01","B1_11","A2_00","B2_00","df", "D1","D2","D3"),
                valor = c(0.17,-0.02 ,0.04    ,-0.003,0.03   ,0.53   ,0.10   ,-0.03  ,0.01   ,0.7   ,8.80,2.36,1.35,1))

m4 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var8", dummy = dummies)
plot(m4$out)
ver <- cbind(ipc,m4$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d4y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d4y 

resp4 <- (m4$out[,"epsilon"] - mean(m4$out[,"epsilon"]))/sd(m4$out[,"epsilon"])
d4ep <- dygraph(resp4) %>%
  dySeries("V1", label = "resíduo", color = "#AD2D1F", strokeWidth = 2) %>%
  dyShading(from = -3, to = 3, axis = "y", color = "#F0DFDB")
d4ep

# ACF
acf(resp4, lag.max = 96, ci.col = "red", main = "", xaxt = "n")
axis(1, at = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9)*12)
title(main = list("Função de autocorrelação nos resíduos", font = 1, cex = 0.9 ))
grid(nx = NA, ny = NULL)

# # previsão fora da amostra simples 
# prev4 <- betategarch_forecasting(out = m4, y = ipc, type = "mean-var8", h = 4)
# prev4. <- ts(c(tail(ipc,1), prev4[,1]), end = end(prev4), freq = 12)
# ts.plot(ipc, prev4., col = 1:2)
# 
# # previsão dentro da amostra
# data <- data.frame(ano = c(rep(2016,12),rep(2017,7)), mes = c(1:12,1:7))
# initial <- m4$otimizados$par
# prev <- NULL
# for(i in 1:nrow(data)){
#   y0 <- window(ipc, end = c(data$ano[i],data$mes[i]), freq = 12)
#   d1 <- BETS.dummy(start = start(y0), end = end(y0), year = 2002, month = 11)
#   d2 <- BETS.dummy(start = start(y0), end = end(y0), year = 2000, month = 7)
#   d3 <- BETS.dummy(start = start(y0), end = end(y0), year = 2015, month = 1)
#   d4 <- BETS.dummy(start = start(y0), end = end(y0), year = 2016, month = 1)
#   dummies <- cbind(d1,d2,d3,d4)
#   
#   m <- betategarch_estimation(y0, initial = initial, dummy = dummies, type = "mean-var8")
#   if(!(m$otimizados$message %in% c("relative convergence (4)", "X-convergence (3)"))){
#     m <- betategarch_estimation(y0, initial = c(0.17,-0.02 ,0.04,-0.003,0.03,0.53,0.10,-0.03,0.01,0.7,8.80,2.36,1.35,1), 
#                                 type = "mean-var8", dummy = dummies)
#     message(m$otimizados$par)
#   }
#   prev <- rbind(prev, betategarch_forecasting(out = m, y = y0, type = "mean-var8", h = 1))
#   initial <- m$otimizados$par
# }
# prev_m4 <- na.omit(cbind(ipc, ts(prev[,1], start = c(data$ano[2], data$mes[2]), freq = 12)))
# colnames(prev_m4) <- c("ipc","prev")
# ts.plot(prev_m4, col = c(1,"orangered"), lty = c(3,1))
# ERRO(prev_m4[,"ipc"],prev_m4[,"prev"])

# MODELO 5 ------------------------------------------------------------------------------------
# estimar modelo: variância + média, s(t-1) s(t-11) s(t-12) f(t-1) f(t-11) s/ dummies

# OUTLIERS: JULHO DE 2000 E NOVEMBRO DE 2002
d1 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2002, month = 11)
d2 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2000, month = 7)
d3 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2015, month = 1)
d4 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2016, month = 1)
dummies <- cbind(d1,d2,d3,d4)

d <- data.frame(param = c("w1","w2"  ,"A1_00","A1_01","A1_11","A1_12","B1_00","B1_01","B1_11","A2_00","B2_00","df"),
                valor = c(0.17,-0.55 ,0.01   ,-0.01  ,0.01   ,-0.02   ,0.3    ,0.1   ,0.08   ,0.16   ,0.6  ,8 ))

m5 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var9")
plot(m5$out)
ver <- cbind(ipc,m5$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d5y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d5y 

resp5 <- (m5$out[,"epsilon"] - mean(m5$out[,"epsilon"]))/sd(m5$out[,"epsilon"])
d5ep <- dygraph(resp5) %>%
  dySeries("V1", label = "resíduo", color = "#AD2D1F", strokeWidth = 2) %>%
  dyShading(from = -3, to = 3, axis = "y", color = "#F0DFDB")
d5ep

# ACF
acf(resp5, lag.max = 96, ci.col = "red", main = "", xaxt = "n")
axis(1, at = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9)*12)
title(main = list("Função de autocorrelação nos resíduos", font = 1, cex = 0.9 ))
grid(nx = NA, ny = NULL)

# MODELO 6 ------------------------------------------------------------------------------------
# estimar modelo: variância + média, s(t-1) s(t-11) s(t-12) f(t-1) f(t-11) s/ dummies

# OUTLIERS: JULHO DE 2000 E NOVEMBRO DE 2002
d1 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2002, month = 11)
d2 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2003, month = 1)
d3 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2000, month = 7)
d4 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2015, month = 1)
d5 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2016, month = 1)
d6 <- BETS.dummy(start = start(ipc), end = end(ipc), year = 2010, month = 1)
dummies <- cbind(d1,d2,d3,d4,d5,d6)

d <- data.frame(param = c("w1","w2"  ,"A1_00","A1_01","A1_11","A1_12","B1_00","B1_01","B1_11","A2_00","B2_00","df","D1","D2","D3","D4","D5","D6"),
                valor = c(0.12,-0.03 ,0.02   ,0.008  ,0.03   ,-0.02   ,0.4   ,-0.21  ,0.68   ,0.12   ,0.71   ,12  ,2   ,1   ,1   ,1   ,1   ,1))

m6 <- betategarch_estimation(ipc, initial = d[,2], type = "mean-var10", dummy = dummies)
plot(m6$out)
ver <- cbind(ipc,m6$out[,c("f1","sigma")])
colnames(ver) <- c("ipc","f1","sigma")
d6y <- dygraph(ver) %>%
  dySeries("ipc", label = "FGV IPC", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "steelblue", strokeWidth = 2) %>%
  dySeries("sigma", label = "Sigma", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
d6y 

resp6 <- (m6$out[,"epsilon"] - mean(m6$out[,"epsilon"]))/sd(m6$out[,"epsilon"])
d6ep <- dygraph(resp6) %>%
  dySeries("V1", label = "resíduo", color = "#AD2D1F", strokeWidth = 2) %>%
  dyShading(from = -3, to = 3, axis = "y", color = "#F0DFDB") %>%
  dyAxis(valueRange = c(-4,4), name = "y")
d6ep

# ACF
acf(resp6, lag.max = 96, ci.col = "red", main = "", xaxt = "n")
axis(1, at = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9)*12)
title(main = list("Função de autocorrelação nos resíduos", font = 1, cex = 0.9 ))
grid(nx = NA, ny = NULL)


# previsão fora da amostra simples
prev6 <- betategarch_forecasting(out = m6, y = ipc, type = "mean-var10", h = 4)
prev6. <- ts(c(tail(ipc,1), prev6[,1]), end = end(prev6), freq = 12)
ts.plot(ipc, prev6., col = 1:2)

# previsão dentro da amostra
data <- data.frame(ano = c(rep(2016,12),rep(2017,7)), mes = c(1:12,1:7))
initial <- m6$otimizados$par
prev <- NULL
for(i in 1:nrow(data)){
  y0 <- window(ipc, end = c(data$ano[i],data$mes[i]), freq = 12)
  d1 <- BETS.dummy(start = start(y0), end = end(y0), year = 2002, month = 11)
  d2 <- BETS.dummy(start = start(y0), end = end(y0), year = 2003, month = 1)
  d3 <- BETS.dummy(start = start(y0), end = end(y0), year = 2000, month = 7)
  d4 <- BETS.dummy(start = start(y0), end = end(y0), year = 2015, month = 1)
  d5 <- BETS.dummy(start = start(y0), end = end(y0), year = 2016, month = 1)
  d6 <- BETS.dummy(start = start(y0), end = end(y0), year = 2010, month = 1)
  dummies <- cbind(d1,d2,d3,d4,d5,d6)

  m <- betategarch_estimation(y0, initial = initial, dummy = dummies, type = "mean-var10")
  if(!(m$otimizados$message %in% c("relative convergence (4)", "X-convergence (3)"))){
    m <- betategarch_estimation(y0, initial = c(0.12,-0.03,0.02,0.008,0.03,-0.02,0.4,-0.21,0.68,0.12,0.71,12,2,1,1,1,1,1),
                                type = "mean-var10", dummy = dummies)
    message(m$otimizados$par)
  }
  prev <- rbind(prev, betategarch_forecasting(out = m, y = y0, type = "mean-var10", h = 1))
  initial <- m$otimizados$par
}
prev_m6 <- na.omit(cbind(ipc, ts(prev[,1], start = c(data$ano[2], data$mes[2]), freq = 12)))
colnames(prev_m6) <- c("ipc","prev")
ts.plot(prev_m6, col = c(1,"orangered"), lty = c(3,1))
ERRO(prev_m6[,"ipc"],prev_m6[,"prev"])

# EXPORTAR PARA SHINY -----------------------------------------------

resultados <- list(graphs_y = list(d1y = d1y, d2y = d2y, d3y = d3y, d4y = d4y, d5y = d5y, d6y = d6y),
                   graphs_resp = list(d1ep = d1ep, d2ep = d2ep, d3ep = d3ep, d4ep = d4ep, d5ep = d5ep, d6ep = d6ep),
                   residuos = list(resp1 = resp1, resp2 = resp2, resp3 = resp3, resp4 = resp4, resp5 = resp5, resp6 = resp6),
                   previsao = list(fora = prev6., dentro = prev_m6),
                   modelos = list(m6 = m6)
                   )
saveRDS(resultados, "shiny_simulacao/data/novos_modelos.rds")

