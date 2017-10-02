library(betategarch)
library(dygraphs)
source("functions/betategarch_estimation.R")

# exemplo NASDAQ ---------------------------------
data(nasdaq)
y <- ts(nasdaq[,2])

# usando o pacote betategarch - modelo y[t] = exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
m <- tegarch(y, skew = F, asym = F)
summary(m)
yhat1 <- ts(as.vector(fitted(m)), end = end(y), freq = 1)

data <- cbind(y,yhat1)
graph_betategarch <- dygraph(data) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat1", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
graph_betategarch 

# usando a rotina desenvolvida - modelo y[t] = exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
k <- betategarch_estimation(y, initial = c(0.02, 0.05, 0.95,10), type = "var")
yhat2 <- k$out[,"sigma"]
k
graph_rotina <- dygraph(cbind(y,yhat2)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat2", label = "Sigma2", color = "steelblue", strokeWidth = 2) %>%
  dyRangeSelector()
graph_rotina

# exportar pro shiny
# graphs <- list(graph_betategarch = graph_betategarch, graph_rotina = graph_rotina)
# saveRDS(graphs, "shiny_simulacao/data/dygraphs_nasdaq.rds")

# usando a rotina desenvolvida - modelo y[t] = f1 + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
k <- betategarch_estimation(y, initial = c(0.02, 0.05, 0.95,10,5), type = "mean-var")
yhat3 <- k$out[,"sigma"]
k
graph_rotina2 <- dygraph(cbind(y,yhat3)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat3", label = "Sigma2", color = "steelblue", strokeWidth = 2) %>%
  dyRangeSelector()
graph_rotina2

# usando a rotina desenvolvida - modelo y[t] = f1[t] + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
k <- betategarch_estimation(y, initial = c(0.02,0.02,0.05, 0.05, 0.95,0.95,10), type = "mean-var2")
yhat4 <- cbind(y, k$out[,c("f1","sigma")])
colnames(yhat4) <- c("y","f1","f2")
k

graph_rotina3 <- dygraph(yhat4[,1:3]) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("f1", label = "Média", color = "orangered", strokeWidth = 2) %>%
  dySeries("f2", label = "Sigma2", color = "steelblue", strokeWidth = 2) %>%
  dyRangeSelector()
graph_rotina3

todos_sigmas <- cbind(yhat1, yhat2, yhat3)

graph_sigmas <- dygraph(cbind(yhat1,yhat2,yhat3)) %>% 
  dySeries("yhat1", label = "betategarch", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat2", label = "rotina_var", color = "steelblue", strokeWidth = 2) %>%
  dySeries("yhat3", label = "rotina_varmean", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
graph_sigmas
