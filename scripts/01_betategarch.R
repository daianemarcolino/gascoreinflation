library(betategarch)
library(GAS)
source("functions/betategarchGAS_v1.R")

# exemplo NASDAQ ---------------------------------
data(nasdaq)
y <- ts(nasdaq[,2])

# usando o pacote betategarch
m <- tegarch(y, asym = F, skew = F)
summary(m)
yhat <- ts(as.vector(fitted(m)), end = end(y), freq = 1)

graph_betategarch <- dygraph(cbind(y,yhat)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))
graph_betategarch 

# usando a rotina desenvolvida
k <- betategarchGAS(y, initial.1 = c(10,2), initial.2 = c(1, 0.5, 0.3), type = "var")
yhat <- k[[1]]

graph_rotina <- dygraph(cbind(y,yhat)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))
graph_rotina

# usando o GAS
GASSpec <- UniGASSpec(Dist = "std", ScalingType = "Inv",
                      GASPar = list(location = FALSE, scale = TRUE, shape = FALSE))
Fit <- UniGASFit(GASSpec, y)
yhat <- ts(getFilteredParameters(Fit)[,"scale"], end = end(y), freq = frequency(y))

graph_gas <- dygraph(cbind(y,yhat)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))
graph_gas

# salvar pra usar no shiny
# saveRDS(list(graph_gas = graph_gas, graph_rotina = graph_rotina, graph_betategarch = graph_betategarch),
#         "shiny_simulacao/data/dygraphs_nasdaq.rds")


# exemplo ipc ---------------------------------
y <- window(BETS.get(191), start = c(2002,1), freq = 12)

# usando o pacote betategarch
m <- tegarch(y, asym = F, skew = F)
summary(m)
yhat <- ts(as.vector(fitted(m)), end = end(y), freq = 12)

graph_betategarch <- dygraph(cbind(y,yhat)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))
graph_betategarch 

# usando a rotina desenvolvida
k <- betategarchGAS(y, type = "var")
yhat <- k[[1]]

graph_rotina <- dygraph(cbind(y,yhat)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))
graph_rotina

# usando o GAS
GASSpec <- UniGASSpec(Dist = "std", ScalingType = "Inv",
                      GASPar = list(location = FALSE, scale = TRUE, shape = FALSE))
Fit <- UniGASFit(GASSpec, y)
yhat <- ts(getFilteredParameters(Fit)[,"scale"], end = end(y), freq = frequency(y))

graph_gas <- dygraph(cbind(y,yhat)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))
graph_gas

# salvar pra usar no shiny
saveRDS(list(graph_gas = graph_gas, graph_rotina = graph_rotina, graph_betategarch = graph_betategarch),
        "shiny_simulacao/data/dygraphs_nasdaq.rds")

m1 <- betategarchGAS(nasdaq[,2])
ts.plot(cbind(m1$y,m1$f2), col = 1:2)
hist(m1$f2)

library(BETS)
BETS.search("IPC")
ipc <- window(BETS.get(191), start = c(2002,1), freq = 12)
plot(ipc)

k <- (betategarchGAS(ipc, type = "var"))
ts.plot(ipc, k[[1]], col = 1:2)
plot(k[[1]])
k[[2]]
k[[3]]

ts.plot(ipc, k[[1]], col = 1:2)

k <-betategarchGAS(ipc, type = "mean-var")
plot(k[[1]][,1])
ts.plot(ipc, k[[1]][,1], col = 1:2)
ts.plot(ipc, k[[1]][,2], col = 1:2)
k[[2]]
m <- tegarch(ipc, asym = F, skew = F)
m$initial.values
m$lower
m$upper
ts(as.vector(fitted(m)))



k <-betategarchGAS(ipc, type = "var")
plot(k[[1]])
k[[2]]
dygraph(cbind(ipc,k[[1]])) %>% dyRangeSelector()
m <- tegarch(ipc, asym = F, skew = F)
summary(m)
yhat <- ts(as.vector(fitted(m)), end = end(ipc), freq = 12)
dygraph(cbind(ipc,yhat)) %>% dyRangeSelector()
dygraph(cbind(ipc,k[[1]],yhat)) %>% dyRangeSelector()
