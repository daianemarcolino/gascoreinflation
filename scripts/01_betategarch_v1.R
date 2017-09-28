library(betategarch)
library(GAS)
library(dygraphs)
source("betategarchGAS_v1.R")

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

