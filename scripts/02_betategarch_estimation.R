library(betategarch)
library(dygraphs)
source("functions/betategarch_estimation.R")

# exemplo NASDAQ ---------------------------------
data(nasdaq)
y <- ts(nasdaq[,2])

# usando o pacote betategarch
m <- tegarch(y, skew = F, asym = F)
summary(m)
yhat1 <- ts(as.vector(fitted(m)), end = end(y), freq = 1)

data <- cbind(y,yhat1)
graph_betategarch <- dygraph(data) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat1", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector()
graph_betategarch 

# usando a rotina desenvolvida
k <- betategarchGAS(y, initial = c(0.02, 0.05, 0.95,10), type = "var")
yhat2 <- k[[1]]
k
graph_rotina <- dygraph(cbind(y,yhat2)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat2", label = "Sigma2", color = "steelblue", strokeWidth = 2) %>%
  dyRangeSelector()
graph_rotina

# exportar pro shiny
graphs <- list(graph_betategarch = graph_betategarch, graph_rotina = graph_rotina)
saveRDS(graphs, "shiny_simulacao/data/dygraphs_nasdaq.rds")
