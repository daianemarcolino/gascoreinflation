library(betategarch)
library(GAS)
library(dygraphs)
source("functions/betategarchGAS_v2.R")

# exemplo NASDAQ ---------------------------------
data(nasdaq)
y <- ts(nasdaq[,2])

# usando o pacote betategarch
m <- tegarch(y, skew = F, asym = F)
summary(m)
fits <- fitted(m, verbose = T)
plot(fits[,"lambda"])
yhat <- ts(as.vector(fitted(m)), end = end(y), freq = 1)

data <- cbind(y,yhat)
graph_betategarch <- dygraph(data) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))
graph_betategarch 

# usando a rotina desenvolvida
k <- betategarchGAS(y, initial = c(0.02, 0.05, 0.95,10), type = "var")
yhat <- k[[1]]
k
graph_rotina <- dygraph(cbind(y,yhat)) %>% 
  dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))
graph_rotina

w1 <- m$par[1]
A1 <- m$par[3]
B1 <- m$par[2]
df <- m$par[4]

w1 <-  k$otimizados$par[1]
A1 <-  k$otimizados$par[3]
B1 <-  k$otimizados$par[2]
df <-  k$otimizados$par[4]

f2 <- w1/(1-B1)
N <- length(y)
u <- NULL
for(t in 1:N){
  u[t] <- ((df + 1)*y[t]^2) / (df*exp(2*f2[t]) + y[t]^2) - 1
  f2[t+1] <- w1 + A1*((df + 3)/(2*df))*u[t] + B1*f2[t]
}

f2 <- ts(f2[-1], start = start(y), freq = frequency(y))
loglik <-  N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + y^2/(df*exp(2*f2)))) + (N/2)*log((df-2)/df)

oi <- cbind((ts(exp(f2[-c(1:3)]))), y)
colnames(oi) <- c("yhat","y")
dygraph(oi)  %>% dySeries("y", label = "NASDAQ", color = "darkgrey", strokeWidth = 1) %>%
  dySeries("yhat", label = "Sigma2", color = "orangered", strokeWidth = 2) %>%
  dyRangeSelector(dateWindow = c(tail(as.Date(y),200)[1],tail(as.Date(y),1)))

