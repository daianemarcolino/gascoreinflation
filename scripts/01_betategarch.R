library(betategarch)

data(nasdaq)
y <- zoo(nasdaq[,"nasdaqret"], order.by=as.Date(nasdaq[,"day"], "%Y-%m-%d"))
plot(y, main="The Nasdaq 100 index (daily)", xlab="", ylab="Log-return in %")
nasdaq1comp <- tegarch(y)
nasdaq1stdev <- data.frame(fitted(nasdaq1comp))[,1]
plot(nasdaq1stdev, main="", ylab="1-comp: St.dev.", xlab="")
y=c(as.data.frame(y))
ts.plot(cbind(nasdaq[,2], c(nasdaq1stdev)), col = 1:2)
ts.plot(c(y),c(nasdaq1stdev))
nasdaq1comp <- tegarch(y, asym = F, skew = F, components = 1)
nasdaq1comp$lower
nasdaq1comp$upper
nasdaq1comp$par
nasdaq1comp$objective

fit <- c(nasdaq1stdev)

library(GAS)
GASSpec <- UniGASSpec(Dist = "std", ScalingType = "Inv",
                      GASPar = list(location = FALSE, scale = TRUE, shape = FALSE))
GASSpec
Fit <- UniGASFit(GASSpec, nasdaq[,2])
Fit

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

data(nasdaq)
y <- ts(nasdaq[,2])
k <-betategarchGAS(y, type = "var")
plot(k[[1]])
dygraph(cbind(y,k[[1]])) %>% dyRangeSelector()
m <- tegarch(y, asym = F, skew = F)
summary(m)
yhat <- ts(as.vector(fitted(m)), end = end(y), freq = 1)
dygraph(cbind(y,yhat)) %>% dyRangeSelector()
dygraph(cbind(y,log(k[[1]]),yhat)) %>% dyRangeSelector()

k <-betategarchGAS(ipc, type = "var")
plot(k[[1]])
k[[2]]
dygraph(cbind(ipc,k[[1]])) %>% dyRangeSelector()
m <- tegarch(ipc, asym = F, skew = F)
summary(m)
yhat <- ts(as.vector(fitted(m)), end = end(ipc), freq = 12)
dygraph(cbind(ipc,yhat)) %>% dyRangeSelector()
dygraph(cbind(ipc,k[[1]],yhat)) %>% dyRangeSelector()
