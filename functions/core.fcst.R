core.fcst <- function(y, core, lag = 1, h = 24){
  
  data <- na.omit(cbind(y,core))
  y <- data[,1]
  core <- data[,2]

  # inflação (pi) defasada
  pi_lag <- lag(y, k = -lag)
  
  # nucleo (pic) defasado
  pic_lag <- lag(core, k = -lag)
  
  # unir
  tudo <- na.omit(cbind(y = y - pi_lag, x = pic_lag - pi_lag))
  
  # análise h meses fora da amostra
  fcst <- integer(h)
  for(i in h:1){
    tudo0 <- head(tudo, nrow(tudo) - i + 1)
    x_fcst <- tail(tudo0,1)[,"x"]
    reg <- lm(y ~ x, data = head(tudo0,nrow(tudo0)-1))
    fcst[i] <- predict(reg, x_fcst)
  }
  fcst <- ts(fcst[h:1], end = end(tudo), freq = 12)
  
  ts.plot(na.omit(cbind(tudo[,"y"],fcst)), col = c(1,2), lty = c(3,1))
  legend("bottomleft", legend = c("y","fcst"), col = c(1,2), lty = c(3,1), bty = "n")
  
 
  # output
  list(fcst = fcst)
  
}

