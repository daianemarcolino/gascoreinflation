diag.dcs <- function(y, out){

  df <- out$otimizados$par[4]
  ncp <- ts(rowSums(out$out[,c("mu","gamma")]), start = start(out$out), freq = 12)
  ep <- out$out[,"epsilon"]
  # resíduo quantílico
  
  Ft <- pt(y, df = df, ncp = ncp)
  rq <- qnorm(Ft)
  plot(ep,rq, col = 2:1)
  ts.plot(rq,ep, col = 2:1)
  
  # gráfico resíduo de pearson vs. resíduo quantílico
  #ts.plot(y,rq, lty = c(3,1), col = c(1,"dodgerblue"))
  legend("top", legend = c("RP","RQ"), lty = c(3,1), col = c(1,"dodgerblue"), cex = 0.8, bty = "n")
  
  x <- data.frame(matrix(NA, ncol = 7, nrow = 2))
  colnames(x) <- c("Assimetria","Curtose","Média","Mediana","Desvio-padrão","AD.stat","AD.pvalue")
  rownames(x) <- c("y","rq")
  
  x["y",] <- c(mean((ep - mean(ep))^3/sd(ep)^3), mean((ep - mean(ep))^4/sd(ep)^4),
               mean(ep), median(ep), sd(ep), nortest::ad.test(ep)$statistic, nortest::ad.test(ep)$p.value)
  x["rq",] <- c(mean((rq - mean(rq))^3/sd(rq)^3), mean((rq - mean(rq))^4/sd(rq)^4),
                mean(rq), median(rq), sd(rq), nortest::ad.test(rq)$statistic, nortest::ad.test(rq)$p.value)
  x

  r1 <- rt(n = length(y), ncp = mean(ncp), df = df)
  r2 <- rnorm(n = length(rq))
  par(mfrow = c(1,2))
  qqplot(r1,y, main = "Resíduo de Pearson")#, xlim = c(-4.5,4.5), ylim = c(-4.5,4.5))
  qqline(y, distribution = function(x) qt(x, ncp = mean(ncp), df = df))
  qqplot(c(rq),r2, main = "Resíduo quantílico")#, xlim = c(-3,3), ylim = c(-3,3))
  qqline(c(r2))
  
  hist(rq)
  # output
  print(round(out,4))
  list(stats = round(out,4), resid.q = rq)
}