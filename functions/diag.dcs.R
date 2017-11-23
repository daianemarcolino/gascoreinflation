diag.dcs <- function(y, df){
  # y = bsm0$out[,"epsilon"]
  # df = bsm0$otimizados$par[4]
  
  # resíduo quantílico
  Ft <- pt(y, df = df)
  rq <- qnorm(Ft)
  
  # gráfico resíduo de pearson vs. resíduo quantílico
  ts.plot(y,rq, lty = c(3,1), col = c(1,"dodgerblue"))
  legend("top", legend = c("RP","RQ"), lty = c(3,1), col = c(1,"dodgerblue"), cex = 0.8, bty = "n")
  
  #
  
  out <- data.frame(matrix(NA, ncol = 7, nrow = 2))
  colnames(out) <- c("Assimetria","Curtose","Média","Mediana","Desvio-padrão","AD.stat","AD.pvalue")
  rownames(out) <- c("y","rq")
  
  out["rq",] <- c(mean((rq - mean(rq))^3/sd(rq)^3), mean((rq - mean(rq))^4/sd(rq)^4),
               mean(rq), median(rq), sd(rq), nortest::ad.test(rq)$statistic, nortest::ad.test(rq)$p.value)
  
  out["y",] <- c(mean((y - mean(y))^3/sd(y)^3), mean((y - mean(y))^4/sd(y)^4),
               mean(y), median(y), sd(y), nortest::ad.test(y)$statistic, nortest::ad.test(y)$p.value)
  
  set.seed(11112017)
  r1 <- rt(n = length(y), df = df)
  r2 <- rnorm(n = length(rq))
  par(mfrow = c(1,2))
  qqplot(y,r1, main = "Resíduo de Pearson")#, xlim = c(-4.5,4.5), ylim = c(-4.5,4.5))
  qqline(r1)
  qqplot(rq,r2, main = "Resíduo quantílico")#, xlim = c(-3,3), ylim = c(-3,3))
  qqline(r2)
  
  # output
  print(round(out,4))
  list(stats = round(out,4), resid.q = rq)
}