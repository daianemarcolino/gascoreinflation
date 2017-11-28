diag.dcs <- function(y, out, psi = F, dummy = F){
  
  # y: serie y
  # out: saida da fç dcs_fk_estimation
  df <- out$otimizados$par[4]
  if(psi){
    if(dummy){
      ncp <- ts(rowSums(out$out[,c("mu","gamma","psi","dummy")]), start = start(out$out), freq = 12)
    }else{
      ncp <- ts(rowSums(out$out[,c("mu","gamma","psi")]), start = start(out$out), freq = 12)
    }
  }else{
    if(dummy){
      ncp <- ts(rowSums(out$out[,c("mu","gamma","dummy")]), start = start(out$out), freq = 12)
    }else{
      ncp <- ts(rowSums(out$out[,c("mu","gamma")]), start = start(out$out), freq = 12)
    }
   
  }
  
  
  ep <- out$out[,"epsilon"]
  
  # resíduo quantílico
  Ft <- pt(y, df = df, ncp = ncp)
  rq <- qnorm(Ft)
  
  par(mfrow = c(2,3))
  ts.plot(y, main = "y")
  ts.plot(ep, main = "epsilon")
  ts.plot(rq, main = "Resíduo Quantílico")
  hist(y, main = "y")
  hist(ep, main = "epsilon")
  hist(rq, main = "Resíduo Quantílico")
  par(mfrow = c(1,1))
  
  x <- data.frame(matrix(NA, ncol = 7, nrow = 3))
  colnames(x) <- c("Assimetria","Curtose","Média","Mediana","Desvio-padrão","AD.stat","AD.pvalue")
  rownames(x) <- c("y","epsilon","res_quantilico")
  
  x["y",] <- c(mean((y - mean(y))^3/sd(y)^3), mean((y - mean(y))^4/sd(y)^4),
               mean(y), median(y), sd(y), nortest::ad.test(y)$statistic, nortest::ad.test(y)$p.value)
  x["res_quantilico",] <- c(mean((rq - mean(rq))^3/sd(rq)^3), mean((rq - mean(rq))^4/sd(rq)^4),
                            mean(rq), median(rq), sd(rq), nortest::ad.test(rq)$statistic, nortest::ad.test(rq)$p.value)
  x["epsilon",] <- c(mean((ep - mean(ep))^3/sd(ep)^3), mean((ep - mean(ep))^4/sd(ep)^4),
                     mean(ep), median(ep), sd(ep), nortest::ad.test(ep)$statistic, nortest::ad.test(ep)$p.value)
  
  r1 <- rt(n = length(y), ncp = ncp, df = df)
  r2 <- rt(n = length(ep), df = df)
  r3 <- rnorm(n = length(rq))
  
  par(mfrow = c(1,3))
  qqplot(y = y, x = r1, main = "y", ylab = "y", xlab = "t-student, mean = mu + gamma")
  qqline(y = y, distribution = function(x) qt(x, ncp = mean(ncp), df = df), lwd = 2)
  qqplot(y = ep, x = r2, main = "epsilon", ylab = "epsilon", xlab = "t-student, mean = 0")
  qqline(y = ep, distribution = function(x) qt(x, df = df), lwd = 2)
  qqplot(y = rq, x = r3, main = "resíduo quantílico", xlab = "normal(0,1)")
  qqline(y = rq, lwd = 2)
  par(mfrow = c(1,1))
  
  # output
  list(stats = round(x,4), res_quantilico = rq)
}