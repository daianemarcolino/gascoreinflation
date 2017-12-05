diag.dcs <- function(out, type = "t"){
  
  # out: saida da fç dcs_fk_estimation

  # resíduo de pearson
  ep <- out$out[,"epsilon"]
  
  if(type == "t"){
    
    # graus de liberdade
    df <- out$otimizados$par[4]
    
    Ft <- pt(ep, df = df, lower.tail = T)
    rq <- qnorm(Ft)
    
    par(mfrow = c(2,2))
    ts.plot(ep, main = "epsilon")
    ts.plot(rq, main = "Resíduo Quantílico")
    hist(ep, main = "epsilon")
    hist(rq, main = "Resíduo Quantílico")
    par(mfrow = c(1,1))
    
    x <- data.frame(matrix(NA, ncol = 9, nrow = 2))
    colnames(x) <- c("Assimetria","Curtose","Média","Mediana","Desvio-padrão","AD.stat","AD.pvalue","JB.stat","JB.pvalue")
    rownames(x) <- c("epsilon","res_quantilico")
    
    x["res_quantilico",] <- c(TSA::skewness(rq), TSA::kurtosis(rq) + 3, mean(rq), median(rq), sd(rq),
                              nortest::ad.test(rq)$statistic,  nortest::ad.test(rq)$p.value,
                              tseries::jarque.bera.test(rq)$statistic, tseries::jarque.bera.test(rq)$p.value)
    x["epsilon",] <- c(TSA::skewness(ep), TSA::kurtosis(ep) + 3, mean(ep), median(ep), sd(ep), 
                       nortest::ad.test(ep)$statistic, nortest::ad.test(ep)$p.value,
                       tseries::jarque.bera.test(ep)$statistic, tseries::jarque.bera.test(ep)$p.value)
    
    set.seed(123)
    r1 <- rt(length(ep), df = df)
    r2 <- rnorm(length(rq))
    
    par(mfrow = c(1,2))
    qqplot(y = ep, x = r1, main = "epsilon", ylab = "epsilon", xlab = "t-student, mean = 0", ylim = c(round(min(ep,r1)),round(max(ep,r1))), xlim = c(round(min(ep,r1)),round(max(ep,r1))))
    abline(a = 0, b = 1, lwd = 2)
    qqplot(y = rq, x = r2, main = "resíduo quantílico", xlab = "normal(0,1)", ylim = c(round(min(rq,r2)),round(max(rq,r2))), xlim = c(round(min(rq,r2)),round(max(rq,r2))))
    abline(a = 0, b = 1, lwd = 2)
    par(mfrow = c(1,1))
    
    # output
    list(stats = round(x,2), res_quantilico = rq)
    
  }else if(type == "norm"){
    
    par(mfrow = c(1,3))
    ts.plot(ep, main = "epsilon")
    hist(ep, main = "epsilon")
    
    
    x <- data.frame(matrix(NA, ncol = 9, nrow = 1))
    colnames(x) <- c("Assimetria","Curtose","Média","Mediana","Desvio-padrão","AD.stat","AD.pvalue","JB.stat","JB.pvalue")
    rownames(x) <- c("epsilon")
    
    x["epsilon",] <- c(TSA::skewness(ep), TSA::kurtosis(ep) + 3, mean(ep), median(ep), sd(ep), 
                       goftest::ad.test(ep, "pnorm")$statistic, goftest::ad.test(ep, "pnorm")$p.value,
                       tseries::jarque.bera.test(ep)$statistic, tseries::jarque.bera.test(ep)$p.value)
    

    r1 <- rnorm(length(ep))

    qqplot(y = ep, x = r1, main = "qqplot", xlab = "normal(0,1)", ylim = c(round(min(ep,r1)),round(max(ep,r1))), xlim = c(round(min(ep,r1)),round(max(ep,r1))))
    abline(a = 0, b = 1, lwd = 2)
    par(mfrow = c(1,1))
    
    # output
    list(stats = round(x,2))
  }

}



# if(psi){
#   if(dummy){
#     ncp <- ts(rowSums(out$out[,c("mu","gamma","psi","dummy")]), start = start(out$out), freq = 12)
#   }else{
#     ncp <- ts(rowSums(out$out[,c("mu","gamma","psi")]), start = start(out$out), freq = 12)
#   }
# }else{
#   if(dummy){
#     ncp <- ts(rowSums(out$out[,c("mu","gamma","dummy")]), start = start(out$out), freq = 12)
#   }else{
#     ncp <- ts(rowSums(out$out[,c("mu","gamma")]), start = start(out$out), freq = 12)
#   }
# }