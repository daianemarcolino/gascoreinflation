betategarch_forecasting <- function(out, y, type = "mean-var6", dummy = NULL, h = 12){
  # out: saída da função betategarch_estimation
  # y: série a ser prevista
  # type: modelo
 
  if(type == "mean-var6"){ # modelo y[t] = f1[t] + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v) c/ dummy
    par <- out$otimizados$par

    w1 <- par[1]
    w2 <- par[2]
    A1_00 <- par[3]
    A1_01 <- par[4]
    A1_11 <- par[5]
    B1_00 <- par[6]
    B1_01 <- par[7]
    A2_00 <- par[8]
    B2_00 <- par[9]
    df <- par[10]
    # D1 <- par[11]
    # D2 <- par[12]
    # D3 <- par[13]
    # D4 <- par[14]
    # dummy <- ts(eval(parse(text=paste0("rbind(dummy,", paste0(rep(0,h), collapse = ","), ")"))), start = start(dummy), freq = frequency(dummy))
    
    # séries dentro da amostra
    f1 <- c(out$out[,"f1"], rep(NA,h))
    f2 <- c(out$out[,"f2"], rep(NA,h))
    u1 <- c(out$out[,"u1"], rep(NA,h))
    u2 <- c(out$out[,"u2"], rep(NA,h))
    y <- c(window(y, start = start(out$out), freq = frequency(y)), rep(NA,h))
    N <- nrow(out$out)
    
    # cálculo das séries pro futuro da amostra
    for(t in 1:h){
      u1[N+t-1] <- (df + 1)*((y[N+t-1] - f1[N+t-1])/(df*exp(2*f2[N+t-1])))/(1 + (y[N+t-1] - f1[N+t-1])^2/(df*exp(2*f2[N+t-1])))
      u2[N+t-1] <- (((df + 1)*(y[N+t-1] - f1[N+t-1])^2) / (df*exp(2*f2[N+t-1]) + (y[N+t-1] - f1[N+t-1])^2)) - 1
      f1[N+t] <- w1 + A1_00*u1[N+t-1] + A1_01*u1[N+t-2] + A1_11*u1[N+t-12] + B1_00*f1[N+t-1] + B1_01*f1[N+t-2] # + D1*dummy[N+t,1] + D2*dummy[N+t,2] + D3*dummy[N+t,3] + D4*dummy[N+t,4]
      f2[N+t] <- w2 + A2_00*u2[N+t-1] + B2_00*f2[N+t-1]
      y[N+t] <- f1[N+t]
    }
    
    f1 <- ts(f1, start = start(out$out), freq = frequency(out$out))
    f2 <- ts(f2, start = start(out$out), freq = frequency(out$out))
    u1 <- ts(u1, start = start(out$out), freq = frequency(out$out))
    u2 <- ts(u2, start = start(out$out), freq = frequency(out$out))
    y <- ts(y, start = start(out$out), freq = frequency(out$out))
    
    output <- cbind(y,f1,f2,exp(f2),u1,u2)
    colnames(output) <- c("ychapeu", "f1", "f2", "sigma", "u1", "u2")
  
    # output
    tail(output,h)
    }
}