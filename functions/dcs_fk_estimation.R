dcs_fk_estimation <- function(y, initial = NULL, type = "mean", dummy = NULL){

 if(type == "mean"){ # modelo y[t] = f1[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
    # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
    otimizar <- function(y, par){
      w1 <- par[1]
      A1 <- par[2]
      B1 <- par[3]
      df <- par[4]
      f2 <- par[5]
      f1 <- w1/(1-B1)
      N <- length(y)
      u1 <- NULL
      for(t in 1:N){
        u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2)))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2)))
        f1[t+1] <- w1 + A1*u1[t] + B1*f1[t]
      }
      
      f1 <- ts(f1[-(N+1)], start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
      -loglik
    }
    
    
    otimizados <- nlminb(start = c(initial), objective = otimizar, y = y, 
                         lower = c(-Inf,-Inf,-1,2,-Inf), upper = c(Inf,Inf,1,Inf,Inf), control = list(eval.max = 10000, iter.max = 10000))
    
    w1 <- otimizados$par[1]
    A1 <- otimizados$par[2]
    B1 <- otimizados$par[3]
    df <- otimizados$par[4]
    f2 <- otimizados$par[5]
    f1 <- w1/(1-B1)
    N <- length(y)
    u1 <- NULL
    
    for(t in 1:N){
      u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2)))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2)))
      f1[t+1] <- w1 + A1*u1[t] + B1*f1[t]
    }
    
    f1 <- ts(f1[-(N+1)], start = start(y), freq = frequency(y))
    u1 <- ts(c(u1[1:(N-1)], NA), end = end(y), freq =  frequency(y))
    
    
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
    epsilon <- (y - f1)/exp(f2)
    out <- cbind(f1, f2, exp(f2), epsilon, u1)
    colnames(out) <- c("f1", "f2", "sigma","epsilon","u1")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
 }else if(type == "cn"){ # modelo y[t] = f1[t] + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
    # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
    otimizar <- function(y, par){
      w1 <- par[1]
      A1 <- par[2]
      df <- par[3]
      f2 <- par[4]
      f1 <- c(par[5],y[1])
      N <- length(y)
      u1 <- NULL
      for(t in 1:N){
        u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2)))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2)))
        f1[t+1] <-  w1 + A1*u1[t] + f1[t]
      }
      
      f1 <- ts(f1[-(N+1)], start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
      -loglik
    }
    
    
    otimizados <- nlminb(start = c(initial), objective = otimizar, y = y, 
                         lower = c(-Inf,-Inf,2,-Inf,-Inf), upper = c(Inf,Inf,Inf,Inf,Inf), control = list(eval.max = 10000, iter.max = 10000))
    
    w1 <- otimizados$par[1]
    A1 <- otimizados$par[2]
    df <- otimizados$par[3]
    f2 <- otimizados$par[4]
    f1 <- c(otimizados$par[5],y[1])
    N <- length(y)
    u1 <- NULL
    
    for(t in 1:N){
      u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2)))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2)))
      f1[t+1] <- w1 + A1*u1[t] + f1[t]
    }
    
    f1 <- ts(f1[-(N+1)], start = start(y), freq = frequency(y))
    u1 <- ts(c(u1[1:(N-1)], NA), end = end(y), freq =  frequency(y))
    
    
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
    epsilon <- (y - f1)/exp(f2)
    out <- cbind(f1, f2, exp(f2), epsilon, u1)
    colnames(out) <- c("f1", "f2", "sigma","epsilon","u1")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
 }else if(type == "cn2"){ # modelo y[t] = f1[t] + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
   # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
   
   dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 11),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 12)
   )
   
   # initial dummies
   
   y_filter <- y[1:13]
   dummy_filter <- dummy[1:13,]
   t <- 1:13
   dados <- data.frame(cbind(y_filter,t,dummy_filter))
   colnames(dados) <- c("y","t",paste0("D",1:12))
   m <- lm(y ~ ., data = dados)
   
   otimizar <- function(y, par, dummy){
     gamma1 <- par[1]
     gamma2 <- par[2]
     gamma3 <- par[3]
     gamma4 <- par[4]
     gamma5 <- par[5]
     gamma6 <- par[6]
     gamma7 <- par[7]
     gamma8 <- par[8]
     gamma9 <- par[9]
     gamma10 <- par[10]
     gamma11 <- par[11]
     gamma12 <- 1 - (gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + gamma7 + gamma8 + gamma9 + 
                       gamma10 + gamma11)
     df <- par[12]
     f2 <- par[13]
     mu <- par[14]
     beta <- par[15]
     k1 <- par[16]
     k2 <- par[17]
     
     N <- length(y)
     u1 <- NULL
     gamma0 <- NULL
     
     for(t in 1:N){
       gamma0[t] <- gamma1*dummy[t,1] + gamma2*dummy[t,2] + gamma3*dummy[t,3] + gamma4*dummy[t,4] +
         dummy[t,5] + gamma6*dummy[t,6] + gamma7*dummy[t,7] + gamma8*dummy[t,8] + 
         gamma9*dummy[t,9] + gamma10*dummy[t,10] + gamma11*dummy[t,11] + gamma12*dummy[t,12]
       u1[t] <- (df + 1)*((y[t] - mu[t] - gamma0[t])/(df*exp(2*f2)))/(1 + (y[t] - mu[t] - gamma0[t])^2/(df*exp(2*f2)))
       beta[t+1] <- beta[t] + k2*u1[t]
       mu[t+1] <- mu[t] + beta[t] + k1*u1[t] 
       
     }
     
     mu <- ts(mu[-(N+1)], start = start(y), freq = frequency(y))
     gamma0 <- ts(gamma0[-(N+1)], start = start(y), freq = frequency(y))
     beta <- ts(beta[-(N+1)], start = start(y), freq = frequency(y))
     
     # loglik
     loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma0)^2/(df*exp(2*f2))))
     -loglik
   }
   
   otimizados <- nlminb(start = c(as.vector(m$coefficients[3:13]), 10,0.5,0.5,0.5,0.4,0.5), objective = otimizar, y = y, dummy = dummy,
                        lower = c(rep(-Inf,11),2,-Inf,-Inf,-Inf,0,0), upper = c(rep(Inf,11),Inf,Inf,Inf,Inf,1,1), control = list(eval.max = 10000, iter.max = 10000))
   otimizados
   gamma1 <- otimizados$par[1]
   gamma2 <- otimizados$par[2]
   gamma3 <- otimizados$par[3]
   gamma4 <- otimizados$par[4]
   gamma5 <- otimizados$par[5]
   gamma6 <- otimizados$par[6]
   gamma7 <- otimizados$par[7]
   gamma8 <- otimizados$par[8]
   gamma9 <- otimizados$par[9]
   gamma10 <- otimizados$par[10]
   gamma11 <- otimizados$par[11]
   gamma12 <- 1 - (gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + gamma7 + gamma8 + gamma9 + 
                     gamma10 + gamma11)
   df <- otimizados$par[12]
   f2 <- otimizados$par[13]
   mu <- otimizados$par[14]
   beta <- otimizados$par[15]
   k1 <- otimizados$par[16]
   k2 <- otimizados$par[17]
   
   N <- length(y)
   u1 <- NULL
   gamma0 <- NULL
   
   for(t in 1:N){
     gamma0[t] <- gamma1*dummy[t,1] + gamma2*dummy[t,2] + gamma3*dummy[t,3] + gamma4*dummy[t,4] +
       dummy[t,5] + gamma6*dummy[t,6] + gamma7*dummy[t,7] + gamma8*dummy[t,8] + 
       gamma9*dummy[t,9] + gamma10*dummy[t,10] + gamma11*dummy[t,11] + gamma12*dummy[t,12]
     u1[t] <- (df + 1)*((y[t] - mu[t] - gamma0[t])/(df*exp(2*f2)))/(1 + (y[t] - mu[t] - gamma0[t])^2/(df*exp(2*f2)))
     beta[t+1] <- beta[t] + k2*u1[t]
     mu[t+1] <- mu[t] + beta[t] + k1*u1[t] 
   }
   
   mu <- ts(mu[-(N+1)], start = start(y), freq = frequency(y))
   gamma0 <- ts(gamma0[-(N+1)], start = start(y), freq = frequency(y))
   beta <- ts(beta[-(N+1)], start = start(y), freq = frequency(y))
   u1 <- ts(u1, end = end(y), freq =  frequency(y))
   loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma0)^2/(df*exp(2*f2))))
   epsilon <- (y - mu - gamma0)/exp(f2)
   out <- cbind(mu, beta, gamma0, exp(f2), epsilon, u1)
   colnames(out) <- c("mu", "beta", "gamma", "sigma","epsilon","u1")
   
   # output
   print(otimizados)
   invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
   
 }else if(type == "cn3"){ # modelo y[t] = f1[t] + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
   # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
   
   dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 11),
                  BETS::BETS.dummy(start = start(y), end = end(y), month = 12)
   )
   
   # initial dummies
   
   y_filter <- y[1:13]
   dummy_filter <- dummy[1:13,]
   t <- 1:13
   dados <- data.frame(cbind(y_filter,t,dummy_filter))
   colnames(dados) <- c("y","t",paste0("D",1:12))
   m <- lm(y ~ ., data = dados)
   
   otimizar <- function(y, par, dummy){
     gamma1 <- par[1]
     gamma2 <- par[2]
     gamma3 <- par[3]
     gamma4 <- par[4]
     gamma5 <- par[5]
     gamma6 <- par[6]
     gamma7 <- par[7]
     gamma8 <- par[8]
     gamma9 <- par[9]
     gamma10 <- par[10]
     gamma11 <- par[11]
     gamma12 <- 1 - (gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + gamma7 + gamma8 + gamma9 + 
                       gamma10 + gamma11)
     df <- par[12]
     w1 <- par[13]
     w2 <- par[14]
     f1 <- par[15]
     f2 <- par[16]
     f3 <- par[17]
     beta <- par[18]
     a11 <- par[19]
     a12 <- par[20]
     a13 <- par[21]
     a2 <- par[22]
     b2 <- par[23]
     
     N <- length(y)
     mu = u1 = u2 <- NULL
     
     for(t in 1:N){
       
       mu[t] <- w1 + f1[t] + f2[t]
       
       u1[t] <- (df + 1)*((y[t] - mu[t])/(df*exp(2*f3[t])))/(1 + (y[t] - mu[t])^2/(df*exp(2*f3[t])))
       u2[t] <- (((df + 1)*(y[t] - mu[t])^2) / (df*exp(2*f3[t]) + (y[t] - mu[t])^2)) - 1
       
       f1[t+1] <- f1[t] + beta[t] + a11*u1[t]
       beta[t+1] <- beta[t] + a12*u1[t]
       if(t != N){
         f2[t+1] <- gamma1*dummy[t+1,1] + gamma2*dummy[t+1,2] + gamma3*dummy[t+1,3] + gamma4*dummy[t+1,4] +
           dummy[t+1,5] + gamma6*dummy[t+1,6] + gamma7*dummy[t+1,7] + gamma8*dummy[t+1,8] + 
           gamma9*dummy[t+1,9] + gamma10*dummy[t+1,10] + gamma11*dummy[t+1,11] + gamma12*dummy[t+1,12] +
           a13*u1[t]
       }
       f3[t+1] <- w2 + a2*u2[t] + b2*f3[t]
     }
     
     mu <- ts(mu, start = start(y), freq = frequency(y))
     f1 <- ts(f1[-(N+1)], start = start(y), freq = frequency(y))
     f2 <- ts(f2, start = start(y), freq = frequency(y))
     beta <- ts(beta[-(N+1)], start = start(y), freq = frequency(y))
     f3 <- ts(f3[-(N+1)], start = start(y), freq = frequency(y))
     
     # loglik
     loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f3) - ((df + 1)/2)*sum(log(1 + (y - mu)^2/(df*exp(2*f3))))
     -loglik
   }

   data <- data.frame(param = c("df","w1","w2","f1","f2","f3","beta","a11","a12","a13","a2","b2"),
                      value = c(   3,   1,   1,   2, 0.5, 0.5,  0.01,  0.5,  0.5,  0.5, 0.5,0.7),
                      lower = c(   2,-Inf,-Inf,-Inf,-Inf,-Inf,  -Inf, -Inf, -Inf, -Inf,-Inf,-1),
                      upper = c( Inf, Inf, Inf, Inf, Inf, Inf,   Inf,  Inf,  Inf,  Inf, Inf, 1))   
   otimizados <- nlminb(start = c(as.vector(m$coefficients[3:13]), data$value), objective = otimizar, y = y, dummy = dummy,
                        lower = c(rep(-Inf,11), data$lower), upper = c(rep(Inf,11), data$upper), control = list(eval.max = 10000, iter.max = 10000))
   #otimizados
   gamma1 <- otimizados$par[1]
   gamma2 <- otimizados$par[2]
   gamma3 <- otimizados$par[3]
   gamma4 <- otimizados$par[4]
   gamma5 <- otimizados$par[5]
   gamma6 <- otimizados$par[6]
   gamma7 <- otimizados$par[7]
   gamma8 <- otimizados$par[8]
   gamma9 <- otimizados$par[9]
   gamma10 <- otimizados$par[10]
   gamma11 <- otimizados$par[11]
   gamma12 <- 1 - (gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6 + gamma7 + gamma8 + gamma9 + 
                     gamma10 + gamma11)
   df <- otimizados$par[12]
   w1 <- otimizados$par[13]
   w2 <- otimizados$par[14]
   f1 <- otimizados$par[15]
   f2 <- otimizados$par[16]
   f3 <- otimizados$par[17]
   beta <- otimizados$par[18]
   a11 <- otimizados$par[19]
   a12 <- otimizados$par[20]
   a13 <- otimizados$par[21]
   a2 <- otimizados$par[22]
   b2 <- otimizados$par[23]
   
   N <- length(y)
   mu = u1 = u2 <- NULL
   
   for(t in 1:N){
     
     mu[t] <- w1 + f1[t] + f2[t]
     
     u1[t] <- (df + 1)*((y[t] - mu[t])/(df*exp(2*f3[t])))/(1 + (y[t] - mu[t])^2/(df*exp(2*f3[t])))
     u2[t] <- (((df + 1)*(y[t] - mu[t])^2) / (df*exp(2*f3[t]) + (y[t] - mu[t])^2)) - 1
     
     f1[t+1] <- f1[t] + beta[t] + a11*u1[t]
     beta[t+1] <- beta[t] + a12*u1[t]
     if(t != N){
       f2[t+1] <- gamma1*dummy[t+1,1] + gamma2*dummy[t+1,2] + gamma3*dummy[t+1,3] + gamma4*dummy[t+1,4] +
         dummy[t+1,5] + gamma6*dummy[t+1,6] + gamma7*dummy[t+1,7] + gamma8*dummy[t+1,8] + 
         gamma9*dummy[t+1,9] + gamma10*dummy[t+1,10] + gamma11*dummy[t+1,11] + gamma12*dummy[t+1,12] +
         a13*u1[t]
     }
     f3[t+1] <- w2 + a2*u2[t] + b2*f3[t]
   }
   
   mu <- ts(mu, start = start(y), freq = frequency(y))
   f1 <- ts(f1[-(N+1)], start = start(y), freq = frequency(y))
   f2 <- ts(f2, start = start(y), freq = frequency(y))
   beta <- ts(beta[-(N+1)], start = start(y), freq = frequency(y))
   f3 <- ts(f3[-(N+1)], start = start(y), freq = frequency(y))
   loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f3) - ((df + 1)/2)*sum(log(1 + (y - mu)^2/(df*exp(2*f3))))
   epsilon <- (y - mu)/exp(f3)
   out <- cbind(mu, f1, f2, beta, f3, exp(f3), epsilon)
   colnames(out) <- c("mu","f1","f2","beta", "f3", "sigma","epsilon")
   
   # output
   print(otimizados)
   invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
   
 }
  
}

# ts.plot(y,out[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência + Sazonalidade")
# ts.plot(y,out[,"f1"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
# ts.plot(y,out[,"f2"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
# ts.plot(y,out[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")
# plot(out[,c("f1","f2")], col = 1:2, lwd = c(1,2), main = "Tendência e Sazonalidade")
