dcs_fk_estimation <- function(y, initial = NULL, type = "BSM", dummy = NULL){
  
  if(type == "BSM"){ 
    
    # modelo:
    # y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
    # mu[t+1] = mu[t] + beta[t] + a1*u1[t] 
    # beta[t+1] = beta[t] + a2*u1[t] 
    # gamma[t+1] = gamma[t] + k*u[t] 
    # > estimação via ML para densidade condicional de y t-student com média constante e variância constante no tempo
    
    
    # OLS - initial gamma's
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
                   BETS::BETS.dummy(start = start(y), end = end(y), month = 11)
    )
    dados <- data.frame(cbind(y,dummy))[1:13,]
    dados <- cbind(dados,1:13)
    colnames(dados) <- c("y",paste0("D",1:11),"trend")
    m <- lm(y ~ ., data = dados)
    initial_ols <- as.vector(m$coefficients)
    
    otimizar <- function(y, par, initial_ols){
      #par = c(0.5,0.2,0.3,5)
      N <- length(y)
      k1 <- par[1]
      k2 <- k1^2/(2-k1)
      ks <- par[2]
      f2 <- par[3]
      df <- par[4]
      beta <- initial_ols[1]
      mu <- initial_ols[13]
      alpha <- matrix(NA, ncol = 12, nrow = N)
      alpha[1,] <- c(initial_ols[-c(1,13)], - sum(initial_ols[-c(1,13)]))
      j <- cycle(y)
      gamma <- alpha[1,j[1]]
      u1 <- NULL
      
      for(t in 1:(N-1)){ 
        u1[t] <- (df + 1)*((y[t] - mu[t] - gamma[t])/(df*exp(2*f2)))/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
        mu[t+1] <- mu[t] + beta[t] + k1*u1[t] 
        beta[t+1] <- beta[t] + k2*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      } 
      t <- N
      u1[t] <- (df + 1)*((y[t] - mu[t] - gamma[t])/(df*exp(2*f2)))/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      
      mu <- ts(mu, end = end(y), freq = frequency(y))
      beta <- ts(beta, end = end(y), freq = frequency(y))
      gamma <- ts(gamma, end = end(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
      -loglik
    }
    
    otimizados <- nlminb(start = c(0.3,0.1,3,4), objective = otimizar, y = y, initial_ols = initial_ols, 
                         lower = c(0,-Inf,-Inf,2), upper = c(1,Inf,Inf,Inf), control = list(eval.max = 2000, iter.max = 2000))
    
    
    N <- length(y)
    k1 <- otimizados$par[1]
    k2 <- k1^2/(2-k1)
    ks <- otimizados$par[2]
    f2 <- otimizados$par[3]
    df <- otimizados$par[4]
    beta <- initial_ols[1]
    mu <- initial_ols[13]
    alpha <- matrix(NA, ncol = 12, nrow = N)
    alpha[1,] <- c(initial_ols[-c(1,13)], - sum(initial_ols[-c(1,13)]))
    j <- cycle(y)
    gamma <- alpha[1,j[1]]
    u1 <- NULL
    
    for(t in 1:(N-1)){ 
      u1[t] <- (df + 1)*((y[t] - mu[t] - gamma[t])/(df*exp(2*f2)))/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      mu[t+1] <- mu[t] + beta[t] + k1*u1[t] 
      beta[t+1] <- beta[t] + k2*u1[t]
      k <- rep(-ks/11,12) 
      k[j[t]] <- ks
      alpha[t+1,] <- alpha[t,] + k*u1[t]
      gamma[t+1] <- alpha[t+1,j[t+1]]
    } 
    t <- N
    u1[t] <- (df + 1)*((y[t] - mu[t] - gamma[t])/(df*exp(2*f2)))/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
    
    mu <- ts(mu, end = end(y), freq = frequency(y))
    beta <- ts(beta, end = end(y), freq = frequency(y))
    gamma <- ts(gamma, end = end(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
    epsilon <- (y - mu - gamma)/exp(f2)
    out <- cbind(mu, beta, gamma, f2, exp(f2), epsilon)
    colnames(out) <- c("mu","beta","gamma","f2","sigma","epsilon")
  
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
