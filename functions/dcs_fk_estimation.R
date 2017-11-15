dcs_fk_estimation <- function(y, initial = NULL, type = "BSM1", initial.optim = T){
  
  # BSM1: mu, beta, gamma, variancia constante
  # BSM2: mu, gamma, variancia constante
  # BSM_artigo: personalizado - Caivano, Harvey e Luati (2016)
  
  if(type == "BSM1"){ 
    # BSM COM BETA -----
    # modelo:
    # y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
    # mu[t+1] = mu[t] + beta[t] + a1*u1[t] 
    # beta[t+1] = beta[t] + a2*u1[t] 
    # gamma[t+1] = gamma[t] + k*u[t] 
    # > estimação via ML para densidade condicional de y t-student com média constante e variância constante no tempo
    
    otimizar <- function(y, par, initial_beta, initial_mu, initial_gamma){
      
      N <- length(y)
      k1 <- par[1]
      k2 <- k1^2/(2-k1)
      ks <- par[2]
      f2 <- par[3]
      df <- par[4]
      if(initial.optim){
        beta <- par[5]
        mu <- par[6]
        alpha <- matrix(NA, ncol = 12, nrow = N)
        alpha[1,] <- c(par[7:17], - sum(par[7:17]))
      }else{
        beta <- initial_beta
        mu <- initial_mu
        alpha <- matrix(NA, ncol = 12, nrow = N)
        alpha[1,] <- c(initial_gamma, - sum(initial_gamma))  
      }
      j <- cycle(y)
      gamma <- alpha[1,j[1]]
      u1 <- NULL
      
      for(t in 1:(N-1)){ 
        u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
        mu[t+1] <- mu[t] + beta[t] + k1*u1[t] 
        beta[t+1] <- beta[t] + k2*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      } 
      t <- N
      u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      
      mu <- ts(mu, end = end(y), freq = frequency(y))
      beta <- ts(beta, end = end(y), freq = frequency(y))
      gamma <- ts(gamma, end = end(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
      -loglik
    }
    
    otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y, 
                         initial_beta = initial$beta, initial_mu = initial$mu, initial_gamma = initial$gamma, 
                         lower = initial$par$lower, upper = initial$par$upper, control = list(eval.max = 10000, iter.max = 10000))
    
    N <- length(y)
    k1 <- otimizados$par[1]
    k2 <- k1^2/(2-k1)
    ks <- otimizados$par[2]
    f2 <- otimizados$par[3]
    df <- otimizados$par[4]
    if(initial.optim){
      beta <- otimizados$par[5]
      mu <- otimizados$par[6]
      alpha <- matrix(NA, ncol = 12, nrow = N)
      alpha[1,] <- c(otimizados$par[7:17], - sum(otimizados$par[7:17]))
    }else{
      beta <- initial$beta
      mu <- initial$mu
      alpha <- matrix(NA, ncol = 12, nrow = N)
      alpha[1,] <- c(initial$gamma, - sum(initial$gamma))
    }
    j <- cycle(y)
    gamma <- alpha[1,j[1]]
    u1 <- NULL
    
    for(t in 1:(N-1)){ 
      u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      mu[t+1] <- mu[t] + beta[t] + k1*u1[t] 
      beta[t+1] <- beta[t] + k2*u1[t]
      k <- rep(-ks/11,12) 
      k[j[t]] <- ks
      alpha[t+1,] <- alpha[t,] + k*u1[t]
      gamma[t+1] <- alpha[t+1,j[t+1]]
    } 
    t <- N
    u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
    
    mu <- ts(mu, end = end(y), freq = frequency(y))
    beta <- ts(beta, end = end(y), freq = frequency(y))
    gamma <- ts(gamma, end = end(y), freq = frequency(y))
    u1 <- ts(u1, end = end(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
    epsilon <- (y - mu - gamma)/exp(f2)
    nu <- (y - mu - gamma)
    score <- (df + 1)/(df*exp(2*f2))*u1
    out <- cbind(mu, beta, gamma, f2, exp(f2), epsilon, nu, score, u1)
    colnames(out) <- c("mu","beta","gamma","f2","sigma","epsilon","nu","score","u")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }else if(type == "BSM2"){ 
    # BSM SEM BETA -----
    
    # modelo:
    # y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
    # mu[t+1] = mu[t] + k1*u1[t]
    # gamma[t+1] = gamma[t] + k*u[t] 
    # > estimação via ML para densidade condicional de y t-student com média constante e variância constante no tempo
    
    
    otimizar <- function(y, par, initial_mu, initial_gamma){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[2]
      f2 <- par[3]
      df <- par[4]
      if(initial.optim){
        mu <- par[5]
        alpha <- matrix(NA, ncol = 12, nrow = N)
        alpha[1,] <- c(par[6:16], - sum(par[6:16]))
      }else{
        mu <- initial_mu
        alpha <- matrix(NA, ncol = 12, nrow = N)
        alpha[1,] <- c(initial_gamma, - sum(initial_gamma))  
      }
      j <- cycle(y)
      gamma <- alpha[1,j[1]]
      u1 <- NULL
      
      for(t in 1:(N-1)){ 
        u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
        mu[t+1] <- mu[t] + k1*u1[t] 
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      } 
      t <- N
      u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      
      mu <- ts(mu, end = end(y), freq = frequency(y))
      gamma <- ts(gamma, end = end(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
      -loglik
    }
    
    otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y, 
                         initial_mu = initial$mu, initial_gamma = initial$gamma, 
                         lower = initial$par$lower, upper = initial$par$upper, control = list(eval.max = 10000, iter.max = 10000))
    
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[2]
    f2 <- otimizados$par[3]
    df <- otimizados$par[4]
    if(initial.optim){
      mu <- otimizados$par[5]
      alpha <- matrix(NA, ncol = 12, nrow = N)
      alpha[1,] <- c(otimizados$par[6:16], - sum(otimizados$par[6:16]))
    }else{
      mu <- initial$mu
      alpha <- matrix(NA, ncol = 12, nrow = N)
      alpha[1,] <- c(initial$gamma, - sum(initial$gamma))
    }
    j <- cycle(y)
    gamma <- alpha[1,j[1]]
    u1 <- NULL
    
    for(t in 1:(N-1)){ 
      u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      mu[t+1] <- mu[t] + k1*u1[t]
      k <- rep(-ks/11,12) 
      k[j[t]] <- ks
      alpha[t+1,] <- alpha[t,] + k*u1[t]
      gamma[t+1] <- alpha[t+1,j[t+1]]
    } 
    t <- N
    u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
    
    mu <- ts(mu, end = end(y), freq = frequency(y))
    gamma <- ts(gamma, end = end(y), freq = frequency(y))
    u1 <- ts(u1, end = end(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
    epsilon <- (y - mu - gamma)/exp(f2)
    nu <- (y - mu - gamma)
    score <- ((df + 1)/(df*exp(2*f2)))*u1
    out <- cbind(mu, gamma, f2, exp(f2), epsilon, nu, score, u1)
    colnames(out) <- c("mu","gamma","f2","sigma","epsilon", "nu", "score","u")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }else if(type == "BSM_artigo"){ 
    # ARTIGO -----  
    otimizar <- function(y, par){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[2]
      f2 <- par[3]
      df <- par[4]
      mu <- par[5]
      alpha <- matrix(NA, ncol = 12, nrow = N)
      alpha[1,] <- c(par[6:16], - sum(par[6:16]))
      j <- cycle(y)
      gamma <- alpha[1,j[1]]
      u1 <- NULL
      
      for(t in 1:(N-1)){ 
        u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      } 
      t <- N
      u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      
      mu <- ts(mu, end = end(y), freq = frequency(y))
      gamma <- ts(gamma, end = end(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
      -loglik
    }
    
    initial <- data.frame(param = c("k1","ks","f2","df","mu"   ,"g1"   ,"g2"   ,"g3"   ,"g4"  ,"g5"  ,"g6"  ,"g7"  ,"g8"  ,"g9"  ,"g10" ,"g11"),
                          value = c(0.5 ,0.5     ,0.5 ,3   ,15.1018,-0.5795,-0.4761,-0.1780,0.1733,0.1475,0.2648,0.5791,0.4566,0.3157,0.1167,-0.3862),
                          lower = c(0,-Inf,-Inf,2  ,-Inf, rep(-Inf,11)),
                          upper = c(1, Inf, Inf,Inf, Inf, rep(Inf,11)))
    
    
    otimizados <- nlminb(start = initial$value, objective = otimizar, y = y,
                         lower = initial$lower, upper = initial$upper, control = list(eval.max = 10000, iter.max = 10000))
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[2]
    f2 <- otimizados$par[3]
    df <- otimizados$par[4]
    mu <- otimizados$par[5]
    alpha <- matrix(NA, ncol = 12, nrow = N)
    alpha[1,] <- c(otimizados$par[6:16], - sum(otimizados$par[6:16]))
    j <- cycle(y)
    gamma <- alpha[1,j[1]]
    u1 <- NULL
    
    for(t in 1:(N-1)){ 
      u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      mu[t+1] <- mu[t] + k1*u1[t]
      k <- rep(-ks/11,12) 
      k[j[t]] <- ks
      alpha[t+1,] <- alpha[t,] + k*u1[t]
      gamma[t+1] <- alpha[t+1,j[t+1]]
    } 
    t <- N
    u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
    
    mu <- ts(mu, end = end(y), freq = frequency(y))
    gamma <- ts(gamma, end = end(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
    epsilon <- (y - mu - gamma)/exp(f2)
    nu <- (y - mu - gamma)
    u1 <- ts(u1, end = end(y), freq = frequency(y))
    score <- (df + 1)/(df*exp(2*f2))*u1
    out <- cbind(mu, gamma, f2, exp(f2), epsilon, nu, score, u1)
    colnames(out) <- c("mu","gamma","f2","sigma","epsilon","nu","score","u")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }
  
  # # smoothing 
  # pseudo.y <- out[,"mu"] + out[,"gamma"] + out[,"score"]
  # 
  # 
  # fk <- StructTS(pseudo.y, type = "BSM")
  # ts.plot(pseudo.y,fk$fitted[,"slope"])
  
}

# # OLS - initial gamma's
# dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 11)
# )
# dados <- data.frame(cbind(y,dummy))[1:13,]
# dados <- cbind(dados,1:13)
# colnames(dados) <- c("y",paste0("D",1:11),"trend")
# m <- lm(y ~ ., data = dados)
# initial_ols <- as.vector(m$coefficients)
# otimo <- data.frame(param = c("k1","ks","f2","df","mu"   ,"g1"   ,"g2"   ,"g3"   ,"g4"  ,"g5"  ,"g6"  ,"g7"  ,"g8"  ,"g9"  ,"g10" ,"g11"),
#                      value = otimizados$par)
# otimo
# modelo:
# y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
# mu[t+1] = mu[t] + k1*u1[t] 
# gamma[t+1] = gamma[t] + k*u[t] 
# > estimação via ML para densidade condicional de y t-student com média constante e variância constante no tempo


# OLS - initial gamma's
# dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 11)
# )
# 
# # dummies de outliers
# outliers <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 3, year = 2002),
#                   BETS::BETS.dummy(start = start(y), end = end(y), month = 4, year = 2002),
#                   BETS::BETS.dummy(start = start(y), end = end(y), month = 8, year = 2002),
#                   BETS::BETS.dummy(start = start(y), end = end(y), month = 4, year = 2005),
#                   BETS::BETS.dummy(start = start(y), end = end(y), month = 4, year = 2010),
#                   BETS::BETS.dummy(start = start(y), end = end(y), to = c(2008,4))#, to = end(y))
# )
# 
# dados <- data.frame(cbind(y,dummy,outliers))
# colnames(dados) <- c("y",paste0("D",1:11),paste0("outliers_",1:5),"level_break")
# m <- lm(y ~ ., data = dados)
# initial_ols <- as.vector(m$coefficients)
# ts.plot(y,out[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência + Sazonalidade")
# ts.plot(y,out[,"f1"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
# ts.plot(y,out[,"f2"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
# ts.plot(y,out[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")
# plot(out[,c("f1","f2")], col = 1:2, lwd = c(1,2), main = "Tendência e Sazonalidade")
