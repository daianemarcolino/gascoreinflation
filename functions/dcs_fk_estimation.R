dcs_fk_estimation <- function(y, initial = NULL, type = "BSM1", outlier = F, otimo = T){
  
  # BSM1: mu, beta, gamma, variancia constante
  # BSM2: mu, gamma, variancia constante
  # BSM_artigo: personalizado - Caivano, Harvey e Luati (2016)
  
  if(type == "BSM_artigo"){ 
    # ARTIGO -----  
    otimizar <- function(y, par, par2){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[2]
      f2 <- par[3]
      df <- par[4]
      mu <- par2[1]
      alpha <- matrix(NA, ncol = 12, nrow = N)
      alpha[1,] <- c(par2[2:12], - sum(par2[2:12]))
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
    
    
    otimizados <- nlminb(start = initial$value[1:4], objective = otimizar, y = y, par2 = initial$value[-c(1:4)],
                         lower = initial$lower[1:4], upper = initial$upper[1:4], control = list(eval.max = 10000, iter.max = 10000))
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[2]
    f2 <- otimizados$par[3]
    df <- otimizados$par[4]
    mu <- initial$value[5]
    alpha <- matrix(NA, ncol = 12, nrow = N)
    alpha[1,] <- c(initial$value[6:16], - sum(initial$value[6:16]))
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
    b <- ((y - mu - gamma)^2/(df*exp(2*f2)))/(1 + (y - mu - gamma)^2/(df*exp(2*f2)))
    out <- cbind(mu, gamma, f2, exp(f2), epsilon, nu, score, u1, b)
    colnames(out) <- c("mu","gamma","f2","sigma","epsilon","nu","score","u","b")
    
    # output
    #print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }else if(type == "BSM1"){ 
    # DCS COM BETA VARIANTE NO TEMPO (t) -----
    # modelo:
    # y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
    # mu[t+1] = mu[t] + beta[t] + k1*u1[t] 
    # beta[t+1] = beta[t] + k2*u1[t] 
    # gamma[t+1] = gamma[t] + ks*u[t] 
    # > estimação via ML para densidade condicional de y t-student com variância constante no tempo
    
    otimizar <- function(y, par){
      
      N <- length(y)
      k1 <- par[1]
      k2 <- par[2]
      ks <- par[3]
      f2 <- par[4]
      df <- par[5]
      beta <- par[6]
      mu <- par[7]
      alpha <- matrix(NA, ncol = 12, nrow = N + 1)
      alpha[1,] <- c(par[8:18], - sum(par[8:18]))
      j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
      gamma <- alpha[1,j[1]]
      u1 <- NULL
      
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
        mu[t+1] <- mu[t] + beta[t] + k1*u1[t] 
        beta[t+1] <- beta[t] + k2*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      } 
      
      # transformar componentes em série temporal
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(beta, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
      -loglik
    }
    
    if(otimo == T){
      otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y,
                           lower = initial$par$lower, upper = initial$par$upper, control = list(eval.max = 10000, iter.max = 10000))
    }else{
      otimizados <- list()
      otimizados$par <- initial$par$value
    }
    
    N <- length(y)
    k1 <- otimizados$par[1]
    k2 <- otimizados$par[2]
    ks <- otimizados$par[3]
    f2 <- otimizados$par[4]
    df <- otimizados$par[5]
    beta <- otimizados$par[6]
    mu <- otimizados$par[7]
    alpha <- matrix(NA, ncol = 12, nrow = N + 1)
    alpha[1,] <- c(otimizados$par[8:18], - sum(otimizados$par[8:18]))
    j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
    gamma <- alpha[1,j[1]]
    u1 <- NULL
    
    for(t in 1:N){ 
      u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      mu[t+1] <- mu[t] + beta[t] + k1*u1[t] 
      beta[t+1] <- beta[t] + k2*u1[t]
      k <- rep(-ks/11,12) 
      k[j[t]] <- ks
      alpha[t+1,] <- alpha[t,] + k*u1[t]
      gamma[t+1] <- alpha[t+1,j[t+1]]
    } 
    
    # transformar componentes em série temporal
    mu <- ts(mu, start = start(y), freq = frequency(y))
    beta <- ts(beta, start = start(y), freq = frequency(y))
    gamma <- ts(gamma, start = start(y), freq = frequency(y))
    u1 <- ts(u1, end = end(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
    epsilon <- (y - mu - gamma)/exp(f2)
    nu <- (y - mu - gamma)
    score <- (df + 1)/(df*exp(2*f2))*u1
    b <- ((y - mu - gamma)^2/(df*exp(2*f2)))/(1 + (y - mu - gamma)^2/(df*exp(2*f2)))
    out <- cbind(mu, beta, gamma, f2, exp(f2), epsilon, nu, score, u1, b)
    colnames(out) <- c("mu","beta","gamma","f2","sigma","epsilon","nu","score","u","b")
    
    # output
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }else if(type == "BSM2"){ 
    # DCS SEM BETA (t) -----
    # modelo:
    # y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
    # mu[t+1] = mu[t] + k1*u1[t] 
    # gamma[t+1] = gamma[t] + ks*u[t] 
    # > estimação via ML para densidade condicional de y t-student com variância constante no tempo
    
    
    otimizar <- function(y, par){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[3]
      f2 <- par[4]
      df <- par[5]
      mu <- par[7]
      alpha <- matrix(NA, ncol = 12, nrow = N + 1)
      alpha[1,] <- c(par[8:18], - sum(par[8:18]))
      j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
      gamma <- alpha[1,j[1]]
      u1 <- NULL
      
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      } 
      
      # transformar componentes em série temporal
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(NA, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
      -loglik
    }
    
    if(otimo == T){
      otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y,
                           lower = initial$par$lower, upper = initial$par$upper, control = list(eval.max = 10000, iter.max = 10000))
    }else{
      otimizados <- list()
      otimizados$par <- initial$par$value
    }
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[3]
    f2 <- otimizados$par[4]
    df <- otimizados$par[5]
    mu <- otimizados$par[7]
    alpha <- matrix(NA, ncol = 12, nrow = N + 1)
    alpha[1,] <- c(otimizados$par[8:18], - sum(otimizados$par[8:18]))
    j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
    gamma <- alpha[1,j[1]]
    u1 <- NULL
    
    for(t in 1:N){ 
      u1[t] <- (y[t] - mu[t] - gamma[t])/(1 + (y[t] - mu[t] - gamma[t])^2/(df*exp(2*f2)))
      mu[t+1] <- mu[t] + k1*u1[t]
      k <- rep(-ks/11,12) 
      k[j[t]] <- ks
      alpha[t+1,] <- alpha[t,] + k*u1[t]
      gamma[t+1] <- alpha[t+1,j[t+1]]
    } 
    
    # transformar componentes em série temporal
    mu <- ts(mu, start = start(y), freq = frequency(y))
    beta <- ts(NA, start = start(y), freq = frequency(y))
    gamma <- ts(gamma, start = start(y), freq = frequency(y))
    u1 <- ts(u1, end = end(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma)^2/(df*exp(2*f2))))
    epsilon <- (y - mu - gamma)/exp(f2)
    nu <- (y - mu - gamma)
    score <- (df + 1)/(df*exp(2*f2))*u1
    b <- ((y - mu - gamma)^2/(df*exp(2*f2)))/(1 + (y - mu - gamma)^2/(df*exp(2*f2)))
    out <- cbind(mu, beta, gamma, f2, exp(f2), epsilon, nu, score, u1, b)
    colnames(out) <- c("mu","beta","gamma","f2","sigma","epsilon","nu","score","u","b")
    
    # output
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }else if(type == "BSM3"){ 
    # DCS SEM BETA COM PSI (t) -----
    # modelo:
    # y[t] = mu[t] + gamma[t] + psi[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
    # mu[t+1] = mu[t] + k1*u1[t] 
    # gamma[t+1] = gamma[t] + ks*u[t] 
    # psi[t+1] = phi*psi[t] + k3*u1[t]  
    # > estimação via ML para densidade condicional de y t-student com variância constante no tempo
    
    
    otimizar <- function(y, par){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[3]
      f2 <- par[4]
      df <- par[5]
      mu <- par[7]
      alpha <- matrix(NA, ncol = 12, nrow = N + 1)
      alpha[1,] <- c(par[8:18], - sum(par[8:18]))
      j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
      gamma <- alpha[1,j[1]]
      psi <- par[19]
      phi <- par[20]
      k3 <- par[21]
      u1 <- NULL
      
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])/(1 + (y[t] - mu[t] - gamma[t] - psi[t])^2/(df*exp(2*f2)))
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
        psi[t+1] <- phi*psi[t] + k3*u1[t]
      } 
      
      # transformar componentes em série temporal
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(NA, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      psi <- ts(psi, start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma - psi)^2/(df*exp(2*f2))))
      -loglik
    }
    
    if(otimo == T){
      otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y,
                           lower = initial$par$lower, upper = initial$par$upper, control = list(eval.max = 10000, iter.max = 10000))
    }else{
      otimizados <- list()
      otimizados$par <- initial$par$value
    }
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[3]
    f2 <- otimizados$par[4]
    df <- otimizados$par[5]
    mu <- otimizados$par[7]
    alpha <- matrix(NA, ncol = 12, nrow = N + 1)
    alpha[1,] <- c(otimizados$par[8:18], - sum(otimizados$par[8:18]))
    j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
    gamma <- alpha[1,j[1]]
    psi <- otimizados$par[19]
    phi <- otimizados$par[20]
    k3 <- otimizados$par[21]
    u1 <- NULL
    
    for(t in 1:N){ 
      u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])/(1 + (y[t] - mu[t] - gamma[t] - psi[t])^2/(df*exp(2*f2)))
      mu[t+1] <- mu[t] + k1*u1[t]
      k <- rep(-ks/11,12) 
      k[j[t]] <- ks
      alpha[t+1,] <- alpha[t,] + k*u1[t]
      gamma[t+1] <- alpha[t+1,j[t+1]]
      psi[t+1] <- phi*psi[t] + k3*u1[t]
    } 
    
    # transformar componentes em série temporal
    mu <- ts(mu, start = start(y), freq = frequency(y))
    beta <- ts(NA, start = start(y), freq = frequency(y))
    gamma <- ts(gamma, start = start(y), freq = frequency(y))
    psi <- ts(psi, start = start(y), freq = frequency(y))

    u1 <- ts(u1, end = end(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma - psi)^2/(df*exp(2*f2))))
    epsilon <- (y - mu - gamma - psi)/exp(f2)
    nu <- (y - mu - gamma - psi)
    b <- ((y - mu - gamma - psi)^2/(df*exp(2*f2)))/(1 + (y - mu - gamma - psi)^2/(df*exp(2*f2)))
    score <- ((df + 1)/(df*exp(2*f2)))*u1
    out <- cbind(mu, beta, gamma, psi, f2, exp(f2), epsilon, nu, score, u1, b)
    colnames(out) <- c("mu","beta","gamma","psi","f2","sigma","epsilon","nu","score","u","b")
    
    # output
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }else if(type == "BSM1_normal"){ 
    # DCS COM BETA VARIANTE NO TEMPO (normal) -----
    # modelo:
    # y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ normal
    # mu[t+1] = mu[t] + beta[t] + k1*u1[t] 
    # beta[t+1] = beta[t] + k2*u1[t] 
    # gamma[t+1] = gamma[t] + ks*u[t] 
    # > estimação via ML para densidade condicional de y t-student com variância constante no tempo
    
    otimizar <- function(y, par, Dummy){
      
      N <- length(y)
      k1 <- par[1]
      k2 <- par[2]
      ks <- par[3]
      f2 <- par[4]
      beta <- par[6]
      mu <- par[7]
      alpha <- matrix(NA, ncol = 12, nrow = N + 1)
      alpha[1,] <- c(par[8:18], - sum(par[8:18]))
      j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
      gamma <- alpha[1,j[1]]
      if(outlier){
        d <- par[19:(19+ncol(data.frame(Dummy))-1)]
      }else{
        d <- NA
      }
      u1 <- NULL
      
      if(outlier){
        for(t in 1:N){ 
          u1[t] <- (y[t] - mu[t] - gamma[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))
          mu[t+1] <- beta[t] + mu[t] + k1*u1[t]
          beta[t+1] <- beta[t] + k2*u1[t]
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
        }
        
        mu <- ts(mu, start = start(y), freq = frequency(y))
        beta <- ts(beta, start = start(y), freq = frequency(y))
        gamma <- ts(gamma, start = start(y), freq = frequency(y))
        dum <- ts(t(d %*% t(Dummy)), end = end(y), freq = frequency(y))
        
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - dum)^2)/exp(2*f2)
        -loglik
      }else{
        for(t in 1:N){ 
          u1[t] <- (y[t] - mu[t] - gamma[t])
          mu[t+1] <- beta[t] + mu[t] + k1*u1[t]
          beta[t+1] <- beta[t] + k2*u1[t]
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
        } 
        
        mu <- ts(mu, start = start(y), freq = frequency(y))
        beta <- ts(beta, start = start(y), freq = frequency(y))
        gamma <- ts(gamma, start = start(y), freq = frequency(y))
        
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma)^2)/exp(2*f2)
        -loglik
      }
      
    }
    
    if(otimo == T){
      otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y,
                           Dummy = initial$Dummy,
                           lower = initial$par$lower, upper = initial$par$upper, control = list(eval.max = 10000, iter.max = 10000))
    }else{
      otimizados <- list()
      otimizados$par <- initial$par$value
    }
    
    N <- length(y)
    k1 <- otimizados$par[1]
    k2 <- otimizados$par[2]
    ks <- otimizados$par[3]
    f2 <- otimizados$par[4]
    beta <- otimizados$par[6]
    mu <- otimizados$par[7]
    alpha <- matrix(NA, ncol = 12, nrow = N + 1)
    alpha[1,] <- c(otimizados$par[8:18], - sum(otimizados$par[8:18]))
    j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
    gamma <- alpha[1,j[1]]
    if(outlier){
      d <- otimizados$par[19:(19+ncol(data.frame(initial$Dummy))-1)]
    }else{ 
      d <- NA
    }
    u1 <- NULL
    
    if(outlier){
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))
        mu[t+1] <- beta[t] + mu[t] + k1*u1[t]
        beta[t+1] <- beta[t] + k2*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      }
      
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(beta, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      dum <- ts(t(d %*% t(initial$Dummy)), end = end(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - dum)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma - dum)/exp(f2)
      nu <- (y - mu - gamma - dum)
      score <- u1/exp(2*f2)
      out <- cbind(mu, beta, gamma, f2, exp(f2), epsilon, nu, score, u1,dum)
      colnames(out) <- c("mu","beta","gamma", "f2","sigma","epsilon", "nu", "score","u","dummy")
      
    }else{
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t])
        mu[t+1] <- beta[t] + mu[t] + k1*u1[t]
        beta[t+1] <- beta[t] + k2*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      }
      
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(beta, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma)/exp(f2)
      nu <- (y - mu - gamma)
      score <- u1/exp(2*f2)
      out <- cbind(mu, beta, gamma, f2, exp(f2), epsilon, nu, score, u1)
      colnames(out) <- c("mu","beta","gamma", "f2","sigma","epsilon", "nu", "score","u")
    }
    
    # output
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  
    
  }else if(type == "BSM2_normal"){ 
    # DCS SEM BETA (normal) -----
    # modelo:
    # y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ normal
    # mu[t+1] = mu[t] + k1*u1[t]
    # gamma[t+1] = gamma[t] + ks*u[t] 
    # > estimação via ML para densidade condicional de y t-student com variância constante no tempo
    
    otimizar <- function(y, par, Dummy){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[3]
      f2 <- par[4]
      mu <- par[7]
      alpha <- matrix(NA, ncol = 12, nrow = N + 1)
      alpha[1,] <- c(par[8:18], - sum(par[8:18]))
      j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
      gamma <- alpha[1,j[1]]
      if(outlier){
        d <- par[19:(19+ncol(data.frame(Dummy))-1)]
      }else{
        d <- NA
      }
      u1 <- NULL
      
      if(outlier){
        for(t in 1:N){ 
          u1[t] <- (y[t] - mu[t] - gamma[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))
          mu[t+1] <- mu[t] + k1*u1[t]
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
        }
        
        mu <- ts(mu, start = start(y), freq = frequency(y))
        beta <- ts(NA, start = start(y), freq = frequency(y))
        gamma <- ts(gamma, start = start(y), freq = frequency(y))
        dum <- ts(t(d %*% t(Dummy)), end = end(y), freq = frequency(y))
        
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - dum)^2)/exp(2*f2)
        -loglik
      }else{
        for(t in 1:N){ 
          u1[t] <- (y[t] - mu[t] - gamma[t])
          mu[t+1] <- mu[t] + k1*u1[t]
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
        } 
        
        mu <- ts(mu, start = start(y), freq = frequency(y))
        beta <- ts(NA, start = start(y), freq = frequency(y))
        gamma <- ts(gamma, start = start(y), freq = frequency(y))
        
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma)^2)/exp(2*f2)
        -loglik
      }
      
    }
    
    if(otimo == T){
      otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y,
                           Dummy = initial$Dummy,
                           lower = initial$par$lower, upper = initial$par$upper, control = list(eval.max = 10000, iter.max = 10000))
    }else{
      otimizados <- list()
      otimizados$par <- initial$par$value
    }
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[3]
    f2 <- otimizados$par[4]
    mu <- otimizados$par[7]
    alpha <- matrix(NA, ncol = 12, nrow = N + 1)
    alpha[1,] <- c(otimizados$par[8:18], - sum(otimizados$par[8:18]))
    j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
    gamma <- alpha[1,j[1]]
    if(outlier){
      d <- otimizados$par[19:(19+ncol(data.frame(initial$Dummy))-1)]
    }else{ 
      d <- NA
    }
    u1 <- NULL
    
    if(outlier){
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      }
      
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(NA, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      dum <- ts(t(d %*% t(initial$Dummy)), end = end(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - dum)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma - dum)/exp(f2)
      nu <- (y - mu - gamma - dum)
      score <- u1/exp(2*f2)
      out <- cbind(mu, beta, gamma, f2, exp(f2), epsilon, nu, score, u1,dum)
      colnames(out) <- c("mu","beta","gamma", "f2","sigma","epsilon", "nu", "score","u","dummy")
      
    }else{
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t])
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      }
      
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(NA, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma)/exp(f2)
      nu <- (y - mu - gamma)
      score <- u1/exp(2*f2)
      out <- cbind(mu, beta, gamma, f2, exp(f2), epsilon, nu, score, u1)
      colnames(out) <- c("mu","beta","gamma", "f2","sigma","epsilon", "nu", "score","u")
    }
    
    # output
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik)) 
  }else if(type == "BSM3_normal"){ 
    # DCS SEM BETA COM PSI (normal) -----
    # modelo:
    # y[t] = mu[t] + gamma[t] + psi[t] + exp(f2)*epsilon[t], epsilon[t] ~ normal
    # mu[t+1] = mu[t] + k1*u1[t]
    # gamma[t+1] = gamma[t] + ks*u[t] 
    # psi[t+1] = phi*psi[t] + k3*u1[t]  
    # > estimação via ML para densidade condicional de y t-student com variância constante no tempo
    
    otimizar <- function(y, par, Dummy){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[3]
      f2 <- par[4]
      mu <- par[7]
      alpha <- matrix(NA, ncol = 12, nrow = N + 1)
      alpha[1,] <- c(par[8:18], - sum(par[8:18]))
      j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
      gamma <- alpha[1,j[1]]
      psi <- par[19]
      phi <- par[20]
      k3 <- par[21]
      if(outlier){
        d <- par[22:(22+ncol(data.frame(Dummy))-1)]
      }else{
        d <- NA
      }
      u1 <- NULL
      
      if(outlier){
        for(t in 1:N){ 
          u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))
          mu[t+1] <- mu[t] + k1*u1[t]
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
          psi[t+1] <- phi*psi[t] + k3*u1[t]
        }
        
        mu <- ts(mu, start = start(y), freq = frequency(y))
        beta <- ts(NA, start = start(y), freq = frequency(y))
        gamma <- ts(gamma, start = start(y), freq = frequency(y))
        psi <- ts(psi, start = start(y), freq = frequency(y))
        dum <- ts(t(d %*% t(Dummy)), end = end(y), freq = frequency(y))
        
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - psi - dum)^2)/exp(2*f2)
        -loglik
      }else{
        for(t in 1:N){ 
          u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])
          mu[t+1] <- mu[t] + k1*u1[t]
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
          psi[t+1] <- phi*psi[t] + k3*u1[t]
        } 
        
        mu <- ts(mu, start = start(y), freq = frequency(y))
        beta <- ts(NA, start = start(y), freq = frequency(y))
        gamma <- ts(gamma, start = start(y), freq = frequency(y))
        psi <- ts(psi, start = start(y), freq = frequency(y))
        
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - psi)^2)/exp(2*f2)
        -loglik
      }
      
    }
    
    if(otimo == T){
      otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y,
                           Dummy = initial$Dummy,
                           lower = initial$par$lower, upper = initial$par$upper, control = list(eval.max = 10000, iter.max = 10000))
    }else{
      otimizados <- list()
      otimizados$par <- initial$par$value
    }
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[3]
    f2 <- otimizados$par[4]
    mu <- otimizados$par[7]
    alpha <- matrix(NA, ncol = 12, nrow = N + 1)
    alpha[1,] <- c(otimizados$par[8:18], - sum(otimizados$par[8:18]))
    j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
    gamma <- alpha[1,j[1]]
    psi <- otimizados$par[19]
    phi <- otimizados$par[20]
    k3 <- otimizados$par[21]
    if(outlier){
      d <- otimizados$par[22:(22+ncol(data.frame(initial$Dummy))-1)]
    }else{ 
      d <- NA
    }
    u1 <- NULL
    
    if(outlier){
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
        psi[t+1] <- phi*psi[t] + k3*u1[t]
      }
      
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(NA, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      psi <- ts(psi, start = start(y), freq = frequency(y))
      dum <- ts(t(d %*% t(initial$Dummy)), end = end(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - psi - dum)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma - psi - dum)/exp(f2)
      nu <- (y - mu - gamma - psi - dum)
      score <- u1/exp(2*f2)
      out <- cbind(mu, beta, gamma, psi, f2, exp(f2), epsilon, nu, score, u1,dum)
      colnames(out) <- c("mu","beta","gamma","psi", "f2","sigma","epsilon", "nu", "score","u","dummy")
      
    }else{
      for(t in 1:N){ 
        u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
        psi[t+1] <- phi*psi[t] + k3*u1[t]
      }
      
      mu <- ts(mu, start = start(y), freq = frequency(y))
      beta <- ts(NA, start = start(y), freq = frequency(y))
      gamma <- ts(gamma, start = start(y), freq = frequency(y))
      psi <- ts(psi, start = start(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - psi)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma)/exp(f2)
      nu <- (y - mu - gamma - psi)
      score <- u1/exp(2*f2)
      out <- cbind(mu, beta, gamma, psi, f2, exp(f2), epsilon, nu, score, u1)
      colnames(out) <- c("mu","beta","gamma","psi", "f2","sigma","epsilon", "nu", "score","u")
    }
    
    # output
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik)) 
    
    
  }else if(type == "BSM2_beta_norm"){ 
    # DCS COM BETA ESTATICO (NORMAL) -----
    
    # modelo:
    # y[t] = mu[t] + gamma[t] + exp(f2)*epsilon[t], epsilon[t] ~ normal(0,1)
    # mu[t+1] = beta + mu[t] + k1*u1[t]
    # gamma[t+1] = gamma[t] + k*u[t] 
    # > estimação via ML para densidade condicional de y normal com média variante e variância constante no tempo
    
    
    otimizar <- function(y, par, initial_gamma, Dummy){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[2]
      f2 <- par[3]
      mu <- par[4]
      beta <- par[5]
      alpha <- matrix(NA, ncol = 12, nrow = N)
      #alpha[1,] <- c(initial_gamma, - sum(initial_gamma))  
      alpha[1,] <- c(par[6:16], - sum(par[6:16]))  
      if(outlier){
        d <- par[17:(17+ncol(data.frame(Dummy))-1)]
      }else{
        d <- NA
      }
      j <- cycle(y)
      gamma <- alpha[1,j[1]]
      u1 <- NULL
      
      if(outlier){
        for(t in 1:(N-1)){ 
          u1[t] <- (y[t] - mu[t] - gamma[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))
          mu[t+1] <- beta + mu[t] + k1*u1[t] 
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
        } 
        t <- N
        u1[t] <- (y[t] - mu[t] - gamma[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))
        
        mu <- ts(mu, end = end(y), freq = frequency(y))
        gamma <- ts(gamma, end = end(y), freq = frequency(y))
        dum <- ts(t(d %*% t(Dummy)), end = end(y), freq = frequency(y))
        
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - dum)^2)/exp(2*f2)
        -loglik
      }else{
        for(t in 1:(N-1)){ 
          u1[t] <- (y[t] - mu[t] - gamma[t])
          mu[t+1] <- beta + mu[t] + k1*u1[t] 
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
        } 
        t <- N
        u1[t] <- (y[t] - mu[t] - gamma[t])
        
        mu <- ts(mu, end = end(y), freq = frequency(y))
        gamma <- ts(gamma, end = end(y), freq = frequency(y))
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma)^2)/exp(2*f2)
        -loglik
      }
      
    }
    
    if(otimo == T){
      otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y, 
                           initial_gamma = initial$gamma, Dummy = initial$Dummy,
                           lower = initial$par$lower, upper = initial$par$upper, 
                           control = list(eval.max = 10000, iter.max = 10000))
      
    }else{
      otimizados <- list()
      otimizados$par = initial$par$value
    }
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[2]
    f2 <- otimizados$par[3]
    mu <- otimizados$par[4]
    beta <- otimizados$par[5]
    alpha <- matrix(NA, ncol = 12, nrow = N)
    #alpha[1,] <- c(initial$gamma, - sum(initial$gamma))  
    alpha[1,] <- c(otimizados$par[6:16], - sum(otimizados$par[6:16]))  
    if(outlier){
      d <- otimizados$par[17:(17+ncol(data.frame(initial$Dummy))-1)]
    }else{ 
      d <- NA
    }
    j <- cycle(y)
    gamma <- alpha[1,j[1]]
    u1 <- NULL
    
    if(outlier){
      for(t in 1:(N-1)){
        u1[t] <- (y[t] - mu[t] - gamma[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))
        mu[t+1] <- beta + mu[t] + k1*u1[t] 
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      } 
      t <- N
      u1[t] <- (y[t] - mu[t] - gamma[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))
      mu <- ts(mu, end = end(y), freq = frequency(y))
      gamma <- ts(gamma, end = end(y), freq = frequency(y))
      dum <- ts(t(d %*% t(initial$Dummy)), end = end(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - dum)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma - dum)/exp(f2)
      nu <- (y - mu - gamma - dum)
      score <- u1/exp(2*f2)
      out <- cbind(mu, gamma, f2, exp(f2), epsilon, nu, score, u1,dum)
      colnames(out) <- c("mu","gamma", "f2","sigma","epsilon", "nu", "score","u","dummy")
    }else{
      for(t in 1:(N-1)){ 
        u1[t] <- (y[t] - mu[t] - gamma[t])
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
      } 
      t <- N
      u1[t] <- (y[t] - mu[t] - gamma[t])
      
      mu <- ts(mu, end = end(y), freq = frequency(y))
      gamma <- ts(gamma, end = end(y), freq = frequency(y))
      u1 <- ts(u1, end = end(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma)/exp(f2)
      nu <- (y - mu - gamma)
      score <- u1/exp(2*f2)
      out <- cbind(mu, gamma, f2, exp(f2), epsilon, nu, score, u1)
      colnames(out) <- c("mu","gamma", "f2","sigma","epsilon", "nu", "score","u")
    }
    
    # output
    #print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
    
  }else if(type == "BSM2_beta_psi_norm"){ 
    # DCS COM BETA ESTATICO E PSI (NORMAL) -----
    
    # modelo:
    # y[t] = mu[t] + gamma[t] + psi[t] + exp(f2)*epsilon[t], epsilon[t] ~ normal(0,1)
    # mu[t+1] = beta + mu[t] + k1*u1[t]
    # gamma[t+1] = gamma[t] + k*u[t] 
    # psi[t+1] = phi*psi[t] + k3*u[t]
    # > estimação via ML para densidade condicional de y t-student com média constante e variância constante no tempo
    
    
    otimizar <- function(y, par, initial_gamma, Dummy){
      
      N <- length(y)
      k1 <- par[1]
      ks <- par[2]
      f2 <- par[3]
      mu <- par[4]
      beta <- par[5]
      psi <- par[6]
      phi <- par[7]
      k3 <- par[8]
      alpha <- matrix(NA, ncol = 12, nrow = N)
      #alpha[1,] <- c(initial_gamma, - sum(initial_gamma))  
      alpha[1,] <- c(par[9:19], - sum(par[9:19]))  
      if(outlier){
        d <- par[20:(20+ncol(data.frame(Dummy))-1)]
      }else{
        d <- NA
      }
      j <- cycle(y)
      gamma <- alpha[1,j[1]]
      u1 <- NULL
      
      if(outlier){
        for(t in 1:(N-1)){ 
          u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))
          mu[t+1] <- beta + mu[t] + k1*u1[t] 
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
          psi[t+1] <- phi*psi[t] + k3*u1[t]
        } 
        t <- N
        u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))
        
        mu <- ts(mu, end = end(y), freq = frequency(y))
        gamma <- ts(gamma, end = end(y), freq = frequency(y))
        psi <- ts(psi, end = end(y), freq = frequency(y))
        dum <- ts(t(d %*% t(Dummy)), end = end(y), freq = frequency(y))
        
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - psi - dum)^2)/exp(2*f2)
        -loglik
      }else{
        for(t in 1:(N-1)){ 
          u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])
          mu[t+1] <- beta + mu[t] + k1*u1[t] 
          k <- rep(-ks/11,12) 
          k[j[t]] <- ks
          alpha[t+1,] <- alpha[t,] + k*u1[t]
          gamma[t+1] <- alpha[t+1,j[t+1]]
          psi[t+1] <- phi*psi[t] + k3*u1[t]
        } 
        t <- N
        u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])
        
        mu <- ts(mu, end = end(y), freq = frequency(y))
        gamma <- ts(gamma, end = end(y), freq = frequency(y))
        psi <- ts(psi, end = end(y), freq = frequency(y))
        # loglik
        loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - psi)^2)/exp(2*f2)
        -loglik
      }
      
    }
    
    if(otimo == T){
      otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y, 
                           initial_gamma = initial$gamma, Dummy = initial$Dummy,
                           lower = initial$par$lower, upper = initial$par$upper, 
                           control = list(eval.max = 10000, iter.max = 10000))
    }else{
      otimizados <- list()
      otimizados$par = initial$par$value
    }
    
    N <- length(y)
    k1 <- otimizados$par[1]
    ks <- otimizados$par[2]
    f2 <- otimizados$par[3]
    mu <- otimizados$par[4]
    beta <- otimizados$par[5]
    psi <- otimizados$par[6]
    phi <- otimizados$par[7]
    k3 <- otimizados$par[8]
    alpha <- matrix(NA, ncol = 12, nrow = N)
    #alpha[1,] <- c(initial$gamma, - sum(initial$gamma))  
    alpha[1,] <- c(otimizados$par[9:19], - sum(otimizados$par[9:19]))  
    if(outlier){
      d <- otimizados$par[20:(20+ncol(data.frame(initial$Dummy))-1)]
    }else{ 
      d <- NA
    }
    j <- cycle(y)
    gamma <- alpha[1,j[1]]
    u1 <- NULL
    
    if(outlier){
      for(t in 1:(N-1)){
        u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))
        mu[t+1] <- beta + mu[t] + k1*u1[t] 
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
        psi[t+1] <- phi*psi[t] + k3*u1[t]
      } 
      t <- N
      u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))
      mu <- ts(mu, end = end(y), freq = frequency(y))
      gamma <- ts(gamma, end = end(y), freq = frequency(y))
      psi <- ts(psi, end = end(y), freq = frequency(y))
      dum <- ts(t(d %*% t(initial$Dummy)), end = end(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - psi - dum)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma - psi - dum)/exp(f2)
      nu <- (y - mu - gamma - psi - dum)
      score <- u1/exp(2*f2)
      out <- cbind(mu, gamma, psi, psi + mu, f2, exp(f2), epsilon, nu, score, u1,dum)
      colnames(out) <- c("mu","gamma","psi", "psi_mu", "f2","sigma","epsilon", "nu", "score","u","dummy")
    }else{
      for(t in 1:(N-1)){ 
        u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])
        mu[t+1] <- mu[t] + k1*u1[t]
        k <- rep(-ks/11,12) 
        k[j[t]] <- ks
        alpha[t+1,] <- alpha[t,] + k*u1[t]
        gamma[t+1] <- alpha[t+1,j[t+1]]
        psi[t+1] <- phi*psi[t] + k3*u1[t]
      } 
      t <- N
      u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])
      
      mu <- ts(mu, end = end(y), freq = frequency(y))
      gamma <- ts(gamma, end = end(y), freq = frequency(y))
      psi <- ts(psi, end = end(y), freq = frequency(y))
      u1 <- ts(u1, end = end(y), freq = frequency(y))
      loglik <- -(N/2)*log(2*pi*exp(2*f2)) - 0.5*sum((y - mu - gamma - psi)^2)/exp(2*f2)
      epsilon <- (y - mu - gamma - psi)/exp(f2)
      nu <- (y - mu - gamma - psi)
      score <- u1/exp(2*f2)
      out <- cbind(mu, gamma, psi, psi + mu, f2, exp(f2), epsilon, nu, score, u1)
      colnames(out) <- c("mu","gamma","psi", "psi_mu", "f2","sigma","epsilon", "nu", "score","u")
    }
    
    # output
    #print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
    
  }
  
}

# }else if(type == "BSM2_beta_psi"){ 
#   # DCS COM BETA ESTATICO E PSI -----
#   
#   # modelo:
#   # y[t] = mu[t] + gamma[t] + psi[t] + exp(f2)*epsilon[t], epsilon[t] ~ t(v)
#   # mu[t+1] = beta + mu[t] + k1*u1[t]
#   # gamma[t+1] = gamma[t] + k*u[t] 
#   # psi[t+1] = phi*psi[t] + k3*u[t]
#   # > estimação via ML para densidade condicional de y t-student com média constante e variância constante no tempo
#   
#   
#   otimizar <- function(y, par, initial_gamma, Dummy){
#     
#     N <- length(y)
#     k1 <- par[1]
#     ks <- par[2]
#     f2 <- par[3]
#     df <- par[4]
#     mu <- par[5]
#     beta <- par[6]
#     psi <- par[7]
#     phi <- par[8]
#     k3 <- par[9]
#     alpha <- matrix(NA, ncol = 12, nrow = N+1)
#     #alpha[1,] <- c(initial_gamma, - sum(initial_gamma))  
#     alpha[1,] <- c(par[10:20], - sum(par[10:20]))  
#     if(outlier){
#       d <- par[21:(21+ncol(data.frame(Dummy))-1)]
#     }else{
#       d <- NA
#     }
#     j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
#     gamma <- alpha[1,j[1]]
#     u1 <- NULL
#     
#     if(outlier){
#       for(t in 1:N){ 
#         u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))/(1 + (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))^2/(df*exp(2*f2)))
#         mu[t+1] <- beta + mu[t] + k1*u1[t] 
#         k <- rep(-ks/11,12) 
#         k[j[t]] <- ks
#         alpha[t+1,] <- alpha[t,] + k*u1[t]
#         gamma[t+1] <- alpha[t+1,j[t+1]]
#         psi[t+1] <- phi*psi[t] + k3*u1[t]
#       } 
#       # t <- N
#       # u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))/(1 + (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))^2/(df*exp(2*f2)))
#       
#       mu <- ts(mu, start = start(y), freq = frequency(y))
#       gamma <- ts(gamma, start = start(y), freq = frequency(y))
#       psi <- ts(psi, start = start(y), freq = frequency(y))
#       dum <- ts(t(d %*% t(Dummy)), start = start(y), freq = frequency(y))
#       
#       # loglik
#       loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma - psi - dum)^2/(df*exp(2*f2))))
#       -loglik
#     }else{
#       for(t in 1:(N)){ 
#         u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])/(1 + (y[t] - mu[t] - gamma[t] - psi[t])^2/(df*exp(2*f2)))
#         mu[t+1] <- beta + mu[t] + k1*u1[t] 
#         k <- rep(-ks/11,12) 
#         k[j[t]] <- ks
#         alpha[t+1,] <- alpha[t,] + k*u1[t]
#         gamma[t+1] <- alpha[t+1,j[t+1]]
#         psi[t+1] <- phi*psi[t] + k3*u1[t]
#       } 
#       t <- N
#       u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])/(1 + (y[t] - mu[t] - gamma[t] - psi[t])^2/(df*exp(2*f2)))
#       
#       mu <- ts(mu, start = start(y), freq = frequency(y))
#       gamma <- ts(gamma, start = start(y), freq = frequency(y))
#       psi <- ts(psi, start = start(y), freq = frequency(y))
#       # loglik
#       loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma - psi)^2/(df*exp(2*f2))))
#       -loglik
#     }
#     
#   }
#   
#   
#   if(otimo == T){
#     otimizados <- nlminb(start = initial$par$value, objective = otimizar, y = y, 
#                          initial_gamma = initial$gamma, Dummy = initial$Dummy,
#                          lower = initial$par$lower, upper = initial$par$upper, 
#                          control = list(eval.max = 10000, iter.max = 10000))
#   }else{
#     otimizados <- list()
#     otimizados$par = initial$par$value
#   }
#   
#   N <- length(y)
#   k1 <- otimizados$par[1]
#   ks <- otimizados$par[2]
#   f2 <- otimizados$par[3]
#   df <- otimizados$par[4]
#   mu <- otimizados$par[5]
#   beta <- otimizados$par[6]
#   psi <- otimizados$par[7]
#   phi <- otimizados$par[8]
#   k3 <- otimizados$par[9]
#   alpha <- matrix(NA, ncol = 12, nrow = N+1)
#   #alpha[1,] <- c(initial$gamma, - sum(initial$gamma))  
#   alpha[1,] <- c(otimizados$par[10:20], - sum(otimizados$par[10:20]))  
#   if(outlier){
#     d <- otimizados$par[21:(21+ncol(data.frame(initial$Dummy))-1)]
#   }else{ 
#     d <- NA
#   }
#   j <- cycle(ts(c(y,NA), start = start(y), freq = 12))
#   gamma <- alpha[1,j[1]]
#   u1 <- NULL
#   
#   if(outlier){
#     for(t in 1:(N)){ 
#       u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))/(1 + (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))^2/(df*exp(2*f2)))
#       mu[t+1] <- beta + mu[t] + k1*u1[t] 
#       k <- rep(-ks/11,12) 
#       k[j[t]] <- ks
#       alpha[t+1,] <- alpha[t,] + k*u1[t]
#       gamma[t+1] <- alpha[t+1,j[t+1]]
#       psi[t+1] <- phi*psi[t] + k3*u1[t]
#     } 
#     # t <- N
#     # u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))/(1 + (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(initial$Dummy)[t,]))^2/(df*exp(2*f2)))
#     
#     mu <- ts(mu, start = start(y), freq = frequency(y))
#     gamma <- ts(gamma, start = start(y), freq = frequency(y))
#     psi <- ts(psi, start = start(y), freq = frequency(y))
#     dum <- ts(t(d %*% t(initial$Dummy)), start = start(y), freq = frequency(y))
#     loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma - psi - dum)^2/(df*exp(2*f2))))
#     epsilon <- (y - mu - gamma - psi - dum)/exp(f2)
#     nu <- (y - mu - gamma - psi - dum)
#     b <- ((y - mu - gamma - psi - dum)^2/(df*exp(2*f2)))/(1 + (y - mu - gamma - psi - dum)^2/(df*exp(2*f2)))
#     score <- ((df + 1)/(df*exp(2*f2)))*u1
#     out <- cbind(mu, gamma, psi, psi + mu, f2, exp(f2), epsilon, nu, score, u1, b,dum)
#     colnames(out) <- c("mu","gamma","psi", "psi_mu", "f2","sigma","epsilon", "nu", "score","u", "b","dummy")
#   }else{
#     for(t in 1:(N)){ 
#       u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])/(1 + (y[t] - mu[t] - gamma[t] - psi[t])^2/(df*exp(2*f2)))
#       mu[t+1] <- mu[t] + k1*u1[t]
#       k <- rep(-ks/11,12) 
#       k[j[t]] <- ks
#       alpha[t+1,] <- alpha[t,] + k*u1[t]
#       gamma[t+1] <- alpha[t+1,j[t+1]]
#       psi[t+1] <- phi*psi[t] + k3*u1[t]
#     } 
#     # t <- N
#     # u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t])/(1 + (y[t] - mu[t] - gamma[t] - psi[t])^2/(df*exp(2*f2)))
#     
#     mu <- ts(mu, start = start(y), freq = frequency(y))
#     gamma <- ts(gamma, start = start(y), freq = frequency(y))
#     psi <- ts(psi, start = start(y), freq = frequency(y))
#     u1 <- ts(u1, start = start(y), freq = frequency(y))
#     loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma - psi)^2/(df*exp(2*f2))))
#     epsilon <- (y - mu - gamma - psi)/exp(f2)
#     nu <- (y - mu - gamma - psi)
#     b <- ((y - mu - gamma - psi)^2/(df*exp(2*f2)))/(1 + (y - mu - gamma - psi)^2/(df*exp(2*f2)))
#     score <- ((df + 1)/(df*exp(2*f2)))*u1
#     out <- cbind(mu, gamma, psi, psi + mu, f2, exp(f2), epsilon, nu, score, u1, b)
#     colnames(out) <- c("mu","gamma","psi", "psi_mu", "f2","sigma","epsilon", "nu", "score","u", "b")
#   }
#   
#   # output
#   #print(otimizados)
#   invisible(list(out = out, otimizados = otimizados, loglik = -loglik))