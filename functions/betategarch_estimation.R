betategarch_estimation <- function(y, initial = NULL, type = "var"){
  
  
  if(type == "var"){ # modelo y[t] = exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
    # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
    otimizar <- function(y, par){
      w1 <- par[1]
      A1 <- par[2]
      B1 <- par[3]
      df <- par[4]
      f2 <- w1/(1-B1)
      N <- length(y)
      u <- NULL
      for(t in 1:N){
        u[t] <- (((df + 1)*y[t]^2) / (df*exp(2*f2[t]) + y[t]^2)) - 1
        f2[t+1] <- w1 + A1*u[t] + B1*f2[t]
      }
      
      f2 <- ts(f2[-(N+1)], start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + y^2/(df*exp(2*f2))))
      -loglik
    }
    
    otimizados <- nlminb(start = c(initial), objective = otimizar, y = y, 
                         lower = c(-Inf,-Inf,-1,2), upper = c(Inf,Inf,1,Inf), control = list(eval.max = 10000, iter.max = 10000))
    
    w1 <- otimizados$par[1]
    A1 <- otimizados$par[2]
    B1 <- otimizados$par[3]
    df <- otimizados$par[4]
    
    f2 <- w1/(1-B1)
    N <- length(y)
    u <- NULL
    for(t in 1:N){
      u[t] <- (((df + 1)*y[t]^2) / (df*exp(2*f2[t]) + y[t]^2)) - 1
      f2[t+1] <- w1 + A1*u[t] + B1*f2[t]
    }
    f2 <- ts(f2[-(N+1)], start = start(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + y^2/(df*exp(2*f2))))
    epsilon <- y / exp(f2)
    out <- cbind(f2, exp(f2), epsilon)
    colnames(out) <- c("f2", "sigma","epsilon")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }else if(type == "mean-var"){ # modelo y[t] = f1 + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
    # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
    otimizar <- function(y, par){
      w1 <- par[1]
      A1 <- par[2]
      B1 <- par[3]
      df <- par[4]
      f1 <- par[5]
      f2 <- w1/(1-B1)
      N <- length(y)
      u <- NULL
      for(t in 1:N){
        u[t] <- (((df + 1)*(y[t] - f1)^2) / (df*exp(2*f2[t]) + (y[t] - f1)^2)) - 1
        f2[t+1] <- w1 + A1*u[t] + B1*f2[t]
      }
      
      f2 <- ts(f2[-(N+1)], start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
      -loglik
    }
    
    otimizados <- nlminb(start = c(initial), objective = otimizar, y = y, 
                         lower = c(-Inf,-Inf,-1,2,-Inf), upper = c(Inf,Inf,1,Inf,Inf), control = list(eval.max = 10000, iter.max = 10000))
    
    w1 <- otimizados$par[1]
    A1 <- otimizados$par[2]
    B1 <- otimizados$par[3]
    df <- otimizados$par[4]
    f1 <- otimizados$par[5]
    f2 <- w1/(1-B1)
    N <- length(y)
    u <- NULL
    for(t in 1:N){
      u[t] <- (((df + 1)*(y[t] - f1)^2) / (df*exp(2*f2[t]) + (y[t] - f1)^2)) - 1
      f2[t+1] <- w1 + A1*u[t] + B1*f2[t]
    }
    
    f2 <- ts(f2[-(N+1)], start = start(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
    epsilon <- (y - f1)/exp(f2)
    out <- cbind(f1, f2, exp(f2), epsilon)
    colnames(out) <- c("f1", "f2", "sigma","epsilon")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))

  }else if(type == "mean-var2"){ # modelo y[t] = f1[t] + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
    # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
    otimizar <- function(y, par){
      w1 <- par[1]
      w2 <- par[2]
      A1 <- par[3]
      A2 <- par[4]
      B1 <- par[5]
      B2 <- par[6]
      df <- par[7]
      f1 <- w1/(1-B1)
      f2 <- w2/(1-B2)
      N <- length(y)
      u1 <- NULL
      u2 <- NULL
      for(t in 1:N){
        u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2[t])))
        u2[t] <- (((df + 1)*(y[t] - f1[t])^2) / (df*exp(2*f2[t]) + (y[t] - f1[t])^2)) - 1
        f1[t+1] <- w1 + A1*u1[t] + B1*f1[t]
        f2[t+1] <- w2 + A2*u2[t] + B2*f2[t]
      }
      
      f1 <- ts(f1[-(N+1)], start = start(y), freq = frequency(y))
      f2 <- ts(f2[-(N+1)], start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
      -loglik
    }
    
    
    otimizados <- nlminb(start = c(initial), objective = otimizar, y = y, 
                         lower = c(-Inf,-Inf,-Inf,-Inf,-1,-1,2), upper = c(Inf,Inf,Inf,Inf,1,1,Inf), control = list(eval.max = 10000, iter.max = 10000))
    
    w1 <- otimizados$par[1]
    w2 <- otimizados$par[2]
    A1 <- otimizados$par[3]
    A2 <- otimizados$par[4]
    B1 <- otimizados$par[5]
    B2 <- otimizados$par[6]
    df <- otimizados$par[7]
    f1 <- w1/(1-B1)
    f2 <- w2/(1-B2)
    N <- length(y)
    u1 <- NULL
    u2 <- NULL
    for(t in 1:N){
      u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2[t])))
      u2[t] <- (((df + 1)*(y[t] - f1[t])^2) / (df*exp(2*f2[t]) + (y[t] - f1[t])^2)) - 1
      f1[t+1] <- w1 + A1*u1[t] + B1*f1[t]
      f2[t+1] <- w2 + A2*u2[t] + B2*f2[t]
    }
    f1 <- ts(f1[-(N+1)], start = start(y), freq = frequency(y))
    f2 <- ts(f2[-(N+1)], start = start(y), freq = frequency(y))
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
    epsilon <- (y - f1)/exp(f2)
    out <- cbind(f1, f2, exp(f2), epsilon)
    colnames(out) <- c("f1", "f2", "sigma","epsilon")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
    
  }else if(type == "mean-var3"){ # modelo y[t] = f1[t] + exp(f2[t])*epsilon[t], epsilon[t] ~ t(v)
    # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
    otimizar <- function(y, par){
      w1 <- par[1]
      w2 <- par[2]
      A1_00 <- par[3]
      A1_01 <- par[4]
      A1_11 <- 0
      A1_12 <- 0
      B1_00 <- par[5]
      B1_01 <- par[6]
      B1_11 <- 0
      B1_12 <- 0
      A2_00 <- par[7]
      B2_00 <- par[8]
      df <- par[9]
      f1 <- rep(w1/(1 - B1_00),13)
      f2 <- rep(w2/(1 - B2_00),13)
      N <- length(y)
      u1 <- NULL
      u2 <- NULL
      for(t in 1:12){
      u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2[t])))
      u2[t] <- (((df + 1)*(y[t] - f1[t])^2) / (df*exp(2*f2[t]) + (y[t] - f1[t])^2)) - 1
      }
      for(t in 13:N){
        u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2[t])))
        u2[t] <- (((df + 1)*(y[t] - f1[t])^2) / (df*exp(2*f2[t]) + (y[t] - f1[t])^2)) - 1
        f1[t+1] <- w1 + A1_00*u1[t] + A1_01*u1[t-1] + A1_11*u1[t-11] + A1_12*u1[t-12] + B1_00*f1[t] + B1_01*f1[t-1] + B1_11*f1[t-11] + B1_12*f1[t-12]
        f2[t+1] <- w2 + A2_00*u2[t] + B2_00*f2[t]
      }
      
      f1 <- ts(f1[13:N], start = start(y), freq = frequency(y))
      f2 <- ts(f2[13:N], start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
      -loglik
    }
    
    
    otimizados <- nlminb(start = c(initial), objective = otimizar, y = y, 
                         lower = c(-Inf,-Inf,-Inf,-Inf,-1,-1,-Inf,-1,2), upper = c(Inf,Inf,Inf,Inf,1,1,Inf,1,Inf), control = list(eval.max = 10000, iter.max = 10000))
    
    w1 <- otimizados$par[1]
    w2 <- otimizados$par[2]
    A1_00 <- otimizados$par[3]
    A1_01 <- otimizados$par[4]
    A1_11 <- 0
    A1_12 <- 0
    B1_00 <- otimizados$par[5]
    B1_01 <- otimizados$par[6]
    B1_11 <- 0
    B1_12 <- 0
    A2_00 <- otimizados$par[7]
    B2_00 <- otimizados$par[8]
    df <- otimizados$par[9]
    f1 <- rep(w1/(1 - B1_00),13)
    f2 <- rep(w2/(1 - B2_00),13)
    N <- length(y)
    u1 <- NULL
    u2 <- NULL
    for(t in 1:12){
      u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2[t])))
      u2[t] <- (((df + 1)*(y[t] - f1[t])^2) / (df*exp(2*f2[t]) + (y[t] - f1[t])^2)) - 1
    }
    for(t in 13:N){
      u1[t] <- (df + 1)*((y[t] - f1[t])/(df*exp(2*f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2[t])))
      u2[t] <- (((df + 1)*(y[t] - f1[t])^2) / (df*exp(2*f2[t]) + (y[t] - f1[t])^2)) - 1
      f1[t+1] <- w1 + A1_00*u1[t] + A1_01*u1[t-1] + A1_11*u1[t-11] + A1_12*u1[t-12] + B1_00*f1[t] + B1_01*f1[t-1] + B1_11*f1[t-11] + B1_12*f1[t-12]
      f2[t+1] <- w2 + A2_00*u2[t] + B2_00*f2[t]
    }
    f1 <- ts(f1[13:N], end = end(y), freq = frequency(y))
    f2 <- ts(f2[13:N], end = end(y), freq = frequency(y))
    
    loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - sum(f2) - ((df + 1)/2)*sum(log(1 + (y - f1)^2/(df*exp(2*f2))))
    epsilon <- (y - f1)/exp(f2)
    out <- cbind(f1, f2, exp(f2), epsilon)
    colnames(out) <- c("f1", "f2", "sigma","epsilon")
    
    # output
    print(otimizados)
    invisible(list(out = out, otimizados = otimizados, loglik = -loglik))
  }
}
