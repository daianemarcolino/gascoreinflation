simularGAS <- function(density = "normal", n = 100, link = FALSE, mean = "variavel", sd = "fixa", seed = 123, mu = NULL, sigma2 = NULL, w = NULL, A1 = NULL, B1 = NULL, df = NULL){
  
  if(density == "normal" & link == FALSE){
    
    if(mean == "variavel" & sd == "fixa"){
      
      # NORMAL: média variável, variância fixa, link FALSE ---------------------------
      
      # > estimação de y com densidade condicional de normal com média variante no tempo e variância fixa
      # parâmetros estáticos: c(sigma2 = var y, w, A1, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n+1, mean = 0, sd = sqrt(sigma2)))
      
      # valor inicial de f1 = 0
      f1 <- 0
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- f1[t] + epsilon[t]
        f1[t+1] <- w + A1[1]*(y[t] - f1[t]) + B1[1]*f1[t]
      }
      
      y <- ts(y[2:(n+1)])
      f1 <- ts(f1[2:(n+1)])
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, mu e epsilon
      return(cbind(y, f1, epsilon))
      
    }else if(mean == "fixa" & sd == "variavel"){
      
      # NORMAL: média fixa, variância variável, link FALSE ---------------------------
      
      # > estimação de y com densidade condicional de normal com média fixa e variância variável no tempo
      # parâmetros estáticos: c(mu = mean y, w, A1, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n+1, mean = 0, sd = 1))
      
      # valor inicial de sigma = 0.1
      f2 <- 0.1
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- mu + sqrt(f2[t])*epsilon[t]
        f2[t+1] <- w + A1[1]*((y[t] -  mu)^2 - f2[t]) + B1[1]*f2[t]
      }
      
      y <- ts(y[2:(n+1)])
      f2 <- ts(f2[2:(n+1)])
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, f2, epsilon))
      
    }else if(mean == "variavel" & sd == "variavel"){
      
      # NORMAL: média fixa, variância variável, link FALSE ---------------------------
      
      # > estimação de y com densidade condicional de normal com média e variância variantes no tempo
      # parâmetros estáticos: c(w1, A1, B1, w2, A2, B2)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n+1, mean = 0, sd = 1))
      
      # valor inicial de mu = 0
      # valor inicial de sigma = 0
      f1 <- 0
      f2 <- 0.1
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- f1[t-1] + sqrt(f2[t-1])*epsilon[t]
        f1[t+1] <- w[1] + A1[1]*(y[t] -  f1[t]) + B1[1]*f1[t]
        f2[t+1] <- w[2] + A1[2]*((y[t] -  f1[t])^2 - f2[t]) + B1[2]*f2[t]
      }
      
      y <- ts(y[2:(n+1)])
      f1 <- ts(f1[2:(n+1)])
      f2 <- ts(f2[2:(n+1)])
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, f1, f2, epsilon))
    }
    
  }else if(density == "normal" & link == TRUE){
    
    if(mean == "variavel" & sd == "fixa"){
      
      # NORMAL: média variável, variância fixa, link TRUE ---------------------------
      
      # > estimação de y com densidade condicional de normal com média variante no tempo e variância fixa
      # parâmetros estáticos: c(sigma2 = var y, w, A1, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n+1, mean = 0, sd = sqrt(sigma2)))
      
      # valor inicial de f1 = 0
      f1 <- 0
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- f1[t] + epsilon[t]
        f1[t+1] <- w + A1[1]*(y[t] - f1[t]) + B1[1]*f1[t]
      }
      
      y <- ts(y[2:(n+1)])
      f1 <- ts(f1[2:(n+1)])
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, mu e epsilon
      return(cbind(y, f1, epsilon))
      
    }else if(mean == "fixa" & sd == "variavel"){
      
      # NORMAL: média fixa, variância variável, link TRUE ---------------------------
      
      # > estimação de y com densidade condicional de normal com média fixa e variância variável no tempo
      # parâmetros estáticos: c(mu = mean y, w, A1, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n+1, mean = 0, sd = 1))
      
      # valor inicial de sigma = 0.1
      f2 <- 0.1
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- mu + exp(f2[t]/2)*epsilon[t]
        f2[t+1] <- w + A1[1]*((y[t] -  mu)^2/exp(f2[t]) - 1) + B1[1]*f2[t]
      }
      
      y <- ts(y[2:(n+1)])
      f2 <- exp(ts(f2[2:(n+1)]))
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, f2, epsilon))
      
    }else if(mean == "variavel" & sd == "variavel"){
      
      # NORMAL: média fixa, variância variável, link TRUE ---------------------------
      
      # > estimação de y com densidade condicional de normal com média e variância variantes no tempo
      # parâmetros estáticos: c(w1, A1, B1, w2, A2, B2)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n+1, mean = 0, sd = 1))
      
      # valor inicial de mu = 0
      # valor inicial de sigma = 0.1
      f1 <- 0
      f2 <- 0.1
      y <- NA
      for(t in 1:(n+1)){
        y[t] <- mu + exp(f2[t]/2)*epsilon[t]
        f1[t+1] <- w + A1[1]*(y[t] - f1[t]) + B1[1]*f1[t]
        f2[t+1] <- w + A1[2]*((y[t] -  mu)^2/exp(f2[t]) - 1) + B1[2]*f2[t]
      }
      
      y <- ts(y[2:(n+1)])
      f1 <- ts(f1[2:(n+1)])
      f2 <- exp(ts(f2[2:(n+1)]))
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, f1, f2, epsilon))
    }
    
  }else if(density == "t" & link == FALSE){
    
    if(mean == "variavel" & sd == "fixa"){
      
      # T: média variável, variância fixa, link FALSE ---------------------------
      
      # > estimação de y com densidade condicional t com média variante no tempo e variância fixa
      # parâmetros estáticos: c(sigma2 = var y, w, A1, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rt(n+1, df = df))
      
      # valor inicial de f1 = 0
      f1 <- 0
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- f1[t] + epsilon[t]
        f1[t+1] <- w + A1[1]*( 1/(df + 1) * (y[t] - f1[t])/(1 + (y[t] - f1[t])^2/df) * beta(1/2, df/2)/beta(3/2, (df + 2)/2)) + B1[1]*f1[t]
      }
      
      y <- ts(y[2:(n+1)])
      f1 <- ts(f1[2:(n+1)])
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, mu e epsilon
      return(cbind(y, f1, epsilon))
      
    }else if(mean == "fixa" & sd == "variavel"){
      
      # T: média fixa, variância fixa, link FALSE ---------------------------
      
      # > estimação de y com densidade condicional t com média fixa e variância variante no tempo
      # theta (parâmetros estáticos): c(mean y, w, A0, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rt(n+1, df = df))
      
      # valor inicial de sigma = 0.1
      f2 <- 0.1
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- mu + sqrt(f2[t])*epsilon[t]
        f2[t+1] <- w + A1[1]*((df + 3)/df * f2[t] * (((df + 1)* ((y[t] - mu)^2/(df*f2[t]))/(1 + (y[t] - mu)^2/(df*f2[t]))) - 1)) + B1[1]*f2[t]
      }
      
      y <- ts(y[2:(n+1)])
      f2 <- ts(f2[2:(n+1)])
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, f2, epsilon))
      
    }else if(mean == "variavel" & sd == "variavel"){
      
      #  T: média variável, variância variável, link FALSE ---------------------------
      
      # > estimação de y com densidade condicional t com média e variância variante no tempo
      # theta (parâmetros estáticos): c(w, A0, B1, w, A0, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rt(n+1, df = df))
      
      # valor inicial de mu = 0
      # valor inicial de sigma = 0.1
      f1 <- 0
      f2 <- 0.1
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- f1[t] + sqrt(f2[t])*epsilon[t]
        f1[t+1] <- w + A1[1]*( 1/(df + 1) * (y[t] - f1[t])/(1 + (y[t] - f1[t])^2/df) * beta(1/2, df/2)/beta(3/2, (df + 2)/2)) + B1[1]*f1[t]
        f2[t+1] <- w + A1[2]*((df + 3)/df * f2[t] * (((df + 1)* ((y[t] - f1[t])^2/(df*f2[t]))/(1 + (y[t] - f1[t])^2/(df*f2[t]))) - 1)) + B1[2]*f2[t]
      }
      
      y <- ts(y[2:(n+1)])
      f1 <- ts(f1[2:(n+1)])
      f2 <- ts(f2[2:(n+1)])
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, f1, f2, epsilon))
    }
    
  }else if(density == "t" & link == TRUE){
    
    if(mean == "variavel" & sd == "fixa"){
      
      #  T: média variável, variância fixa, link TRUE ---------------------------
      
      # > estimação de y com densidade condicional t com média variante no tempo e variância fixa
      # parâmetros estáticos: c(sigma2 = var y, w, A1, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rt(n+1, df = df))
      
      # valor inicial de f1 = 0
      f1 <- 0
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- f1[t] + epsilon[t]
        f1[t+1] <- w + A1[1]*( 1/(df + 1) * (y[t] - f1[t])/(1 + (y[t] - f1[t])^2/df) * beta(1/2, df/2)/beta(3/2, (df + 2)/2)) + B1[1]*f1[t]
      }
      
      y <- ts(y[2:(n+1)])
      f1 <- ts(f1[2:(n+1)])
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, mu e epsilon
      return(cbind(y, f1, epsilon))
      
      
    }else if(mean == "fixa" & sd == "variavel"){
      
      #  T: média fixa, variância variável, link TRUE ---------------------------
      
      # > estimação de y com densidade condicional t com média fixa e variância variante no tempo
      # theta (parâmetros estáticos): c(mean y, w, A0, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rt(n+1, df = df))
      
      # valor inicial de sigma = 0.1
      f2 <- 0.1
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- mu + sqrt(exp(f2[t]))*epsilon[t]
        f2[t+1] <- w + A1[1]*((df + 3)/df * (((df + 1)*((y[t] - mu)^2/(df*exp(f2[t])))/(1 + (y[t] - mu)^2/(df*exp(f2[t])))) - 1)) + B1[1]*f2[t]
      }
      
      y <- ts(y[2:(n+1)])
      f2 <- exp(ts(f2[2:(n+1)]))
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, f2, epsilon))
      
    }else if(mean == "variavel" & sd == "variavel"){
      
      #  T: média variável, variância variável, link TRUE ---------------------------
      
      # > estimação de y com densidade condicional t com média e variância variante no tempo
      # theta (parâmetros estáticos): c(w, A0, B1, w, A0, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rt(n+1, df = df))
      
      # valor inicial de mu = 0
      # valor inicial de sigma = 0.1
      f1 <- 0
      f2 <- 0.1
      y <- NULL
      for(t in 1:(n+1)){
        y[t] <- f1[t] + sqrt(exp(f2[t]))*epsilon[t]
        f1[t+1] <- w + A1[1]*( 1/(df + 1) * (y[t] - f1[t])/(1 + (y[t] - f1[t])^2/df) * beta(1/2, df/2)/beta(3/2, (df + 2)/2)) + B1[1]*f1[t]
        f2[t+1] <- w + A1[2]*((df + 3)/df * (((df + 1)* ((y[t] - f1[t])^2/(df*exp(f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(f2[t])))) - 1)) + B1[2]*f2[t]
      }
      
      y <- ts(y[2:(n+1)])
      f1 <- ts(f1[2:(n+1)])
      f2 <- exp(ts(f2[2:(n+1)]))
      epsilon <- ts(epsilon[2:(n+1)])
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, f1, f2, epsilon))
    }
    
  }
    
}