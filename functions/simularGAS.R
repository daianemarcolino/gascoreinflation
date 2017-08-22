simularGAS <- function(density = "normal", n = 100, mean = "variavel", sd = "fixa", seed = 123, theta = NULL){
  
  if(density == "normal"){
    
    if(mean == "variavel" & sd == "fixa"){
      
      # média variável, variância fixa ---------------------------
      
      # > estimação de y com densidade condicional de normal com média variante no tempo e variância fixa no tempo
      # theta (parâmetros estáticos): c(sigma2, w, A0, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n, mean = 0, sd = sqrt(theta[1])))
      
      # valor inicial de mu = 0
      mu <- 0
      y <- NA
      for(t in 2:n){
        y[t] <- mu[t-1] + epsilon[t]
        mu[t] <- theta[2] + theta[3]*(y[t] -  mu[t-1]) + theta[4]*mu[t-1]
      }
      
      y <- ts(y)
      mu <- ts(mu)
      
      # retornar y, mu e epsilon
      return(cbind(y, mu, epsilon))
      
    }else if(mean == "fixa" & sd == "variavel"){
      
      # média fixa, variância variável ---------------------------
      
      # > estimação de y com densidade condicional de normal com média fixa e variância no tempo
      # theta (parâmetros estáticos): c(mean y, w, A0, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n, mean = 0, sd = 1))
      
      # valor inicial de sigma = 0
      sigma2 <- 0
      y <- NA
      for(t in 2:n){
        y[t] <- theta[1] + sqrt(sigma2[t-1])*epsilon[t]
        sigma2[t] <- theta[2] + theta[3]*((y[t] -  theta[1])^2 - sigma2[t-1])+ theta[4]*sigma2[t-1]
      }
      
      y <- ts(y)
      sigma2 <- ts(sigma2)
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, sigma2, epsilon))
      
    }else if(mean == "variavel" & sd == "variavel"){
      
      # média fixa, variância variável ---------------------------
      
      # > estimação de y com densidade condicional de normal com média fixa e variância no tempo
      # theta (parâmetros estáticos): c(w, A0, B1, w, A0, B1)
      
      # epsilon 
      set.seed(seed)
      epsilon <- ts(rnorm(n, mean = 0, sd = 1))
      
      # valor inicial de mu = 0
      # valor inicial de sigma = 0
      mu <- 0
      sigma2 <- 0
      y <- NA
      for(t in 2:n){
        y[t] <- mu[t-1] + sqrt(sigma2[t-1])*epsilon[t]
        mu[t] <- theta[1] + theta[2]*(y[t] -  mu[t-1]) + theta[3]*mu[t-1]
        sigma2[t] <- theta[4] + theta[5]*((y[t] -  mu[t-1])^2 - sigma2[t-1])+ theta[6]*sigma2[t-1]
      }
      
      y <- ts(y)
      mu <- ts(mu)
      sigma2 <- ts(sigma2)
      
      # retornar y, sigma2 e epsilon
      return(cbind(y, mu, sigma2, epsilon))
      
    }
  }
}