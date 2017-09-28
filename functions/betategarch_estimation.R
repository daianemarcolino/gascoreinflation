betategarchGAS <- function(y, initial = NULL, type = "var"){
  if(type == "var"){
    
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

    # output
    list(f2 = exp(f2), otimizados = otimizados, loglik = -loglik)
  }
}
