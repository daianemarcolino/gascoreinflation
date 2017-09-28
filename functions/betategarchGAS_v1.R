betategarchGAS <- function(y, initial.1 = NULL, initial.2 = NULL, type = "var"){
  if(type == "var"){
    
    # > estimação via ML para densidade condicional de y t-student com média constante e variância variante no tempo
    inicializar <- function(y, par){
      df <- par[1]
      f1 <- 0
      f2 <- par[2]
      N <- length(y)
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)/(sqrt(pi*df)*gamma(df/2))) - ((df+1)/2)*sum(log((1 + (y - f1)^2/(df*exp(2*f2))))) + (N/2)*log((df-2)/df) - N*f2
      -loglik
    }
    
    ci <- nlminb(start = initial.1, objective = inicializar, y = y,
                 lower = c(3,-Inf), upper = c(Inf,Inf))
    
    otimizar <- function(y, par, ci){
      w1 <- par[1]
      A1 <- par[2]
      B1 <- par[3]
      df <- par[4]
      f1 <- 0
      f2 <- ci$par[2]
      N <- length(y)
      
      for(t in 1:N){
        f2[t+1] <- w1 + A1*((df + 3)/(2*df) * (((df + 1)*((y[t] - f1)^2/(df*exp(2*f2[t])))/(1 + (y[t] - f1)^2/(df*exp(2*f2[t])))) - 1)) + B1*f2[t]
      }
      
      f2 <- ts(f2[-1], start = start(y), freq = frequency(y))
      
      # loglik
      loglik <- N*log(gamma((df+1)/2)/(sqrt(pi*df)*gamma(df/2))) - ((df+1)/2)*sum(log((1 + (y - f1)^2/(df*exp(2*f2))))) + (N/2)*log((df-2)/df) - sum(f2)
      -loglik
    }
    
    otimizados <- nlminb(start = c(initial.2, ci$par), objective = otimizar, y = y, ci = ci,
                         lower = c(-Inf,-Inf,-1,2), upper = c(Inf,Inf,1,Inf), control = list(eval.max = 10000, iter.max = 10000))
    
    # otimizados2 <- optim(par = c(initial, ci$par[1]), fn = otimizar, y = y, ci = ci,
    #                      control = list(maxit = 10000), method = "Nelder-Mead")
    # otimizados3 <- optim(par = c(initial, ci$par[1]), fn = otimizar, y = y, ci = ci,
    #                     control = list(maxit = 10000), method = "BFGS")
    
    
    w1 <- otimizados$par[1]
    A1 <- otimizados$par[2]
    B1 <- otimizados$par[3]
    df <- otimizados$par[4]
    f1 <- 0
    f2 <- ci$par[2]
    N <- length(y)
    
    for(t in 1:N){
      f2[t+1] <- w1 + A1*((df + 3)/(2*df) * (((df + 1)*((y[t] - f1)^2/(df*exp(2*f2[t])))/(1 + (y[t] - f1)^2/(df*exp(2*f2[t])))) - 1)) + B1*f2[t]
    }
    
    f2 <- ts(exp(2*f2[-1]), start = start(y), freq = frequency(y))
    list(f2, otimizados, ci)
  }
}