betategarchGAS <- function(y, initial = NULL, type = "mean-var"){
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
    
    ci <- nlminb(start = c(10,1), objective = inicializar, y = y,
                 lower = c(2,-Inf), upper = c(Inf,Inf))
    
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
    
    otimizados <- nlminb(start = c(initial, ci$par), objective = otimizar, y = y, ci = ci,
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
  }else if(type == "mean-var"){
    
    # > estimação via ML para densidade condicional de y t-student com média e variância variantes no tempo
    inicializar <- function(y, par){
      df <- par[1]
      N <- length(y)
      f1 <- par[2]
      f2 <- par[3]
      # loglik
      p1 <- N*log(gamma((df+1)/2)/(sqrt(pi*df)*gamma(df/2)))
      p2 <- -((df+1)/2)*sum(log((1 + (y-f1)^2/(df*exp(2*f2)))))
      p3 <- N/2*log((df-2)/df) - N*f2
      
      loglik <- p1 + p2 + p3
      -loglik
    }
    
    ci <- nlminb(start = c(10,0,1), objective = inicializar, y = y,
                 lower = c(2,-Inf,-Inf), upper = c(Inf,Inf,Inf))
    
    otimizar <- function(y, par, ci){
      w1 <- par[1]
      A1 <- par[2]
      B1 <- par[3]
      w2 <- par[4]
      A2 <- par[5]
      B2 <- par[6]
      df <- par[7]#ci$par[1]
      N <- length(y)
      f1 <- par[8]#ci$par[2]
      f2 <- par[9]#ci$par[3]
      for(t in 1:N){
        f1[t+1] <- w1 + A1*( 1/(df + 1) * (y[t] - f1[t])/(1 + (y[t] - f1[t])^2/(df*exp(2*f2[t]))) * beta(1/2, df/2)/beta(3/2, (df + 2)/2)) + B1*f1[t]
        f2[t+1] <- w2 + A2*((df + 3)/(2*df) * (((df + 1)*((y[t] - f1[t])^2/(df*exp(2*f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(2*f2[t])))) - 1)) + B2*f2[t]
      }
      f1 <- ts(f1[-1], start = start(y), freq = frequency(y))
      f2 <- ts(f2[-1], start = start(y), freq = frequency(y))
      
      # loglik
      p1 <- N*log(gamma((df+1)/2)/(sqrt(pi*df)*gamma(df/2)))
      p2 <- -(df+1)/2*sum(log((1 + (y-f1)^2/(df*exp(2*f2)))))
      p3 <- N/2*log((df-2)/df) - sum(f2)
      
      
      loglik <- p1 + p2 + p3
      -loglik
    }
    
    otimizados <- nlminb(start = c(0.02,0.05,0.95,0.2,0.5,0.5,ci$par), objective = otimizar, y = y, #ci = ci,
                         lower = c(-Inf,-Inf,-1,-Inf,-Inf,-1), upper = c(Inf,Inf,1,Inf,Inf,1))
    
    w1 <- otimizados$par[1]
    A1 <- otimizados$par[2]
    B1 <- otimizados$par[3]
    w2 <- otimizados$par[4]
    A2 <- otimizados$par[5]
    B2 <- otimizados$par[6]
    
    df <- otimizados$par[7]
    f1 <- otimizados$par[8]
    f2 <- otimizados$par[9]
    N <- length(y)
    
    for(t in 1:N){
      f1[t+1] <- w1 + A1*( 1/(df+1)* beta(0.5,df/2)/beta(3/2,(df+2)/2)*(y[t] - f1[t])/(1 + (y[t] - f1[t])^2/(df*exp(f2[t])))) + B1*f1[t]
      f2[t+1] <- w1 + A1*((df + 3)/df * (((df + 1)*((y[t] - f1[t])^2/(df*exp(f2[t])))/(1 + (y[t] - f1[t])^2/(df*exp(f2[t])))) - 1)) + B1*f2[t]
    }
    
    f1 <- ts(f1[-1], start = start(y), freq = frequency(y))
    f2 <- ts(exp(f2[-1]), start = start(y), freq = frequency(y))
    list(cbind(f1,f2), otimizados, ci)
    
    
  }
  
}

# otimizar <- function(y, par){
#   w1 <- par[1]
#   #w2 <- par[2]
#   A1 <- par[2]
#   #A2 <- par[4]
#   B1 <- par[3]
#   #B2 <- par[6]
#   df <- par[4]
#   
#   N <- length(y)
#   #f1 <- 0.001
#   f2 <- par[5]
#   for(t in 1:N){
#     #f1[t+1] <- w1 + A1*( 1/(df+1)* beta(0.5,df/2)/beta(3/2,(df+2)/2)*(y[t] - f1[t])/(1 + (y[t] - f1[t])^2/(df*exp(f2[t])))) + B1*f1[t]
#     f2[t+1] <- w1 + A1*((df + 3)/df * (((df + 1)*((y[t] - 0)^2/(df*exp(f2[t])))/(1 + (y[t] - 0)^2/(df*exp(f2[t])))) - 1)) + B1*f2[t]
#   }
#   
#   #f1 <- ts(f1[-1], start = start(y), freq = frequency(y))
#   f2 <- ts(f2[-1], start = start(y), freq = frequency(y))
#   
#   # loglik
#   p1 <- N*log( gamma((df+1)/2)/(sqrt(pi*df)*gamma(df/2)))
#   p2 <- -(df+1)/2*sum(1 + (y - 0)^2/(df*exp(f2)))
#   p3 <- N*log(sqrt((df-2)/df)) - sum(f2/2)
#   
#   loglik <- p1 + p2 + p3
#   -loglik
# }