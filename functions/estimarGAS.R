estimarGAS <- function(y, density = "normal", mean = "variavel", sd = "fixa", initial = NULL){
  
  
  if(density == "normal"){
    
    if(mean == "variavel" & sd == "fixa"){
      
      # > estimação via ML para densidade condicional de y normal com média variante no tempo e variância fixa no tempo
      
      # y[t] = mu[t-1] + epsilon[t], epsilon[t] ~ N(0,sigma2)
      # ln p(y[t]|*) = -0.5*ln(2*pi*sigma2) - 0.5*(y[t] - mu[t-1])^2
      # nabla[t] = (y[t] - mu[t-1])^2/sigma2
      # matriz_inform[t-1] = 1/sigma2
      # st = y[t] - mu[t-1]
      # f[t] = w + A0*(y[t] - mu[t-1]) + B1*mu[t-1]
      
      otimizar <- function(y, par){
        mu <- 0 
        for(t in 2:length(y)){
          mu[t] <- par[2] + par[3]*(y[t] -  mu[t-1]) + par[4]*mu[t-1]
        }
        mu <- ts(mu, start = start(y), freq = frequency(y))
        loglik <- sum(-0.5*log(2*pi*par[1]) - 0.5*(y - lag(mu,-1))^2/par[1], na.rm = T)
        -loglik
      }
      
      
      
      # initial: c(sigma2, w, A0, B1)
      otimizados <- optim(par = initial, fn = otimizar, y = y, method = "Nelder-Mead", hessian = F)
      #lower = c(0,0,0,0), upper = c(Inf,Inf,Inf,Inf), control = list(trace = 0))
      param <- otimizados$par
      
      mu <- 0 
      for(t in 2:length(y)){
        mu[t] <- param[2] + param[3]*(y[t] -  mu[t-1]) + param[4]*mu[t-1]
      }
      
      list(optim = param,
           ft = ts(mu, end = end(y), freq = frequency(y)))
      
    }else if(mean == "fixa" & sd == "variavel"){
      
      # média fixa, variância variável ---------------------------
      
      # > estimação via ML para densidade condicional de y normal com média fixa e variância variante no tempo
      # theta (parâmetros estáticos): c(mean y, w, A0, B1)
      
      # y[t] = mu + sigma[t-1]*epsilon[t], epsilon[t] ~ N(0,1)
      # ln p(y[t]|*) = -0.5*ln(2*pi*sigma2[t-1]) - 0.5*(y[t] - mu)^2/sigma2[t-1]
      # nabla[t] = -0.5/sigma2[t-1] + 0.5*(y[t] - mu)^2/(sigma2[t-1])^2
      # matriz_inform[t-1] = 1/(sigma2[t-1])^2
      # st = (y[t] - mu)^2 - sigma2[t-1]
      # f[t] = w + A0*((y[t] - mu)^2 - sigma2[t-1]) + B1*sigma2[t-1]
      
      otimizar <- function(y, par){
        sigma2 <- 0
        for(t in 2:n){
          sigma2[t] <- par[2] + par[3]*((y[t] - par[1])^2 - sigma2[t-1])+ par[4]*sigma2[t-1]
        }
        sigma2 <- ts(sigma2, start = start(y), freq = frequency(y))
        loglik <- sum(-0.5*log(2*pi*lag(sigma2,-1)) - 0.5*(y - par[1])^2/lag(sigma2,-1), na.rm = T)
        -loglik
      }
      
      # initial: c(sigma2, w, A0, B1)
      otimizados <- optim(par = initial, fn = otimizar, y = y, method = "Nelder-Mead", hessian = F)
      #lower = c(0,0,0,0), upper = c(Inf,Inf,Inf,Inf), control = list(trace = 0))
      param <- otimizados$par
      
      sigma2 <- 0 
      for(t in 2:length(y)){
        sigma2[t] <- param[2] + param[3]*((y[t] -  param[1])^2 - sigma2[t-1])+ param[4]*sigma2[t-1]
      }
      
      list(optim = param,
           ft = ts(sigma2, end = end(y), freq = frequency(y)))
      
    }else if(mean == "variavel" & sd == "variavel"){
      
      # média variável, variância variável ---------------------------
      
      # > estimação via ML para densidade condicional de y normal com média e variância variáveis no tempo
      # theta (parâmetros estáticos): c(w, A0, B1, w, A0, B1)
      
      # y[t] = mu + sigma[t-1]*epsilon[t], epsilon[t] ~ N(0,1)
      # ln p(y[t]|*) = -0.5*ln(2*pi*sigma2[t-1]) - 0.5*(y[t] - mu[t-1])^2/sigma2[t-1]
      # mu: nabla[t] = (y[t] - mu[t-1])/sigma2[t-1] 
      # sigma: nabla[t] = 0.5/sigma2[t-1]*(-1 + (y[t] - mu[t-1])^2/sigma2[t-1])
      
      # mu: matriz_inform[t-1] = 1/(sigma2[t-1])^2
      # sigma: matriz_inform[t-1] = 0.5/(sigma2[t-1])^2
      
      # mu: st = (y[t] - mu[t-1])^2 
      # sigma: st = (y[t] - mu[t-1])^2 - sigma2[t-1]
      
      # mu: f[t] = w + A0*((y[t] - mu[t-1])^2) + B1*mu[t-1]
      # sigma: f[t] = w + A0*((y[t] - mu[t-1])^2 - sigma2[t-1]) + B1*sigma2[t-1]
      
      
      otimizar <- function(y, par){
        mu <- 0
        sigma2 <- 0
        for(t in 2:n){
          mu[t] <- par[1] + par[2]*((y[t] - mu[t-1])^2)+ par[3]*mu[t-1]
          sigma2[t] <- par[4] + par[5]*((y[t] - mu[t-1])^2 - sigma2[t-1]) + par[6]*sigma2[t-1]
        }
        mu <- ts(mu, start = start(y), freq = frequency(y))
        sigma2 <- ts(sigma2, start = start(y), freq = frequency(y))
        loglik <- sum(-0.5*log(2*pi*lag(sigma2,-1)) - 0.5*(y - mu[t-1])^2/lag(sigma2,-1), na.rm = T)
        -loglik
      }
      
      # initial: c(sigma2, w, A0, B1)
      otimizados <- optim(par = initial, fn = otimizar, y = y, method = "Nelder-Mead", hessian = F)
      param <- otimizados$par
      
      mu <- 0
      sigma2 <- 0 
      for(t in 2:length(y)){
        mu[t] <- param[1] + param[2]*((y[t] - mu[t-1])^2) + param[3]*mu[t-1]
        sigma2[t] <- param[4] + param[5]*((y[t] - mu[t-1])^2 - sigma2[t-1]) + param[6]*sigma2[t-1]
      }
      
      list(optim = param,
           ft = ts(cbind(mu,sigma2), end = end(y), freq = frequency(y)))
    }
    
  } # fim if density = normal
  
}

# nabla <- function(y, par){
#   mu <- 0 
#   for(t in 2:n){
#     mu[t] <- par[2] + par[3]*(y[t] -  mu[t-1]) + par[4]*mu[t-1]
#   }
#   mu <- ts(mu)
#   -sum((y - lag(mu,-1))^2/par[1], na.rm = T)
# }