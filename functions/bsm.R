bsm <- function(y){
  # Kalman Filter - Basic Structural Model
  # y[t] = mu[t] + gamma[t] + epsilon[t]
  # mu[t+1] = mu[t] + beta[t] + erro1[t] 
  # beta[t+1] = beta[t] + erro2[t] 
  # gamma[t+1] = -sum(gamma[t:(t-s+2)]) + erro3[t] 
  
  # y[t] = Z[t] * alpha[t] + epsilon[t]
  # alpha[t+1] = T[t] * alpha[t] + R[t] * neta[t]
  
  #y = pseudo.y
  
  otimizar <- function(y, par, p = 1, m = 12, r = 2){
    
    N <- length(y)
    
    # matrizes estáticas
    #Z.mu <- matrix(c(1,0), nrow = 1, ncol = 2)
    Z.mu <- matrix(c(1), nrow = 1, ncol = 1)
    Z.gamma <- matrix(c(1,rep(0,10)), nrow = 1, ncol = 11)
    #T.mu <- matrix(c(1,0,1,1), nrow = 2, ncol = 2)
    T.mu <- matrix(c(1), nrow = 1, ncol = 1)
    T.gamma <- cbind(matrix(c(rep(-1,10),diag(10)), nrow = 11, ncol = 10, byrow = T),c(-1,rep(0,10)))
    #R.mu <- diag(2)
    R.mu <- diag(1)
    R.gamma <- matrix(c(1,rep(0,10)))
    
    Zt <- t(c(Z.mu,Z.gamma))
    Tt <- as.matrix(Matrix::bdiag(T.mu,T.gamma))
    Rt <- as.matrix(Matrix::bdiag(R.mu,R.gamma))
    #Qt <- matrix(c(par[1],0,0,0,par[2],0,0,0,par[3]), ncol = r, nrow = r)
    Qt <- matrix(c(par[1],0,0,par[2]), ncol = r, nrow = r)
    
    # matrizes dinâmicas
    alphat <- array(NA, dim = c(m,p,N))
    alphatt <- array(NA, dim = c(m,p,N))
    Ft <- array(NA, dim = c(p,p,N))
    Kt <- array(NA, dim = c(m,p,N))
    Pt <- array(NA, dim = c(m,m,N))
    Ptt <- array(NA, dim = c(m,m,N))
    vt <- matrix(NA, ncol = 1, nrow = N)
    Ht <- matrix(par[3], ncol = p, nrow = p)
    
    # condições iniciais
    alphat[,,1] <- matrix(rep(0,m), ncol = 1)
    Pt[,,1] <- 10e7*diag(m)
    
    # filtro de kalman
    for(t in 1:(N-1)){
      vt[t] <- y[t] - Zt %*% matrix(alphat[,,t], ncol = 1)
      Ft[,,t] <- Zt %*% Pt[,,t] %*% t(Zt) + Ht
      Kt[,,t] <- Tt %*% Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t])
      alphatt[,,t] <- alphat[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% vt[t]
      Ptt[,,t] <- Pt[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% Zt %*% Pt[,,t]
      alphat[,,t+1] <- Tt %*% alphat[,,t] + matrix(Kt[,,t], nrow = m, ncol = p) %*% vt[t]
      Pt[,,t+1] <- Tt %*% Pt[,,t] %*% t(Tt - Kt[,,t] %*% Zt) + Rt %*% Qt %*% t(Rt)
    }
    t <- N
    vt[t] <- y[t] - Zt %*% matrix(alphat[,,t], ncol = 1)
    Ft[,,t] <- Zt %*% Pt[,,t] %*% t(Zt) + Ht
    Kt[,,t] <- Tt %*% Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t])
    alphatt[,,t] <- alphat[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% vt[t]
    Ptt[,,t] <- Pt[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% Zt %*% Pt[,,t]
    
    # transformar em st
    vt <- ts(vt, end = end(y), freq = 12)
    Ft <- ts(c(Ft), end = end(y), freq = 12)
    
    # verossimilhança
    loglik <- -N/2 * log(2*pi) - 0.5 * sum(log(Ft)[-c(1:12)] + vt[-c(1:12)]^2/Ft[-c(1:12)])
    -loglik
  }
  
  # initial = data.frame(param = c("sigma_mu","sigma_beta","sigma_gamma","sigma_y"),
  #                      value = c(0.01, 2, 0.45, 0.05),
  #                      lower = c(0,0,0,0),
  #                      upper = c(Inf,Inf,Inf,Inf))
  initial = data.frame(param = c("sigma_mu","sigma_gamma","sigma_y"),
                       value = c(0.01, 0.45, 0.05),
                       lower = c(0.000001,0.00001,0.00001),
                       upper = c(Inf,Inf,Inf))
  otimizados <- nlminb(start = initial$value, objective = otimizar, y = y, 
                       lower = initial$lower, upper = initial$upper, control = list(eval.max = 10000, iter.max = 10000))
  
  p <- 1
  m <- 12
  r <- 2
  N <- length(y)
  # matrizes estáticas
  #Z.mu <- matrix(c(1,0), nrow = 1, ncol = 2)
  Z.mu <- matrix(c(1), nrow = 1, ncol = 1)
  Z.gamma <- matrix(c(1,rep(0,10)), nrow = 1, ncol = 11)
  #T.mu <- matrix(c(1,0,1,1), nrow = 2, ncol = 2)
  T.mu <- matrix(c(1), nrow = 1, ncol = 1)
  T.gamma <- cbind(matrix(c(rep(-1,10),diag(10)), nrow = 11, ncol = 10, byrow = T),c(-1,rep(0,10)))
  #R.mu <- diag(2)
  R.mu <- diag(1)
  R.gamma <- matrix(c(1,rep(0,10)))
  
  Zt <- t(c(Z.mu,Z.gamma))
  Tt <- as.matrix(Matrix::bdiag(T.mu,T.gamma))
  Rt <- as.matrix(Matrix::bdiag(R.mu,R.gamma))
  #Qt <- matrix(c(par[1],0,0,0,par[2],0,0,0,par[3]), ncol = r, nrow = r)
  Qt <- matrix(c(otimizados$par[1],0,0,otimizados$par[2]), ncol = r, nrow = r)
  
  # matrizes dinâmicas
  alphat <- array(NA, dim = c(m,p,N))
  alphatt <- array(NA, dim = c(m,p,N))
  Ft <- array(NA, dim = c(p,p,N))
  Kt <- array(NA, dim = c(m,p,N))
  Pt <- array(NA, dim = c(m,m,N))
  Ptt <- array(NA, dim = c(m,m,N))
  vt <- matrix(NA, ncol = 1, nrow = N)
  Ht <- matrix(otimizados$par[3], ncol = p, nrow = p)
  
  # condições iniciais
  alphat[,,1] <- matrix(rep(0,m), ncol = 1)
  Pt[,,1] <- 10e7*diag(m)
  
  # filtro de kalman
  for(t in 1:(N-1)){
    vt[t] <- y[t] - Zt %*% matrix(alphat[,,t], ncol = 1)
    Ft[,,t] <- Zt %*% Pt[,,t] %*% t(Zt) + Ht
    Kt[,,t] <- Tt %*% Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t])
    alphatt[,,t] <- alphat[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% vt[t]
    Ptt[,,t] <- Pt[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% Zt %*% Pt[,,t]
    alphat[,,t+1] <- Tt %*% alphat[,,t] + matrix(Kt[,,t], nrow = m, ncol = p) %*% vt[t]
    Pt[,,t+1] <- Tt %*% Pt[,,t] %*% t(Tt - Kt[,,t] %*% Zt) + Rt %*% Qt %*% t(Rt)
  }
  t <- N
  vt[t] <- y[t] - Zt %*% matrix(alphat[,,t], ncol = 1)
  Ft[,,t] <- Zt %*% Pt[,,t] %*% t(Zt) + Ht
  Kt[,,t] <- Tt %*% Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t])
  alphatt[,,t] <- alphat[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% vt[t]
  Ptt[,,t] <- Pt[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% Zt %*% Pt[,,t]
  
  # suavizador de kalman
  rt <- array(NA, dim = c(m,1,N))
  Lt <- array(NA, dim = c(m,m,N))
  alphachapeut <- array(NA, dim = c(m,1,N))
  Nt <- array(NA, dim = c(m,m,N))
  Vt <- array(NA, dim = c(m,m,N))
  
  rt[,,N] <- 0
  Nt[,,N] <- 0
  Lt[,,N] <- Tt - Kt[,,N] %*% Zt
   
  for(t in (N-1):1){
    Lt[,,t] <- Tt - Kt[,,t] %*% Zt
    rt[,,t] <- t(Zt) %*% solve(Ft[,,t+1]) %*% vt[t+1] + t(Lt[,,t+1]) %*% rt[,,t+1]
    Nt[,,t] <- t(Zt) %*% solve(Ft[,,t+1]) %*% Zt + t(Lt[,,t+1]) %*% Nt[,,t+1] %*% Lt[,,t+1]
    alphachapeut[,,t+1] <- alphat[,,t+1] + Pt[,,t+1] %*% rt[,,t]
    Vt[,,t+1] <- Pt[,,t+1] + Pt[,,t+1] %*% Nt[,,t] %*% Pt[,,t+1]
  }

  # transformar em st
  vt <- ts(vt, end = end(y), freq = 12)
  Ft <- ts(c(Ft), end = end(y), freq = 12)
  mu <- ts(c(alphat[1,,]), end = end(y), freq = 12)
  gamma <- ts(c(alphat[2,,]), end = end(y), freq = 12)
  mus <- ts(c(alphachapeut[1,,]), end = end(y), freq = 12)
  gammas <- ts(c(alphachapeut[2,,]), end = end(y), freq = 12)
  
  # verossimilhança
  loglik <- -N/2 * log(2*pi) - 0.5 * sum(log(Ft)[-c(1:12)] + vt[-c(1:12)]^2/Ft[-c(1:12)])
  
  filter <- cbind(mu,gamma)
  smooth <- cbind(mus,gammas)
  colnames(smooth) <- c("mu","gamma")
  
  # output
  list(filter = filter, smooth = smooth)
  
}
