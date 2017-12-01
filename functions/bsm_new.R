bsm <- function(y){
  # Kalman Filter - Basic Structural Model
  # y[t] = mu[t] + gamma[t] + epsilon[t]
  # mu[t+1] = mu[t] + beta[t] + erro1[t] 
  # beta[t+1] = beta[t] + erro2[t] 
  # gamma[t+1] = -sum(gamma[t:(t-s+2)]) + erro3[t] 
  
  # y[t] = Z[t] * alpha[t] + epsilon[t]
  # alpha[t+1] = T[t] * alpha[t] + R[t] * neta[t]
  
  #y = pseudo.y
  
  # modelo:
  # y[t] = mu[t] + gamma[t] + psi[t] + epsilon[t], epsilon[t] ~ t(v)
  # mu[t+1] = beta + mu[t] + u1[t]
  # gamma[t+1] = gamma[t] + u2[t] 
  # psi[t+1] = phi*psi[t] + u3[t]
  # > estimação via ML para densidade condicional de y t-student com média constante e variância constante no tempo
  
  


  otimizar <- function(y, par){
    
    # N <- length(y)
    # k1 <- par[1]
    # ks <- par[2]
    # f2 <- par[3]
    # df <- par[4]
    # mu <- par[5]
    # beta <- par[6]
    # psi <- par[7]
    # phi <- par[8]
    # k3 <- par[9]
    # alpha <- matrix(NA, ncol = 12, nrow = N)
    # #alpha[1,] <- c(initial_gamma, - sum(initial_gamma))  
    # alpha[1,] <- c(par[10:20], - sum(par[10:20]))  
    # if(outlier){
    #   d <- par[21:(21+ncol(data.frame(Dummy))-1)]
    # }else{
    #   d <- NA
    # }
    # j <- cycle(y)
    # gamma <- alpha[1,j[1]]
    # u1 <- NULL
    # 
    # if(outlier){
    #   for(t in 1:(N-1)){ 
    #     u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))/(1 + (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))^2/(df*exp(2*f2)))
    #     mu[t+1] <- beta + mu[t] + k1*u1[t] 
    #     k <- rep(-ks/11,12) 
    #     k[j[t]] <- ks
    #     alpha[t+1,] <- alpha[t,] + k*u1[t]
    #     gamma[t+1] <- alpha[t+1,j[t+1]]
    #     psi[t+1] <- phi*psi[t] + k3*u1[t]
    #   } 
    #   t <- N
    #   u1[t] <- (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))/(1 + (y[t] - mu[t] - gamma[t] - psi[t] - as.vector(t(d) %*% as.matrix(Dummy)[t,]))^2/(df*exp(2*f2)))
    #   
    #   mu <- ts(mu, end = end(y), freq = frequency(y))
    #   gamma <- ts(gamma, end = end(y), freq = frequency(y))
    #   psi <- ts(psi, end = end(y), freq = frequency(y))
    #   dum <- ts(t(d %*% t(Dummy)), end = end(y), freq = frequency(y))
    #   
    #   # loglik
    #   loglik <- N*log(gamma((df+1)/2)) - (N/2)*log(pi) - N*log(gamma(df/2)) - (N/2)*log(df) - N*(f2) - ((df + 1)/2)*sum(log(1 + (y - mu - gamma - psi - dum)^2/(df*exp(2*f2))))
    #   -loglik
      
    N <- length(y)
    p <- 1
    m <- 13
    r <- 3
    #par <- rep(0.2,8)
    
    # matrizes
    Z.mu <- matrix(c(1), nrow = 1, ncol = 1)
    Z.gamma <- matrix(c(1,rep(0,10)), nrow = 1, ncol = 11)
    Z.psi <- matrix(c(1), nrow = 1, ncol = 1)
    Zt <- t(c(Z.mu,Z.gamma, Z.psi))
    
    T.mu <- matrix(c(1), nrow = 1, ncol = 1)
    T.gamma <- cbind(matrix(c(rep(-1,10),diag(10)), nrow = 11, ncol = 10, byrow = T),c(-1,rep(0,10)))
    T.psi <- matrix(par[1], nrow = 1, ncol = 1, byrow = T)
    Tt <- as.matrix(Matrix::bdiag(T.mu,T.gamma, T.psi))
    
    # At <- matrix(c(1,rep(0,12)), ncol = 1)
    At <- matrix(c(par[2],rep(0,12)), ncol = 1)
    
    R.mu <- matrix(c(1,rep(0,12)), ncol = 1)
    R.gamma <- matrix(c(0,1,rep(0,11)))
    R.psi <- matrix(c(rep(0,12), 1), ncol = 1)
    Rt <- cbind(R.mu, R.gamma, R.psi)
    
    #Qt <- matrix(c(1,0,0,0,1,0,0,0,1), ncol = 3, nrow = 3)
    Qt <- matrix(c((par[3]),0,0,0,(par[4]),0,0,0,(par[5])), ncol = r, nrow = r)

    # matrizes dinâmicas
    alphat <- array(NA, dim = c(m,p,N))
    alphatt <- array(NA, dim = c(m,p,N))
    Ft <- array(NA, dim = c(p,p,N))
    Kt <- array(NA, dim = c(m,p,N))
    Pt <- array(NA, dim = c(m,m,N))
    Ptt <- array(NA, dim = c(m,m,N))
    vt <- matrix(NA, ncol = 1, nrow = N)
    Ht <- matrix(par[6], ncol = p, nrow = p)
    
    # condições iniciais
    alphat[,,1] <- matrix(c(par[2],rep(0,m-1)), ncol = 1)
    #alphat[,,1] <- matrix(par[7:(7+m-1)], ncol = 1)
    Pt[,,1] <- 10e7*diag(m)
    
    # filtro de kalman
    for(t in 1:(N-1)){
      vt[t] <- y[t] - Zt %*% matrix(alphat[,,t], ncol = 1)
      Ft[,,t] <- Zt %*% Pt[,,t] %*% t(Zt) + Ht
      Kt[,,t] <- Tt %*% Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t])
      alphatt[,,t] <- alphat[,,t] + Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% vt[t]
      Ptt[,,t] <- Pt[,,t] - Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% Zt %*% Pt[,,t]
      alphat[,,t+1] <- At + Tt %*% alphat[,,t] + matrix(Kt[,,t], nrow = m, ncol = p) %*% vt[t]
      Pt[,,t+1] <- Tt %*% Pt[,,t] %*% t(Tt - Kt[,,t] %*% Zt) + Rt %*% Qt %*% t(Rt)
    }
    t <- N
    vt[t] <- y[t] - Zt %*% matrix(alphat[,,t], ncol = 1)
    Ft[,,t] <- Zt %*% Pt[,,t] %*% t(Zt) + Ht
    Kt[,,t] <- Tt %*% Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t])
    alphatt[,,t] <- alphat[,,t] - Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% vt[t]
    Ptt[,,t] <- Pt[,,t] - Pt[,,t] %*% t(Zt) %*% solve(Ft[,,t]) %*% Zt %*% Pt[,,t]
    
    # transformar em st
    vt <- ts(vt, end = end(y), freq = 12)
    Ft <- ts(c(Ft), end = end(y), freq = 12)
    print(par)
    # verossimilhança
    loglik <- -N/2 * log(2*pi) - 0.5 * sum(log(abs(Ft)) + vt^2/Ft)
    -loglik
  }
  # c(0.561017902,-0.040152044, 0.130667845, 0.080210650,-0.044419827,-0.244535510,-0.120733417,-0.206155100,-0.247965665,-0.100009180, 0.074202254, 1.457084370, 1.883646757)
  # c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11), "d1", "d2")
  # 
  # initial = data.frame(param = c("phi","beta","sigma_mu","sigma_gamma","sigma_psi","sigma_y", paste0("alpha",1:13,"[0]")),
  #                      value = c(0.24,0.5,0.01,0.2,0.2,0.2,0.68,0.56,-0.04,0.13,0.08,-0.04,-0.24,-0.12,-0.20,-0.25,-0.1,0.07,-0.72),
  #                      lower = c(-1,-Inf,  0,  0,  0,  0,Inf, rep(-Inf,11),-Inf),
  #                      upper = c( 1, Inf,Inf,Inf,Inf,Inf,Inf,rep(Inf,11),Inf))
  initial = data.frame(param = c("phi","beta","sigma_mu","sigma_gamma","sigma_psi","sigma_y"),
                       value = c(0.05,0.2,1,1,1,0.3),
                       #value = c(0.0398239491, 0.0013323167, 0.1231947059, 0.0018684698, 0.0006363430, 0.0004379238),
                       lower = c(-1,-Inf,  0,  0,  0,  0),
                       upper = c( 1, Inf,Inf,Inf,Inf,Inf))
  otimizados <- nlminb(start = initial$value, objective = otimizar, y = y, 
                       lower = initial$lower, upper = initial$upper, control = list(eval.max = 10000, iter.max = 10000))
  
  otimizados
  N <- length(y)
  p <- 1
  m <- 13
  r <- 3
  #par <- rep(0.2,8)
  
  # matrizes
  Z.mu <- matrix(c(1), nrow = 1, ncol = 1)
  Z.gamma <- matrix(c(1,rep(0,10)), nrow = 1, ncol = 11)
  Z.psi <- matrix(c(1), nrow = 1, ncol = 1)
  Zt <- t(c(Z.mu,Z.gamma, Z.psi))
  
  # T.mu <- matrix(c(1), nrow = 1, ncol = 1)
  # T.gamma <- cbind(matrix(c(rep(-1,10),diag(10)), nrow = 11, ncol = 10, byrow = T),c(-1,rep(0,10)))
  # T.psi <- matrix(1, nrow = 1, ncol = 1, byrow = T)
  T.mu <- matrix(c(1), nrow = 1, ncol = 1)
  T.gamma <- cbind(matrix(c(rep(-1,10),diag(10)), nrow = 11, ncol = 10, byrow = T),c(-1,rep(0,10)))
  T.psi <- matrix(otimizados$par[1], nrow = 1, ncol = 1, byrow = T)
  Tt <- as.matrix(Matrix::bdiag(T.mu,T.gamma, T.psi))
  
  # At <- matrix(c(1,rep(0,12)), ncol = 1)
  At <- matrix(c(otimizados$par[2],rep(0,12)), ncol = 1)
  
  R.mu <- matrix(c(1,rep(0,12)), ncol = 1)
  R.gamma <- matrix(c(0,1,rep(0,11)))
  R.psi <- matrix(c(rep(0,12), 1), ncol = 1)
  Rt <- cbind(R.mu, R.gamma, R.psi)
  
  #Qt <- matrix(c(1,0,0,0,1,0,0,0,1), ncol = 3, nrow = 3)
  Qt <- matrix(c((otimizados$par[3]),0,0,0,(otimizados$par[4]),0,0,0,(otimizados$par[5])), ncol = r, nrow = r)
  
  # matrizes dinâmicas
  alphat <- array(NA, dim = c(m,p,N))
  alphatt <- array(NA, dim = c(m,p,N))
  Ft <- array(NA, dim = c(p,p,N))
  Kt <- array(NA, dim = c(m,p,N))
  Pt <- array(NA, dim = c(m,m,N))
  Ptt <- array(NA, dim = c(m,m,N))
  vt <- matrix(NA, ncol = 1, nrow = N)
  Ht <- matrix(otimizados$par[6], ncol = p, nrow = p)
  
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
    alphat[,,t+1] <- At + Tt %*% alphat[,,t] + matrix(Kt[,,t], nrow = m, ncol = p) %*% vt[t]
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
  psi <- ts(c(alphat[13,,]), end = end(y), freq = 12)
  mus <- ts(c(alphachapeut[1,,]), end = end(y), freq = 12)
  gammas <- ts(c(alphachapeut[2,,]), end = end(y), freq = 12)
  psis <- ts(c(alphachapeut[13,,]), end = end(y), freq = 12)
  
  # verossimilhança
  loglik <- -N/2 * log(2*pi) - 0.5 * sum(log(Ft)[-c(1:12)] + vt[-c(1:m)]^2/Ft[-c(1:m)])
  
  filter <- cbind(mu,gamma)
  smooth <- cbind(mus,gammas)
  colnames(smooth) <- c("mu","gamma")
  
  # output
  list(filter = filter, smooth = smooth)
  
}
