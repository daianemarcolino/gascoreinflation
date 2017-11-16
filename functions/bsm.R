bsm <- function(y, beta = T, iter = 1){
  # Kalman Filter - Basic Structural Model
  # y[t] = mu[t] + gamma[t] + epsilon[t]
  # mu[t+1] = mu[t] + beta[t] + erro1[t] 
  # beta[t+1] = beta[t] + erro2[t] 
  # gamma[t+1] = -sum(gamma[t:(t-s+2)]) + erro3[t] 
  
  # y[t] = Z[t] * alpha[t] + epsilon[t]
  # alpha[t+1] = T[t] * alpha[t] + R[t] * neta[t]
  
  #y = pseudo.y
  
  dat <- t(as.matrix(y))
  n <- dim(dat)[1]
  TT <- dim(dat)[2]
  
  if(beta){
    m <- 13
    Z <- matrix(c(1,0,1,rep(0,10)), nrow = n, ncol = m, byrow = TRUE)
    B <- matrix(c(1,1,rep(0,11),
                  0,1,rep(0,11),
                  0,0,rep(-1,11),
                  0,0,1,rep(0,10),
                  0,0,0,1,rep(0,9),
                  0,0,0,0,1,rep(0,8),
                  0,0,0,0,0,1,rep(0,7),
                  0,0,0,0,0,0,1,rep(0,6),
                  0,0,0,0,0,0,0,1,rep(0,5),
                  0,0,0,0,0,0,0,0,1,rep(0,4),
                  0,0,0,0,0,0,0,0,0,1,rep(0,3),
                  0,0,0,0,0,0,0,0,0,0,1,rep(0,2),
                  0,0,0,0,0,0,0,0,0,0,0,1,rep(0,1)),
                nrow = m, ncol = m, byrow = TRUE)
    Q <- matrix(list(0), nrow = m, ncol = m, byrow = TRUE)
    diag(Q) <- c("q_trend", "q_slope", rep("q_season",11))
    R <- matrix("r_y", nrow = n, ncol = n, byrow = TRUE)
    V0 <- diag(10^7,m)
    x0 = A = U <- "zero"
  }else{
    m <- 12
    Z <- matrix(c(1,1,rep(0,10)), nrow = n, ncol = m, byrow = TRUE)
    B <- matrix(c(1,rep(0,11),
                  0,rep(-1,11),
                  0,1,rep(0,10),
                  0,0,1,rep(0,9),
                  0,0,0,1,rep(0,8),
                  0,0,0,0,1,rep(0,7),
                  0,0,0,0,0,1,rep(0,6),
                  0,0,0,0,0,0,1,rep(0,5),
                  0,0,0,0,0,0,0,1,rep(0,4),
                  0,0,0,0,0,0,0,0,1,rep(0,3),
                  0,0,0,0,0,0,0,0,0,1,rep(0,2),
                  0,0,0,0,0,0,0,0,0,0,1,rep(0,1)),
                nrow = m, ncol = m, byrow = TRUE)
    Q <- matrix(list(0), nrow = m, ncol = m, byrow = TRUE)
    diag(Q) <- c("q_trend", rep("q_season",11))
    R <- matrix("r", nrow = n, ncol = n, byrow = TRUE)
    V0 <- diag(10^7,m)
    x0 = A = U <- "zero"
  }
  
  dfa.model <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0)
  cntl.list <- list(maxit = iter)
  
  model <- MARSS(dat, model = dfa.model, control = cntl.list, method = "BFGS", fun.kf = "MARSSkfss")
  
  if(beta){
    estados <- ts(t(print(model, what = "xtT", silent = T)), start = start(y), freq = 12)[,1:3]
    colnames(estados) <- c("mu","beta","gamma")
  }else{
    estados <- ts(t(print(model, what = "xtT", silent = T)), start = start(y), freq = 12)[,1:2]
    colnames(estados) <- c("mu","gamma")
  }
  
  #output
  return(estados)  
}
