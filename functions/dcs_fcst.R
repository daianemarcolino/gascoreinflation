dcs_fcst <- function(y, out = F, start = c(2018,1), h = 12, initial = NULL, type = "BSM3", outlier = F, m = 2000, dist = "t"){
  
  # out: previsão dentro ou fora da amostra
  # start: inicio da previsao dentro da amostra (out = F)
  # h: número de passos para prever fora da amostra (out = T)
  # initial: parametros usados na função de otimização
  # type: tipo do modelo DCS (ver dcs_fk_estimation)
  # outler: se tem outlier no modelo
  # m: número de replicações para a previsão, default = 2000
  # dist: distribuição do erro
  
  if(!out){
    
    # pos para cortar e estimar
    pos <- which(start[1] == as.numeric(substr(as.Date(y),1,4)) & start[2] == as.numeric(substr(as.Date(y),6,7))) - 1
    
    datas <- list(ano = as.numeric(substr(as.Date(y),1,4))[pos:length(y)],
                  mes = as.numeric(substr(as.Date(y),6,7))[pos:length(y)])
    
    # window y
    y0 <- window(y, end = c(datas$ano[1],datas$mes[1]), freq = 12)
    
    # parametros iniciais
    initial0 <- initial
    if(outlier){
      initial0$Dummy <- window(initial$Dummy, end = c(datas$ano[1],datas$mes[1]), freq = 12)
    }
    
    # primeira e única otimização
    out <- dcs_fk_estimation(y0, initial = initial0, type = type, outlier = outlier, otimo = T)
    
    # guardar parâmetros otimizados
    initial0$par$value <- out$otimizados$par
    
    # definir o número de passos a frente e criar matriz de armazenamento
    k <- length(datas$ano) - 1
    yt <- matrix(NA, nrow = k, ncol = m)
    mut <- matrix(NA, nrow = k, ncol = m)
    
    if(outlier){
      for(i in 1:m){
        out0 <- out
        for(j in 1:k){
          yt[j,i] <- ifelse(dist == "norm", rnorm(1, mean = sum(tail(out0$out[,c("psi","mu","gamma")],1)), sd = tail(out0$out[,"sigma"],1)),
                            sum(tail(out0$out[,c("psi","mu","gamma")],1)) +  tail(out0$out[,"sigma"],1)*rt(1, df = out0$otimizados$par[5]))
          mut[j,i] <- tail(out0$out[,"mu"],1)
          novoy <- ts(c(y0,as.vector(yt[1:j,i])), start = start(y0), freq = 12)
          initial0$Dummy <- window(initial$Dummy, end = c(datas$ano[j+1],datas$mes[j+1]), freq = 12)
          out0 <- dcs_fk_estimation(novoy, initial = initial0, type = type, outlier = outlier, otimo = F)
        }
        message("replication ",i)
      }
    }else{
      for(i in 1:m){
        out0 <- out
        for(j in 1:k){
          yt[j,i] <-  ifelse(dist == "norm", rnorm(1, mean = sum(tail(out0$out[,c("psi","mu","gamma")],1)), sd = tail(out0$out[,"sigma"],1)),
                             sum(tail(out0$out[,c("psi","mu","gamma")],1)) +  tail(out0$out[,"sigma"],1)*rt(1, df = out0$otimizados$par[5]))
          mut[j,i] <- tail(out0$out[,"mu"],1)
          novoy <- ts(c(y0,as.vector(yt[1:j,i])), start = start(y0), freq = 12)
          out0 <- dcs_fk_estimation(novoy, initial = initial0, type = type, outlier = outlier, otimo = F)
        }
        message("replication ",i)
      }
    }
    
    # arrumar dados para exportar
    yt <- ts(yt, start = start, freq = 12)
    colnames(yt) <- paste0("iter",1:ncol(yt))
    mut <- ts(mut, start = start, freq = 12)
    colnames(mut) <- paste0("iter",1:ncol(mut))
    
    y_fcst <- ts(rowMeans(yt), start = start, freq = 12)
    y_upper <- ts(apply(yt, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.975)), start = start, freq = 12)
    y_lower <- ts(apply(yt, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.025)), start = start, freq = 12)
    
    mu_fcst <- ts(rowMeans(mut), start = start, freq = 12)
    mu_upper <- ts(apply(mut, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.975)), start = start, freq = 12)
    mu_lower <- ts(apply(mut, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.025)), start = start, freq = 12)
    
    
    # output
    list(replications_y = yt, fcst_y = cbind(mean = y_fcst, lower = y_lower, upper = y_upper),
         replications_mu = mut, fcst_mu = cbind(mean = mu_fcst, lower = mu_lower, upper = mu_upper))
    
  }else{
    
    # parametros iniciais
    initial0 <- initial
    
    # primeira e única otimização
    out <- dcs_fk_estimation(y, initial = initial0, type = type, outlier = outlier, otimo = T)
    
    # guardar parâmetros otimizados
    initial0$par$value <- out$otimizados$par
    
    # definir o número de passos a frente e criar matriz de armazenamento
    yt <- matrix(NA, nrow = h, ncol = m)
    mut <- matrix(NA, nrow = h, ncol = m)
    
    if(outlier){
      for(i in 1:m){
        out0 <- out
        for(j in 1:h){
          yt[j,i] <- ifelse(dist == "norm", rnorm(1, mean = sum(tail(out0$out[,c("psi","mu","gamma")],1)), sd = tail(out0$out[,"sigma"],1)),
                            sum(tail(out0$out[,c("psi","mu","gamma")],1)) +  tail(out0$out[,"sigma"],1)*rt(1, df = out0$otimizados$par[5]))
          mut[j,i] <- tail(out0$out[,"mu"],1)
          novoy <- ts(c(y,as.vector(yt[1:j,i])), start = start(y), freq = 12)
          initial0$Dummy <- ts(c(initial$Dummy,rep(0,j)), start = start(initial0$Dummy), freq = 12)
          out0 <- dcs_fk_estimation(novoy, initial = initial0, type = type, outlier = outlier, otimo = F)
        }
        message("replication ",i)
      }
    }else{
      for(i in 1:m){
        out0 <- out
        for(j in 1:k){
          yt[j,i] <-  ifelse(dist == "norm", rnorm(1, mean = sum(tail(out0$out[,c("psi","mu","gamma")],1)), sd = tail(out0$out[,"sigma"],1)),
                             sum(tail(out0$out[,c("psi","mu","gamma")],1)) +  tail(out0$out[,"sigma"],1)*rt(1, df = out0$otimizados$par[5]))
          mut[j,i] <- tail(out0$out[,"mu"],1)
          novoy <- ts(c(y,as.vector(yt[1:j,i])), start = start(y), freq = 12)
          out0 <- dcs_fk_estimation(novoy, initial = initial0, type = type, outlier = outlier, otimo = F)
        }
        message("replication ",i)
      }
    }
    
    # arrumar dados para exportar
    yt <- ts(yt, start = start, freq = 12)
    colnames(yt) <- paste0("iter",1:ncol(yt))
    mut <- ts(mut, start = start, freq = 12)
    colnames(mut) <- paste0("iter",1:ncol(mut))
    
    y_fcst <- ts(rowMeans(yt), start = start, freq = 12)
    y_upper <- ts(apply(yt, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.975)), start = start, freq = 12)
    y_lower <- ts(apply(yt, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.025)), start = start, freq = 12)
    
    mu_fcst <- ts(rowMeans(mut), start = start, freq = 12)
    mu_upper <- ts(apply(mut, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.975)), start = start, freq = 12)
    mu_lower <- ts(apply(mut, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.025)), start = start, freq = 12)
    
    
    # output
    list(replications_y = yt, fcst_y = cbind(mean = y_fcst, lower = y_lower, upper = y_upper),
         replications_mu = mut, fcst_mu = cbind(mean = mu_fcst, lower = mu_lower, upper = mu_upper))
    
  }
}  

