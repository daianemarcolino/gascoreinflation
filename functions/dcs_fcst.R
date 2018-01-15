dcs_fcst <- function(y, start = c(2013,1), initial = NULL, type = "BSM2_beta_psi", outlier = F, m = 1000){
  
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
  
  if(outlier){
    for(i in 1:m){
      for(j in 1:k){
        yt[j,i] <- rt(1, df = out$otimizados$par[5])*tail(out$out[,"sigma"],1) + sum(tail(out$out[,c("mu","gamma","psi")],1))
        novoy <- ts(c(y0,as.vector(yt[1:j,i])), start = start(y0), freq = 12)
        initial0$Dummy <- window(initial$Dummy, end = c(datas$ano[j+1],datas$mes[j+1]), freq = 12)
        out <- dcs_fk_estimation(novoy, initial = initial0, type = type, outlier = outlier, otimo = F)
      }
      message("replication ",i)
    }
  }else{
    for(i in 1:m){
      for(j in 1:k){
        yt[j,i] <- rt(1, df = out$otimizados$par[5])*tail(out$out[,"sigma"],1) + sum(tail(out$out[,c("mu","gamma","psi")],1))
        novoy <- ts(c(y0,as.vector(yt[1:j,i])), start = start(y0), freq = 12)
        out <- dcs_fk_estimation(novoy, initial = initial0, type = type, outlier = outlier, otimo = F)
      }
      message("replication ",i)
    }
  }
  
  # arrumar dados para exportar
  yt <- ts(yt, start = start, freq = 12)
  colnames(yt) <- paste0("iter",1:ncol(yt))
  y_fcst <- ts(rowMeans(yt), start = start, freq = 12)
  y_upper <- ts(apply(yt, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.975)), start = start, freq = 12)
  y_lower <- ts(apply(yt, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.025)), start = start, freq = 12)
  
  # output
  list(replications = yt, fcst = cbind(mean = y_fcst, lower = y_lower, upper = y_upper))
}  

