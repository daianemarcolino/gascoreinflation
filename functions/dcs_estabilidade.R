dcs_estabilidade <- function(y, start = c(2013,1), initial = NULL, type = "BSM2_beta_psi", outlier = F){

  # pos para cortar e estimar
  pos <- which(start[1] == as.numeric(substr(as.Date(y),1,4)) & start[2] == as.numeric(substr(as.Date(y),6,7)))
  
  datas <- list(ano = as.numeric(substr(as.Date(y),1,4))[pos:length(y)],
                mes = as.numeric(substr(as.Date(y),6,7))[pos:length(y)])
  
  initial0 <- initial
  core_month <- ts(NA, start = start, end = end(y), freq = 12)
  core_full <- ts(NA, start = start, end = end(y), freq = 12)
  
  for(i in 1:length(datas$ano)){
    # window y
    y0 <- window(y, end = c(datas$ano[i],datas$mes[i]), freq = 12)
    # parametros iniciais
    initial0$Dummy <- window(initial$Dummy, end = c(datas$ano[i],datas$mes[i]), freq = 12)
    # otimização
    out <- dcs_fk_estimation(y0, initial = initial0, type = type, outlier = outlier, otimo = T)
    # guardar parâmetros otimizados
    #initial0$par$value <- out$otimizados$par
    # guardar núcleo no tempo t
    core_month[i] <- tail(out$out[,"mu"],2)[1]
    message("- date: ", datas$ano[i],"/",datas$mes[i])
  }
  
  core_full <- window(out$out[,"mu"], start = start(core_full), end = end(core_full), freq = 12)
  
  # output
  list(out = cbind(core_month, core_full))
}  
 
 