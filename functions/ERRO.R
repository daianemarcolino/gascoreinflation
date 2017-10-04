# medidas de erro

ERRO <- function(y, ychapeu){
  data <- na.omit(cbind(y,ychapeu))
  medidas <- matrix(NA, nrow = 1, ncol = 3)
  colnames(medidas) <- c("RMSE","MAPE","sMAPE")
  medidas[1,"RMSE"] <- sqrt(mean((data[,1] - data[,2])^2))
  medidas[1,"MAPE"] <- mean(abs((data[,1] - data[,2])/data[,1]))*100
  medidas[1,"sMAPE"] <-     (1*nrow(data))*(sum(abs(data[,1] - data[,2]))/sum(abs(data[,1] + data[,2])[-1]))

  #output
  print(medidas)
  invisible(medidas)
}
  