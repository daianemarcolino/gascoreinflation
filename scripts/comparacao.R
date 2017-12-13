core <- dcs_padrao_psi_dummy$out[,"mu"]
core12 <- ((dcs_padrao_psi_dummy$out[,"mu"]/100+1)^12-1)*100

ts.plot(ipc, core, ipc_ma2, col = c(1,2,1), lwd = c(1,2,2), lty = c(3,1,1))
legend("topleft", legend = c("Novo","Triplo filtro"), lty = c(1), lwd = 2, col =2:1, bty = "n")

ts.plot(ipc, geom3(core), ipc_ma3, col = c(1,2,1), lwd = c(1,2,2), lty = c(3,1,1))
legend("topleft", legend = c("Novo","Triplo filtro"), lty = c(1), lwd = 2, col =2:1, bty = "n")

ts.plot(ipc, core, ipc_ma3, col = c(1,2,1), lwd = c(1,2,2), lty = c(3,1,1))
legend("topleft", legend = c("Novo","Triplo filtro"), lty = c(1), lwd = 2, col =2:1, bty = "n")

pseudo.y <- (1 - dcs_padrao_psi_dummy$out[,"b"])*ipc + dcs_padrao_psi_dummy$out[,"b"]*ts(rowSums(dcs_padrao_psi_dummy$out[,c("mu","gamma","psi","dummy")]), start = start(dcs_padrao_psi_dummy$out), freq = 12)
ts.plot(ipc,pseudo.y, col = 1:2)
core3 <- bsm(pseudo.y,type = "BSM2_beta", iter = 100)
ts.plot(ipc,core,core3[,"mu"], col = c(1:2,4))



geom3 <- function(x, anual = F){ 
  # x: ts 
  n <- length(x)
  data <- data.frame(x, y = rep(NA, n), w = rep(NA, n), z = rep(NA, n), row.names = 1:n)
  colnames(data)[1] <- "x"
  data$y <- data$x/100 + 1
  
  for(i in n:3){
    data[i,"w"] <-  data$y[i]*data$y[i-1]*data$y[i-2]
  }
  data[,"w"] <- data[,"w"]^(1/3)
  if(anual){data[,"w"] <- data[,"w"]^12}
  data$z <- (data$w - 1)*100
  st <- ts(data$z, start = start(x), end = end(x), freq = 12)
  st
}