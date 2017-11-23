estab_dcs <- function(y, parametros, h = 24){
  
  params <- c(0.100000000, 0.100000000,-1.315336496, 4.268469972, 0.549590894,-0.001452403, 0.448915171,-0.092701743, 0.141294841, 0.101951308,-0.048848711,-0.265443596,-0.086032216,-0.142166490,-0.212792790,-0.067664267, 0.064508602)
  
  parametros <- list(
    par = data.frame(
      name = c("k1","ks","f2","df","mu[0]","beta", paste0("gamma",1:11)),
      value = params,
      lower = c(0,0,-Inf,4,-Inf,-Inf, rep(-Inf,11)),
      upper = c(0.1,Inf,Inf,Inf,Inf,Inf,rep(Inf,11))
    ),
    gamma = NA
  )
  parametros
  dcs_res <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2_beta")
}