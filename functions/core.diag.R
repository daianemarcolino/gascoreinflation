core.diag <- function(y, core, par, conf = 0.95){
  
  # conf: confidence level (0 < conf < 1)
  # lags: ur.df lag (numeric -> c(y, core, diff_y, diff_core))
  # type: ur.df type (trend, drift, none -> c(y, core, diff_y, diff_core))
  
  data <- na.omit(cbind(y,core))
  y <- data[,1]
  core <- data[,2]
  
  # diag 1: estatísticas descritivas --------
  tabela_desc <- data.frame(matrix(NA, nrow = 2, ncol = 5))
  colnames(tabela_desc) <- c("Mean","Median","Std. Error","CV", "Bias")
  rownames(tabela_desc) <- c("y","core")
  tabela_desc[,1] <- c(mean(y), mean(core))
  tabela_desc[,2] <- c(median(y), median(core))
  tabela_desc[,3] <- c(sd(y), sd(core))
  tabela_desc[,4] <- tabela_desc[,3]/tabela_desc[,1]
  tabela_desc[,5] <- c(NA, tabela_desc[-1,1]-tabela_desc[1,1])
  tabela_desc <- round(tabela_desc,4)
  
  # diag 2: teste do viés --------
  
  # modelo irrestrito
  reg_ur <- lm(y ~ core)
  summary_reg <- summary(reg_ur)
  coef <- summary_reg$coefficients[,"Estimate"]
  sds <- summary_reg$coefficients[,"Std. Error"]
  t_alpha <- coef[1]/sds[1]
  t_beta <- (coef[2] - 1)/sds[2]
  t_test <- qt((1+conf)/2,length(y) - 2)
  
  # teste para alpha
  if(abs(t_alpha) > t_test){
    msg1 <- paste0("   Conclusion: reject H0 (",conf*100,"% confidence level)")
  }else{
    msg1 <- paste0("   Conclusion: don't reject H0 (",conf*100,"% confidence level)")
  }
  
  # teste para beta
  if(abs(t_beta) > t_test){
    msg2 <- paste0("   Conclusion: reject H0 (",conf*100,"% confidence level)")
  }else{
    msg2 <- paste0("   Conclusion: don't reject H0 (",conf*100,"% confidence level)")
  }
  
  # Soma de quadrado dos residuos
  SSR_r <- sum((y - core)^2)
  SSR_ur <- sum(resid(reg_ur)^2)
  F0 <- ((SSR_r - SSR_ur)/2)/(SSR_ur/(length(y) - 2))
  Fc <- qf(0.95, 2,length(y) - 2)
  if(F0 > Fc){ 
    msg3 <- paste0("   Conclusion: reject H0 (",conf*100,"% confidence level)")
  }else{
    msg3 <- paste0("   Conclusion: don't reject H0 (",conf*100,"% confidence level)")
  }
  pvalor <- 1 - pf(F0, 2, length(y) - 2)

  message("- ALPHA TEST -> H0: alpha = 0 \n", 
      paste("   coef:", round(coef[1],4)),
      paste("\n   sd:", round(sds[1],4)), "\n", msg1, " ","\n \n",
      "- BETA TEST -> H0: beta = 1 \n", 
      paste("   coef:", round(coef[2],4)),
      paste("\n   sd:", round(sds[2],4)),  "\n", msg2, " ", "\n \n",
      "- ALPHA E BETA TEST -> H0: alpha = 0 e beta = 1","\n", 
      "   p-value: ", round(pvalor,4), "\n", msg3)
  names(coef)[2] <- deparse(substitute(x))
  
  diag2 <- data.frame(stats = round(c(t_alpha, t_beta, F0),4),
             t_value = round(c(rep(t_test,2), NA),4),
             f_value = round(c(NA,NA,Fc),4))
  diag2$conclusion <- ifelse(abs(diag2$stats) > diag2$t_value, "reject H0", "don't reject H0")
  diag2$conclusion[3] <- ifelse(diag2$stats[3] > diag2$f_value[3], "reject H0", "don't reject H0")
  rownames(diag2) <- c("H0: alpha = 0", "H0: beta = 1", "H0: alpha = 0 e beta = 1")
  
  # acrescentar vies
  tabela_desc$conclusion <- c(NA, ifelse(diag2$conclusion[3] == "don't reject H0", "unbiased", "biased"))
  
  # diag 3: estacionariedade ----
  
  y0 <- window(y, start = c(2008,4), freq = 12)
  y1 <- window(y, end = c(2008,5), freq = 12)
  
  # y completo e partes
  adf10 <- urca::ur.df(y, lags = 12, type = "drift")
  adf11 <- urca::ur.df(y0, lags = 12, type = "drift")
  adf12 <- urca::ur.df(y1, lags = 12, type = "drift")

  msg_urdf10 <- ifelse(adf10@teststat[1] < -2.88, "stationary", "not stationary")
  msg_urdf11 <- ifelse(adf11@teststat[1] < -2.88, "stationary", "not stationary")
  msg_urdf12 <- ifelse(adf12@teststat[1] < -2.88, "stationary", "not stationary")

  # diff y completo e partes
  adf100 <- urca::ur.df(diff(y), lags = 11, type = "none")
  adf111 <- urca::ur.df(diff(y0), lags = 8, type = "none")
  adf122 <- urca::ur.df(diff(y1), lags = 2, type = "none")

  msg_urdf100 <- ifelse(adf100@teststat[1] < -1.95, "stationary", "not stationary")
  msg_urdf111 <- ifelse(adf111@teststat[1] < -1.95, "stationary", "not stationary")
  msg_urdf122 <- ifelse(adf122@teststat[1] < -1.95, "stationary", "not stationary")
  
  core0 <- window(core, start = c(2008,4), freq = 12)
  core1 <- window(core, end = c(2008,5), freq = 12)
  
  # core completo e partes
  adf20 <- urca::ur.df(core, lags = 12, type = "drift")
  adf21 <- urca::ur.df(core0, lags = 10, type = "drift")
  adf22 <- urca::ur.df(core1, lags = 12, type = "drift")

  msg_urdf20 <- ifelse(adf20@teststat[1] < -2.88, "stationary", "not stationary")
  msg_urdf21 <- ifelse(adf21@teststat[1] < -2.88, "stationary", "not stationary")
  msg_urdf22 <- ifelse(adf22@teststat[1] < -2.88, "stationary", "not stationary")
  
  # diff y completo e partes
  adf200 <- urca::ur.df(diff(core), lags = 11, type = "none")
  adf211 <- urca::ur.df(diff(core0), lags = 2, type = "none")
  adf222 <- urca::ur.df(diff(core1), lags = 11, type = "none")

  msg_urdf200 <- ifelse(adf200@teststat[1] < -1.95, "stationary", "not stationary")
  msg_urdf211 <- ifelse(adf211@teststat[1] < -1.95, "stationary", "not stationary")
  msg_urdf222 <- ifelse(adf222@teststat[1] < -1.95, "stationary", "not stationary")
  
  
  out_urdf <- data.frame(y_level = c(msg_urdf10,msg_urdf11,msg_urdf12),
                         y_diff = c(msg_urdf100,msg_urdf111,msg_urdf122),
                         core_level = c(msg_urdf20,msg_urdf21,msg_urdf22),
                         core_diff = c(msg_urdf200,msg_urdf211,msg_urdf222))
  
  # diag 4: cointegracao ----
  coint <- summary(ca.jo(cbind(y,core), type = "eigen"))
  out_coint <- ifelse(coint@teststat[2] > coint@cval[2,2], "coint", "not coint")
  
  # diag 5: atratividade ----
  #data <- window(data, start = c(2007,8), freq = 12)
  reg <- lm(data[,1] ~ data[,2])
  erro <- ts(resid(reg), end = end(data), freq = 12)
  
  dum <- BETS::BETS.dummy(start = start(data), end = end(data), month = 12, year = 2002)
  
  difs <- na.omit(cbind(diff(y), diff(core), lag(erro, k = -1), dum))
  colnames(difs) <- c("y","core","erro_1", "dummy")
  
  modelo_y <- dynlm(difs[,"y"] ~ difs[,"erro_1"] + L(difs[,"y"],1) +
                         L(difs[,"y"],2) + 
                         L(difs[,"y"],3) +
                         L(difs[,"y"],4) + 
                         L(difs[,"y"],5) +
                         L(difs[,"y"],6) +
                         #L(difs[,"y"],7) +
                         #L(difs[,"y"],8) + 
                         #L(difs[,"y"],9) + 
                         #L(difs[,"y"],10) +
                         L(difs[,"y"],11) + 
                         L(difs[,"y"],12))
  summary(modelo_y)
  plot(modelo_y$residuals)
  acf(modelo_y$residuals, lag.max = 36)
  
  modelo_core <- dynlm(difs[,"core"] ~ difs[,"erro_1"] + difs[,"dummy"] +
                         L(difs[,"core"],1) +
                         L(difs[,"core"],2) + 
                         L(difs[,"core"],3) +
                          #L(difs[,"core"],4) + 
                          #L(difs[,"core"],5) +
                          L(difs[,"core"],6) +
                          #L(difs[,"core"],7) +
                          L(difs[,"core"],8) + 
                          L(difs[,"core"],9) + 
                          #L(difs[,"core"],10) +
                          L(difs[,"core"],11) + 
                          L(difs[,"core"],13))
  summary(modelo_core)
  plot(modelo_core$residuals)
  acf(modelo_core$residuals, lag.max = 36)
  # lambda_h < 0 e lambda_c > 0 -> as séries se ajustam uma a outra       
  
  # previsão
  
  pi12 <- lag(y, k = -12)
  pic12 <- lag(core, k = -12)
  d <- na.omit(cbind(y - pi12, pic12 - pi12))
  d <- window(na.omit(cbind(y - pi12, pic12 - pi12)), start = c(2006,4), freq = 12)
  colnames(d) <- c("y","x")
  
  md <- lm(d[,1] ~ d[,2])
  summary(md)
  
  pi24 <- lag(y, k = -24)
  pic24 <- lag(core, k = -24)
  d <- na.omit(cbind(y - pi24, pic24 - pi24))
  d <- window(na.omit(cbind(y - pi24, pic24 - pi24)), start = c(2006,4), freq = 12)
  colnames(d) <- c("y","x")
  
  md <- lm(d[,1] ~ d[,2])
  summary(md)
  
  # output
  list(stats = tabela_desc, bias.test = diag2, 
       urdf = list(y = list(level = adf1, diff = adf3), 
                   core = list(level = adf2, diff = adf4)
       ),
       coint = coint
  ) 
  
  
}

