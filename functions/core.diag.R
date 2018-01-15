core.diag <- function(y, core, conf = 0.95, test = "bias", lags_y = 1:12, lags_core = 1:12, lag_level = c(12,12,12,12,12,12), lag_diff = c(12,12,12,12,12,12)){
  
  # conf: confidence level (0 < conf < 1)
  # lags: ur.df lag (numeric -> c(y, core, diff_y, diff_core))
  # type: ur.df type (trend, drift, none -> c(y, core, diff_y, diff_core))
  
  data <- na.omit(cbind(y,core))
  y <- data[,1]
  core <- data[,2]
  
  # cortar na metade
  n <- length(y)
  y0 <- ts(y[1:round((n/2))], start = start(y), freq = 12)
  y1 <- ts(y[round((n/2)+1):n], end = end(y), freq = 12)
  core0 <- ts(core[1:round((n/2))], start = start(core), freq = 12)
  core1 <- ts(core[round((n/2)+1):n], end = end(core), freq = 12)
  
  if(test == "bias"){
    
    # diag 1: estatísticas descritivas --------
    
    # full
    tabela_desc <- data.frame(matrix(NA, nrow = 2, ncol = 5))
    colnames(tabela_desc) <- c("Mean","Median","Std. Error","CV", "Bias")
    rownames(tabela_desc) <- c("y","core")
    tabela_desc[,1] <- c(mean(y), mean(core))
    tabela_desc[,2] <- c(median(y), median(core))
    tabela_desc[,3] <- c(sd(y), sd(core))
    tabela_desc[,4] <- tabela_desc[,3]/tabela_desc[,1]
    tabela_desc[,5] <- c(NA, tabela_desc[-1,1]-tabela_desc[1,1])
    tabela_desc <- round(tabela_desc,4)
    
    # part1
    tabela_desc0 <- data.frame(matrix(NA, nrow = 2, ncol = 5))
    colnames(tabela_desc0) <- c("Mean","Median","Std. Error","CV", "Bias")
    rownames(tabela_desc0) <- c("y","core")
    tabela_desc0[,1] <- c(mean(y0), mean(core0))
    tabela_desc0[,2] <- c(median(y0), median(core0))
    tabela_desc0[,3] <- c(sd(y0), sd(core0))
    tabela_desc0[,4] <- tabela_desc[,3]/tabela_desc[,1]
    tabela_desc0[,5] <- c(NA, tabela_desc0[-1,1]-tabela_desc0[1,1])
    tabela_desc0 <- round(tabela_desc0,4)
    
    # part2
    tabela_desc1 <- data.frame(matrix(NA, nrow = 2, ncol = 5))
    colnames(tabela_desc1) <- c("Mean","Median","Std. Error","CV", "Bias")
    rownames(tabela_desc1) <- c("y","core")
    tabela_desc1[,1] <- c(mean(y1), mean(core1))
    tabela_desc1[,2] <- c(median(y1), median(core1))
    tabela_desc1[,3] <- c(sd(y1), sd(core1))
    tabela_desc1[,4] <- tabela_desc1[,3]/tabela_desc1[,1]
    tabela_desc1[,5] <- c(NA, tabela_desc1[-1,1]-tabela_desc1[1,1])
    tabela_desc1 <- round(tabela_desc1,4)
    
    # diag 2: teste do viés --------
    
    # modelo irrestrito - full
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
    
    ############ modelo irrestrito - part1
    reg_ur <- lm(y0 ~ core0)
    summary_reg <- summary(reg_ur)
    coef <- summary_reg$coefficients[,"Estimate"]
    sds <- summary_reg$coefficients[,"Std. Error"]
    t_alpha <- coef[1]/sds[1]
    t_beta <- (coef[2] - 1)/sds[2]
    t_test <- qt((1+conf)/2,length(y0) - 2)
    
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
    SSR_r <- sum((y0 - core0)^2)
    SSR_ur <- sum(resid(reg_ur)^2)
    F0 <- ((SSR_r - SSR_ur)/2)/(SSR_ur/(length(y0) - 2))
    Fc <- qf(0.95, 2,length(y0) - 2)
    if(F0 > Fc){ 
      msg3 <- paste0("   Conclusion: reject H0 (",conf*100,"% confidence level)")
    }else{
      msg3 <- paste0("   Conclusion: don't reject H0 (",conf*100,"% confidence level)")
    }
    pvalor <- 1 - pf(F0, 2, length(y0) - 2)
    
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
    tabela_desc0$conclusion <- c(NA, ifelse(diag2$conclusion[3] == "don't reject H0", "unbiased", "biased"))
    
    ############ modelo irrestrito - part2
    reg_ur <- lm(y1 ~ core1)
    summary_reg <- summary(reg_ur)
    coef <- summary_reg$coefficients[,"Estimate"]
    sds <- summary_reg$coefficients[,"Std. Error"]
    t_alpha <- coef[1]/sds[1]
    t_beta <- (coef[2] - 1)/sds[2]
    t_test <- qt((1+conf)/2,length(y1) - 2)
    
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
    SSR_r <- sum((y1 - core1)^2)
    SSR_ur <- sum(resid(reg_ur)^2)
    F0 <- ((SSR_r - SSR_ur)/2)/(SSR_ur/(length(y1) - 2))
    Fc <- qf(0.95, 2,length(y1) - 2)
    if(F0 > Fc){ 
      msg3 <- paste0("   Conclusion: reject H0 (",conf*100,"% confidence level)")
    }else{
      msg3 <- paste0("   Conclusion: don't reject H0 (",conf*100,"% confidence level)")
    }
    pvalor <- 1 - pf(F0, 2, length(y1) - 2)
    
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
    tabela_desc1$conclusion <- c(NA, ifelse(diag2$conclusion[3] == "don't reject H0", "unbiased", "biased"))
    
    # output
    list(out = list(full = tabela_desc,
                    part1 = tabela_desc0,
                    part2 = tabela_desc1)) 
    
  }else if(test == "root"){
    
    # diag 3: estacionariedade ----
    
    
    
    # y completo e partes
    adf10 <- urca::ur.df(y, lags = lag_level[1], type = "drift")
    adf11 <- urca::ur.df(y0, lags = lag_level[2], type = "drift")
    adf12 <- urca::ur.df(y1, lags = lag_level[3], type = "drift")
    
    msg_urdf10 <- ifelse(adf10@teststat[1] < -2.88, "stationary", "not stationary")
    msg_urdf11 <- ifelse(adf11@teststat[1] < -2.88, "stationary", "not stationary")
    msg_urdf12 <- ifelse(adf12@teststat[1] < -2.88, "stationary", "not stationary")
    
    # diff y completo e partes
    adf100 <- urca::ur.df(diff(y), lags = lag_diff[1], type = "none")
    adf111 <- urca::ur.df(diff(y0), lags = lag_diff[2], type = "none")
    adf122 <- urca::ur.df(diff(y1), lags = lag_diff[3], type = "none")
    
    msg_urdf100 <- ifelse(adf100@teststat[1] < -1.95, "stationary", "not stationary")
    msg_urdf111 <- ifelse(adf111@teststat[1] < -1.95, "stationary", "not stationary")
    msg_urdf122 <- ifelse(adf122@teststat[1] < -1.95, "stationary", "not stationary")
    
    # core completo e partes
    adf20 <- urca::ur.df(core, lags = lag_level[4], type = "drift")
    adf21 <- urca::ur.df(core0, lags = lag_level[5], type = "drift")
    adf22 <- urca::ur.df(core1, lags = lag_level[6], type = "drift")
    
    msg_urdf20 <- ifelse(adf20@teststat[1] < -2.88, "stationary", "not stationary")
    msg_urdf21 <- ifelse(adf21@teststat[1] < -2.88, "stationary", "not stationary")
    msg_urdf22 <- ifelse(adf22@teststat[1] < -2.88, "stationary", "not stationary")
    
    # diff y completo e partes
    adf200 <- urca::ur.df(diff(core), lags = lag_diff[4], type = "none")
    adf211 <- urca::ur.df(diff(core0), lags = lag_diff[5], type = "none")
    adf222 <- urca::ur.df(diff(core1), lags = lag_diff[6], type = "none")
    
    msg_urdf200 <- ifelse(adf200@teststat[1] < -1.95, "stationary", "not stationary")
    msg_urdf211 <- ifelse(adf211@teststat[1] < -1.95, "stationary", "not stationary")
    msg_urdf222 <- ifelse(adf222@teststat[1] < -1.95, "stationary", "not stationary")
    
    
    out_urdf <- data.frame(y_level = c(msg_urdf10,msg_urdf11,msg_urdf12),
                           y_diff = c(msg_urdf100,msg_urdf111,msg_urdf122),
                           core_level = c(msg_urdf20,msg_urdf21,msg_urdf22),
                           core_diff = c(msg_urdf200,msg_urdf211,msg_urdf222))
    print(out_urdf)
    out_urtest <- list(
      full = list(
        y = list(level = adf10,
                 diff = adf100),
        core = list(level = adf20,
                    diff = adf200)
      ),
      part1 = list(
        y = list(level = adf11,
                 diff = adf111),
        core = list(level = adf21,
                    diff = adf211)
      ),
      part2 = list(
        y = list(level = adf12,
                 diff = adf122),
        core = list(level = adf22,
                    diff = adf222)
      )
    )
    
    # output
    invisible(list(out = out_urdf, adf = out_urtest) )
    
  }else if(test == "coint"){
    
    # diag 4: cointegracao ----
    coint <- summary(ca.jo(cbind(y,core), type = "eigen"))
    out_coint <- ifelse(coint@teststat[2] > coint@cval[2,2], "cointegrated time series", "non-cointegrated time series")
    print(out_coint)
    # output
    invisible(list(out = coint))
   
  }else if(test == "attract"){
    
    # diag 5: atratividade ----
    erro <- ts(resid(lm(y ~ core)), end = end(core), freq = 12)
    #print(summary(lm(y ~ core)))
    # dum <- cbind(BETS::BETS.dummy(start = start(data), end = end(data), month = 11, year = 2002),
    #              BETS::BETS.dummy(start = start(data), end = end(data), month = 7, year = 2000),
    #              BETS::BETS.dummy(start = start(data), end = end(data), month = 12, year = 2002))
    # 
    # colnames(dum) <- c("dummy1","dummy2","dummy3")
    difs <- na.omit(cbind(diff(y), diff(core), lag(erro, k = -1)))
    colnames(difs) <- c("y","core","erro_1")
     
    modelo_y <- dynlm(difs[,"y"] ~ difs[,"erro_1"] + 
                        L(difs[,"y"],lags_y))
    
    print(summary(modelo_y))
    
    par(mfrow = c(2,2))
    plot((modelo_y$residuals - mean(modelo_y$residuals))/sd(modelo_y$residuals))
    abline(h = c(-3,3), lty = 2, col = 2)
    acf(modelo_y$residuals, lag.max = 36)
     
    modelo_core <- dynlm(difs[,"core"] ~ difs[,"erro_1"] +
                           L(difs[,"core"],lags_core))
    
    print(summary(modelo_core))
    plot((modelo_core$residuals - mean(modelo_core$residuals))/sd(modelo_core$residuals))
    abline(h = c(-3,3), lty = 2, col = 2)
    acf(modelo_core$residuals, lag.max = 36)   
    
    if(summary(modelo_y)$coefficients[2,1] < 0 & summary(modelo_y)$coefficients[2,4] < 0.05
       & summary(modelo_core)$coefficients[2,4] > 0.05){
      message("core attracts inflation")
    }else if(summary(modelo_y)$coefficients[2,1] < 0 & summary(modelo_y)$coefficients[2,4] < 0.05
            & summary(modelo_core)$coefficients[2,1] > 0 & summary(modelo_core)$coefficients[2,4] < 0.05){
      message("core attracts inflation and inflation attracts core")
    }
    
    # output
    invisible(list(out_y = modelo_y, out_core = modelo_core))
  
  }else if(test == "fcst"){
  

    # diag 6: previsão ----------------------

    pi3 <- lag(y, k = -3)
    pic3 <- lag(core, k = -3)
    d <- na.omit(cbind(y - pi3, pic3 - pi3))
    colnames(d) <- c("y","x")
    
    md <- lm(d[,1] ~ d[,2])
    message("h = 3 meses")
    print(summary(md))
    
    pi6 <- lag(y, k = -6)
    pic6 <- lag(core, k = -6)
    d <- na.omit(cbind(y - pi6, pic6 - pi6))
    colnames(d) <- c("y","x")
    
    md <- lm(d[,1] ~ d[,2])
    message("h = 6 meses")
    print(summary(md))
    
    
    pi12 <- lag(y, k = -12)
    pic12 <- lag(core, k = -12)
    d <- na.omit(cbind(y - pi12, pic12 - pi12))
    colnames(d) <- c("y","x")

    md <- lm(d[,1] ~ d[,2])
    message("h = 12 meses")
    print(summary(md))

    pi24 <- lag(y, k = -24)
    pic24 <- lag(core, k = -24)
    d <- na.omit(cbind(y - pi24, pic24 - pi24))
    colnames(d) <- c("y","x")

    md <- lm(d[,1] ~ d[,2])
    message("h = 24 meses")
    print(summary(md))

  }else if(test == "seas"){
    
    # diag 7: sazonalidade
    out <- rbind(qs(seas(core))[1,],
                qs(seas(core0))[1,],
                qs(seas(core1))[1,])
    rownames(out) <- c("full","part1","part2")
    
    list(out = out)
  
  }
   
  
  
  
}

