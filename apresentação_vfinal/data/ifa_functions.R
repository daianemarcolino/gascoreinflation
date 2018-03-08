# funções utilizadas no arquivo 02 Rotina.R

# leitura --------------------------
# arrumar arquivo do IPC
ipcts <- function(x, start = c(1999,1)){
  # x: data.frame
  
  x[,2] <- as.character(x[,2])
  codigos <- data.frame(cod = x$Codigo, descricao = x$Descricao)
  codigos$cod <- paste0("cod_",codigos$cod)
  
  # remover colunas codigo e descricao da tabela
  x2 <- x[,-c(1,2)]
  rownames(x2) <- codigos$cod
  x3 <- t(x2)
  
  # arrumar datas
  ipc.ts <- ts(x3, start = start, freq = 12)
  
  # Retornar resultados
  list(ipc = ipc.ts[,"cod_0"], subitens = ipc.ts[,-1], codigos = codigos)
}

# formatar taxas -------------------------------------------------------
# Acumulado nos últimos 12 meses
acum12 <- function(x){ 
  # x: ts 
  n <- length(x)
  data <- data.frame(x, y = rep(NA, n), w = rep(NA, n), z = rep(NA, n), row.names = 1:n)
  data$y <- data$x/100 + 1
  for(i in n:13){
    data[i,"w"] <-  data$y[i]*data$y[i-1]*data$y[i-2]*
      data$y[i-3]*data$y[i-4]*data$y[i-5]*data$y[i-6]*
      data$y[i-7]*data$y[i-8]*data$y[i-9]*data$y[i-10]*data$y[i-11]
  }
  data$z <- (data$w - 1)*100
  st <- ts(data$z, start = start(x), end = end(x), freq = 12)
  st  
}

acumano <- function(x){ 
  inicio <- start(x)
  fim <- end(x)
  inicio_ok <- c(inicio[1],1)
  fim_ok <- c(fim[1],12)
  x_completo <- ts(NA, start = inicio_ok, end = fim_ok, freq = 12)
  dados <- data.frame(cbind(x,x_completo))
  matriz <- matrix(dados[,1], ncol = 12, byrow = T)#, nrow = 18)
  matriz <- (matriz/100 + 1)
  prod <- (apply(matriz, MARGIN = 1, FUN = cumprod) - 1)*100
  acum <- ts(matrix(prod, ncol = 1), start = start(x_completo), freq = 12)
  acum
}

# Anualizado
anual <- function(x){
  # x: ts
  ((x/100+1)^12-1)*100
}

# Relativo
relativo <- function(x){
  x/100+1
}

# acumulada em três meses + anualizar
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

# núcleos tradicionais --------------

# verificar volatilidade
vola1999 <- function(x, alpha = 2){
  # x: mts
  x1 <- window(x, start = c(1999,1), end = c(2003,12), freq = 12)
  x2 <- window(x, start = c(2004,1), end = c(2009,12), freq = 12)
  x3 <- window(x, start = c(2010,1), end = end(x), freq = 12)
  x4 <- window(x, start = c(2004,1), end = end(x), freq = 12)
  
  desvios <- apply(x, MARGIN = 2, FUN = sd, na.rm = T)
  desvios1 <- apply(x1, MARGIN = 2, FUN = sd, na.rm = T)
  desvios2 <- apply(x2, MARGIN = 2, FUN = sd, na.rm = T)
  desvios3 <- apply(x3, MARGIN = 2, FUN = sd, na.rm = T)
  #desvios4 <- apply(x4, MARGIN = 2, FUN = sd, na.rm = T)
  
  desvios <- (desvios - mean(desvios, na.rm = T))/sd(desvios, na.rm = T)
  desvios1 <- (desvios1 - mean(desvios1, na.rm = T))/sd(desvios1, na.rm = T)
  desvios2 <- (desvios2 - mean(desvios2, na.rm = T))/sd(desvios2, na.rm = T)
  desvios3 <- (desvios3 - mean(desvios3, na.rm = T))/sd(desvios3, na.rm = T)
  #desvios4 <- (desvios4 - mean(desvios4, na.rm = T))/sd(desvios4, na.rm = T)
  
  vola <- data.frame(matrix(NA, ncol = 5, nrow = ncol(x)))
  rownames(vola) <- colnames(x)
  colnames(vola) <- c("Descricao","sd_total","sd_1999_2003","sd_2004_2009","sd_2010_fim")#,"sd_2004_fim")
  vola$Descricao <- codigos$descricao[-1]
  vola$sd_total <- desvios
  vola$sd_1999_2003 <- desvios1
  vola$sd_2004_2009 <- desvios2
  vola$sd_2010_fim <- desvios3
  # vola$sd_2004_fim <- desvios4 
  
  masc1 <- abs(vola$sd_total) > alpha
  masc2 <- abs(vola$sd_1999_2003) > alpha
  masc3 <- abs(vola$sd_2004_2009) > alpha
  masc4 <- abs(vola$sd_2010_fim) > alpha
  #masc5 <- abs(vola$sd_2004_fim) > alpha
  
  vola$SAI1 <- masc1
  vola$SAI2 <- rowSums(cbind(masc2,masc3,masc4), na.rm = T) > 1
  # vola$SAI3 <- masc5
  
  vola <- vola[order(vola$sd_total, decreasing = T),]
  vola
}

# verificar volatilidade
vola <- function(x, order = T, alpha = 2){
  # x: mts
  x1 <- window(x, start = c(2005,1), end = c(2008,12), freq = 12)
  x2 <- window(x, start = c(2009,1), end = end(x), freq = 12)
  
  desvios <- apply(x, MARGIN = 2, FUN = sd, na.rm = T)
  desvios1 <- apply(x1, MARGIN = 2, FUN = sd, na.rm = T)
  desvios2 <- apply(x2, MARGIN = 2, FUN = sd, na.rm = T)
  
  desvios <- (desvios - mean(desvios, na.rm = T))/sd(desvios, na.rm = T)
  desvios1 <- (desvios1 - mean(desvios1, na.rm = T))/sd(desvios1, na.rm = T)
  desvios2 <- (desvios2 - mean(desvios2, na.rm = T))/sd(desvios2, na.rm = T)
  
  vola <- data.frame(matrix(NA, ncol = 4, nrow = ncol(x)))
  rownames(vola) <- colnames(x)
  colnames(vola) <- c("Descricao","sd_total","sd_2005_2008","sd_2009_fim")
  vola$Descricao <- codigos$descricao[-1]
  vola$sd_total <- desvios
  vola$sd_2005_2008 <- desvios1
  vola$sd_2009_fim <- desvios2
  
  masc1 <- abs(vola$sd_total) > alpha 
  masc2 <- abs(vola$sd_2005_2008) > alpha 
  masc3 <- abs(vola$sd_2009_fim) > alpha 
  
  vola$SAI1 <- masc1
  vola$SAI2 <- sum(masc2 + masc3, na.rm = T) > 0
  
  if(order){ vola <- vola[order(vola$sd_total, decreasing = T),] }
  vola
}

core.ex <- function(variacao, pesos, sd){
  # sd: saida de vola1999 ou vola
  saem <- rownames(sd[sd$SAI1 | sd$SAI2, ])
  pesos[,colnames(pesos) %in% saem] <- 0
  pesos <- (pesos/rowSums(pesos, na.rm = T))*100
  
  ts(rowSums(variacao*pesos, na.rm = T)/100, start = start(variacao), freq = 12)
}

core.dp <- function(inf, sub, pesos, j = 48){
  # inf: índice da inflação (variação)
  # sub: variação dos subitens
  # pesos: pesos dos subitens
  # j: janela para o cálculo da volatilidade
  
  pi_rel <- sub - inf
  colnames(pi_rel) <- colnames(sub)
  
  data <- as.Date(paste0(start(pi_rel)[1], "-", start(pi_rel)[2],"-01"))
  data <- seq(data, by = "1 month", length.out = nrow(pi_rel))
  novos_pesos <- sub*0
  
  for(t in nrow(pi_rel):j){
    
    ano1 <- as.numeric(substr(data[t-j+1],1,4))
    mes1 <- as.numeric(substr(data[t-j+1],6,7))
    
    ano2 <- as.numeric(substr(data[t],1,4))
    mes2 <- as.numeric(substr(data[t],6,7))
    
    sub_pi_rel <- window(pi_rel, start = c(ano1,mes1), end = c(ano2,mes2), freq = 12)
    desvios <- apply(sub_pi_rel, MARGIN = 2, FUN = sd, na.rm = T)
    novos_pesos[t,] <- ((1/desvios)/sum(1/desvios, na.rm = T))*pesos[t,]
  }
  
  novos_pesos2 <- window(novos_pesos, start = c(as.numeric(substr(data[j+1],1,4)), as.numeric(substr(data[j+1],6,7))),
                         freq = 12)
  novos_pesos3 <- novos_pesos2/rowSums(novos_pesos2,na.rm = T)
  sub2 <- window(sub, start = start(novos_pesos2), freq = 12)
  core <- ts(rowSums(novos_pesos3*sub2, na.rm = T), start = start(sub2), freq = 12)
  core
}


core.ma <- function(sub, pesos, inf = 20, sup = 20, suave = F){
  # sub: subitens (variação)
  # pesos: pesos dos subitens
  # sup e inf: corte percentual das caudas inferior e superior

  pesos_novos <- pesos*0
  #pesos_novos[is.na(pesos_novos)] <- 0
  #pesos[is.na(pesos)] <- 0
  
  if(suave){
    # itens que serão suavizados
    
    codigos_suave <- c(210101,210103,220101,220103,220105,220111,410301,410305,410307,410309,
                       410319,420105,420107,420111,420114,420115,420116,420117,420118,420123,
                       420125,420126,420131,420133,420135,420137,510101,510103,510105,510107,
                       510153,510155,610101,610103,610105,610107,610109,610111,610115,610303,
                       620501,620503,620505,620507,620509,620903,620909,620913,720101,720105,
                       720115,720301,720305,810101,810105)
   
    codigos_suave <- paste0("cod_",codigos_suave)
    filtro_suave <- sub[,colnames(sub) %in% codigos_suave]
    filtro_suave <- (filtro_suave/100 + 1)^(1/12)
    filtro_suave2 <- filtro_suave
    for(i in 12:nrow(filtro_suave)){
      v <- apply(filtro_suave[(i-11):i,], MARGIN = 2, FUN = prod)
      masc_v <- !is.na(v)
      filtro_suave2[i,masc_v] <- v[masc_v]
    }
    filtro_suave2 <- (filtro_suave2-1)*100
    sub[,colnames(sub) %in% codigos_suave] <- filtro_suave2
  }
  
  for(i in 1:nrow(sub)){
    dados <- data.frame(x = sub[i,], pesos = pesos[i,])
    dados <- dados[order(dados$x),]
    dados$pesos_acum <- cumsum(dados$pesos)
    pos_inf <- min(which(dados$pesos_acum > inf))
    pos_sup <- min(which(dados$pesos_acum > 100-sup))
    
    dados_aux1 <- dados[1:(pos_inf-1),]
    dados_aux3 <- dados[(pos_inf+1):(pos_sup-1),]
    dados_aux5 <- dados[(pos_sup+1):nrow(dados),]
    
    peso_sai_inf <- inf - dados[pos_inf-1, "pesos_acum"]
    peso_entra_inf <- dados[pos_inf, "pesos_acum"] - inf 
    peso_entra_sup <- (100-sup) - dados[pos_sup-1, "pesos_acum"]
    peso_sai_sup <- dados[pos_sup, "pesos_acum"] - (100-sup) 
    
    dados_aux2 <- rbind(dados[pos_inf,],dados[pos_inf,])
    dados_aux4 <- rbind(dados[pos_inf,],dados[pos_inf,])
    
    dados_aux2$pesos <- c(peso_sai_inf, peso_entra_inf)
    dados_aux4$pesos <- c(peso_entra_sup, peso_sai_sup)
    
    rownames(dados_aux2) <-  c(paste0(rownames(dados[pos_inf,]),"_sai"),
                               rownames(dados[pos_inf,]))
    
    rownames(dados_aux4) <-  c(rownames(dados[pos_sup,]),
                               paste0(rownames(dados[pos_sup,]),"_sai"))
    
    dados_aux <- rbind(dados_aux1, dados_aux2,
                       dados_aux3, dados_aux4, dados_aux5)
    dados_aux$pesos_acum <- cumsum(dados_aux$pesos)
    dados_aux$conta <- 0 
    dados_aux[dados_aux$pesos_acum > inf & dados_aux$pesos_acum < 100 - sup | dados_aux$pesos_acum == as.character(100 - sup), "conta"] <- 1
    
    dados_aux[grepl("sai", rownames(dados_aux)), "conta"] <- 0
    
    for(j in rownames(subset(dados_aux, dados_aux$conta == 1))){
      pesos_novos[i,j] <- dados_aux[j,"pesos"]
    }
    
  }
  pesos_novos2 <- pesos_novos/rowSums(pesos_novos,na.rm = T)*100
  
  core <- ts(rowSums(pesos_novos2*sub, na.rm = T)/100, start = start(sub), freq = 12)
  core
}

# núcleos ma para serviços --------------

core.ma.sv <- function(sub, pesos, inf = 20, sup = 20, suave = F){
  # sub: subitens (variação)
  # pesos: pesos dos subitens
  # sup e inf: corte percentual das caudas inferior e superior
  
  pesos_novos <- pesos*0
  #pesos_novos[is.na(pesos_novos)] <- 0
  #pesos[is.na(pesos)] <- 0
  
  if(suave){
    # itens que serão suavizados
    
    codigos_suave <- c(510101,510103,510105,510107,510153,510155)
    
    codigos_suave <- paste0("cod_",codigos_suave)
    filtro_suave <- sub[,colnames(sub) %in% codigos_suave]
    filtro_suave <- (filtro_suave/100 + 1)^(1/12)
    filtro_suave2 <- filtro_suave
    for(i in 12:nrow(filtro_suave)){
      v <- apply(filtro_suave[(i-11):i,], MARGIN = 2, FUN = prod)
      masc_v <- !is.na(v)
      filtro_suave2[i,masc_v] <- v[masc_v]
    }
    filtro_suave2 <- (filtro_suave2-1)*100
    sub[,colnames(sub) %in% codigos_suave] <- filtro_suave2
  }
  
  for(i in 1:nrow(sub)){
    dados <- data.frame(x = sub[i,], pesos = pesos[i,])
    dados <- dados[order(dados$x),]
    dados$pesos_acum <- cumsum(dados$pesos)
    pos_inf <- min(which(dados$pesos_acum > inf))
    pos_sup <- min(which(dados$pesos_acum > 100-sup))
    
    dados_aux1 <- dados[1:(pos_inf-1),]
    dados_aux3 <- dados[(pos_inf+1):(pos_sup-1),]
    dados_aux5 <- dados[(pos_sup+1):nrow(dados),]
    
    peso_sai_inf <- inf - dados[pos_inf-1, "pesos_acum"]
    peso_entra_inf <- dados[pos_inf, "pesos_acum"] - inf 
    peso_entra_sup <- (100-sup) - dados[pos_sup-1, "pesos_acum"]
    peso_sai_sup <- dados[pos_sup, "pesos_acum"] - (100-sup) 
    
    dados_aux2 <- rbind(dados[pos_inf,],dados[pos_inf,])
    dados_aux4 <- rbind(dados[pos_inf,],dados[pos_inf,])
    
    dados_aux2$pesos <- c(peso_sai_inf, peso_entra_inf)
    dados_aux4$pesos <- c(peso_entra_sup, peso_sai_sup)
    
    rownames(dados_aux2) <-  c(paste0(rownames(dados[pos_inf,]),"_sai"),
                               rownames(dados[pos_inf,]))
    
    rownames(dados_aux4) <-  c(rownames(dados[pos_sup,]),
                               paste0(rownames(dados[pos_sup,]),"_sai"))
    
    dados_aux <- rbind(dados_aux1, dados_aux2,
                       dados_aux3, dados_aux4, dados_aux5)
    dados_aux$pesos_acum <- cumsum(dados_aux$pesos)
    dados_aux$conta <- 0 
    dados_aux[dados_aux$pesos_acum > inf & dados_aux$pesos_acum < 100 - sup | dados_aux$pesos_acum == as.character(100 - sup), "conta"] <- 1
    
    dados_aux[grepl("sai", rownames(dados_aux)), "conta"] <- 0
    
    for(j in rownames(subset(dados_aux, dados_aux$conta == 1))){
      pesos_novos[i,j] <- dados_aux[j,"pesos"]
    }
    
  }
  pesos_novos2 <- pesos_novos/rowSums(pesos_novos,na.rm = T)*100
  
  core <- ts(rowSums(pesos_novos2*sub, na.rm = T)/100, start = start(sub), freq = 12)
  list(core = core, pesos = pesos_novos2)
}

# filtro ma para serviços bacen ------------------------------------
core.ma.bacen <- function(sub, pesos, inf = 20, sup = 20, suave = F){
  # sub: subitens (variação)
  # pesos: pesos dos subitens
  # sup e inf: corte percentual das caudas inferior e superior
  
  pesos_novos <- pesos*0
  #pesos_novos[is.na(pesos_novos)] <- 0
  #pesos[is.na(pesos)] <- 0
  
  # if(suave){
  #   # itens que serão suavizados
  #   
  #   codigos_suave <- c(210101,210103,220101,220103,220105,220111,410301,410305,410307,410309,
  #                      410319,420105,420107,420111,420114,420115,420116,420117,420118,420123,
  #                      420125,420126,420131,420133,420135,420137,510101,510103,510105,510107,
  #                      510153,510155,610101,610103,610105,610107,610109,610111,610115,610303,
  #                      620501,620503,620505,620507,620509,620903,620909,620913,720101,720105,
  #                      720115,720301,720305,810101,810105)
  #   
  #   codigos_suave <- paste0("cod_",codigos_suave)
  #   filtro_suave <- sub[,colnames(sub) %in% codigos_suave]
  #   filtro_suave <- (filtro_suave/100 + 1)^(1/12)
  #   filtro_suave2 <- filtro_suave
  #   for(i in 12:nrow(filtro_suave)){
  #     v <- apply(filtro_suave[(i-11):i,], MARGIN = 2, FUN = prod)
  #     masc_v <- !is.na(v)
  #     filtro_suave2[i,masc_v] <- v[masc_v]
  #   }
  #   filtro_suave2 <- (filtro_suave2-1)*100
  #   sub[,colnames(sub) %in% codigos_suave] <- filtro_suave2
  # }
  
  for(i in 1:nrow(sub)){
    dados <- data.frame(x = sub[i,], pesos = pesos[i,])
    dados <- dados[order(dados$x),]
    dados$pesos_acum <- cumsum(dados$pesos)
    pos_inf <- min(which(dados$pesos_acum > inf))
    pos_sup <- min(which(dados$pesos_acum > 100-sup))
    
    dados_aux1 <- dados[1:(pos_inf-1),]
    dados_aux3 <- dados[(pos_inf+1):(pos_sup-1),]
    dados_aux5 <- na.omit(dados[(pos_sup+1):nrow(dados),])
    
    if(pos_inf == 1){
      peso_sai_inf <- inf
    }else{
      peso_sai_inf <- inf - dados[pos_inf-1, "pesos_acum"]
    }
    peso_entra_inf <- dados[pos_inf, "pesos_acum"] - inf 
   
    peso_entra_sup <- (100-sup) - dados[pos_sup-1, "pesos_acum"]
    peso_sai_sup <- dados[pos_sup, "pesos_acum"] - (100-sup) 
    
    dados_aux2 <- rbind(dados[pos_inf,],dados[pos_inf,])
    dados_aux4 <- rbind(dados[pos_inf,],dados[pos_inf,])
    
    dados_aux2$pesos <- c(peso_sai_inf, peso_entra_inf)
    dados_aux4$pesos <- c(peso_entra_sup, peso_sai_sup)
    
    rownames(dados_aux2) <-  c(paste0(rownames(dados[pos_inf,]),"_sai"),
                               rownames(dados[pos_inf,]))
    
    rownames(dados_aux4) <-  c(rownames(dados[pos_sup,]),
                               paste0(rownames(dados[pos_sup,]),"_sai"))
    
    if(pos_inf == 1){
      dados_aux <- rbind(dados_aux2,
                         dados_aux3, dados_aux4, dados_aux5)
    }else{
      dados_aux <- rbind(dados_aux1, dados_aux2,
                         dados_aux3, dados_aux4, dados_aux5)
    }
    
    
    dados_aux$pesos_acum <- cumsum(dados_aux$pesos)
    dados_aux$conta <- 0 
    dados_aux[dados_aux$pesos_acum > inf & dados_aux$pesos_acum < 100 - sup | dados_aux$pesos_acum == as.character(100 - sup), "conta"] <- 1
    
    dados_aux[grepl("sai", rownames(dados_aux)), "conta"] <- 0
    
    for(j in rownames(subset(dados_aux, dados_aux$conta == 1))){
      pesos_novos[i,j] <- dados_aux[j,"pesos"]
    }
    
  }
  pesos_novos2 <- pesos_novos/rowSums(pesos_novos,na.rm = T)*100
  
  core <- ts(rowSums(pesos_novos2*sub, na.rm = T)/100, start = start(sub), freq = 12)
  list(core = core, pesos = pesos_novos2)
}

# avaliar medidas ------------

avaliar.vies <- function(y,x, conf = 0.95, cov = NULL){
  if(is.null(cov)){
  reg <- lm(y ~ x)}else{reg <- lm(y ~ x + cov)}
  reg2 <- summary(reg)
  coef <-  reg2$coefficients[,"Estimate"]
  sds <- reg2$coefficients[,"Std. Error"]
  t_alpha <- coef[1]/sds[1]
  t_beta <- (coef[2] - 1)/sds[2]
  t_test <- qt((1+conf)/2,length(y) - 2)
  
  # teste para alpha
  if(abs(t_alpha) > t_test){msg1 <- paste0("   Conclusão: rejeita-se H0 com ",conf*100,"% de confiança")
  }else{msg1 <- paste0("   Conclusão: não rejeita-se H0 com ",conf*100,"% de confiança")}
  
  # teste para beta
  if(abs(t_beta) > t_test){msg2 <- paste0("   Conclusão: rejeita-se H0 com ",conf*100,"% de confiança")
  }else{msg2 <- paste0("   Conclusão: não rejeita-se H0 com ",conf*100,"% de confiança")}
  
  # teste conjunto
  r <- c(0,1)
  R <- rbind(c(1,0),c(0,1))
  pvalor3 <- round(linearHypothesis(reg,c("(Intercept) = 0", "x = 1"))$`Pr(>F)`[2], digits = 6)
  msg3 <- "H0: alpha = 0 e beta = 1"
  cat("- ALPHA", "   H0: alpha = 0", paste("   coef:", round(coef[1],4)),
      paste("   sd:", round(sds[1],4)), msg1, " ",
      "- BETA", "   H0: beta = 1", paste("   coef:", round(coef[2],4)),
      paste("   sd:", round(sds[2],4)), msg2, " ",
      "- ALPHA E BETA", "   H0: alpha = 0 e beta = 1",
      paste("   p-valor:",pvalor3), " ", sep = "\n")
  names(coef)[2] <- deparse(substitute(x))
  b <- list(coef = coef, sd = sds)
}

REQM <- function(y,x){
  dados <- na.omit(cbind(y,x))
  sqrt(mean((dados[,1] - dados[,2])^2, na.rm = T))
}

MAPE <- function(y,x){
  dados <- cbind(y,x)
  mean(abs((dados[,1] - dados[,2])/dados[,1]), na.rm = T)
}