
# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
library(MARSS)
library(zoo)
library(numDeriv)
source("functions/dcs_fk_estimation.R")
source("functions/dcs_fk_estimation_exercise.R")
source("functions/ERRO.R")
source("functions/bsm.R")
source("functions/diag.dcs.R", encoding = "utf8")

# leitura -----------
ipc <- window(readRDS("data/ipc.rds"), start = c(2001,1), freq = 12)
nucleo_tf <- readRDS("data/nucleo_tf.rds")

# obter parâmetros para entrar na função de otimização
y = ipc
dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 11),
               BETS::BETS.dummy(start = start(y), end = end(y), month = 12)
)

y_filter <- y#[1:13]
dummy_filter <- dummy#[1:13,]
t <- 1:length(y_filter)
dados <- data.frame(cbind(y_filter,t,dummy_filter))
dados <- data.frame(window(cbind(ipc,dummy),start = start(ipc), freq = 12))
colnames(dados) <- c("y",paste0("D",1:12))
md <- lm(y ~ ., data = dados)
data.frame(as.vector(md$coefficients)[2:12])

# initial gamma (só para inicializar o filtro, obtido via regressão linear)
# initial_gamma <- c(0.48070729,-0.05859547, 0.15370156, 0.10032539,-0.03874261,-0.26932386,-0.12546543,-0.18050091,-0.20419881,-0.07149072, 0.06560687)
# initial_gamma <- c(0.0125, 0.0775,-0.1750,-0.0275, 0.2600,-0.2025,-0.1050, 0.7225,-0.1100,-0.5425, 0.0350)
initial_gamma <- c(0.418235294,-0.215294118,-0.011764706,-0.048235294,-0.183529412,-0.445294118,-0.320000000,-0.401764706,-0.374705882,-0.247647059, 0.008823529)

# BSM padrão : mu (beta[t]) + gamma (normal) ------------------------------------

parametros1_normal <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","beta[1|0]","mu[1|0]",paste0("gamma",1:11)          ),
    value = c(0.1 ,0.5 ,0.5 ,5   ,0          ,0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0   ,0   ,0.05,-Inf,-Inf       ,-Inf     ,rep(-Inf,11)                  ),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf        ,Inf      ,rep(Inf,11)                   )
  ),
  Dummy = NA
)

parametros1_normal
dcs1_normal <- dcs_fk_estimation(ipc, initial = parametros1_normal, type = "BSM1_normal", outlier = F, otimo = T, parinitial = F)
dcs1_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs1_normal$hessian)))

# parametros1.2_normal <- list(
#   par = data.frame(
#     name =  c("k1","k2","ks","f2"),
#     value = c(0.1 ,0.5 ,0.5 ,5   ),
#     lower = c(0   ,0   ,0.05,-Inf),
#     upper = c(Inf ,Inf ,Inf ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("beta[1|0]","mu[1|0]",paste0("gamma",1:11)          ),
#     value = dcs1_normal$otimizados$par[-c(1:4)],
#     lower = c(-Inf       ,-Inf     ,rep(-Inf,11)                  ),
#     upper = c(Inf        ,Inf      ,rep(Inf,11)                   )
#   ),
#   Dummy = NA
# )
# 
# parametros1.2_normal
# dcs1.2_normal <- dcs_fk_estimation(ipc, initial = parametros1.2_normal, type = "BSM1_normal", outlier = F, otimo = T, parinitial = T)
# dcs1.2_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs1.2_normal$hessian)))
# cbind(dcs1_normal$otimizados$par[1:4],dcs1.2_normal$otimizados$par)

k <- data.frame(name = parametros1_normal$par$name, 
           #lower = parametros1_normal$par$lower,
           #upper = parametros1_normal$par$upper,
           initial = round(parametros1_normal$par$value,4), 
           otimo = round(dcs1_normal$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs1_normal$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k
ts.plot(ipc,dcs1_normal$out[,"mu"], col = 1:2)

ts.plot(dcs1_normal$out[,"epsilon"], col = 1)
round(dcs1_normal$out[,"epsilon"],2)

diag_dcs1_normal <- diag.dcs(out = dcs1_normal, type = "norm")
diag_dcs1_normal$stats

# adicionar dummy

parametros1_normald <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","beta[1|0]","mu[1|0]",paste0("gamma",1:11)          ,"d1"),
    value = c(0.1 ,0.5 ,0.5 ,5   ,0     ,0        ,as.vector(initial_gamma)[1:11],0   ),
    lower = c(0   ,0   ,0.05,-Inf,-Inf  ,-Inf     ,rep(-Inf,11)                  ,-Inf),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf   ,Inf      ,rep(Inf,11)                   , Inf)
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros1_normald

dcs1_normald <- dcs_fk_estimation(ipc, initial = parametros1_normald, type = "BSM1_normal", outlier = T, otimo = T, parinitial = F)
dcs1_normald$otimizados$par/sqrt(diag(MASS::ginv(dcs1_normald$hessian)))

# parametros1.2_normald <- list(
#   par = data.frame(
#     name =  c("k1","k2","ks","f2"),
#     value = c(0.1 ,0.5 ,0.5 ,5   ),
#     lower = c(0   ,0   ,0.05,-Inf),
#     upper = c(Inf ,Inf ,Inf ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("beta[1|0]","mu[1|0]",paste0("gamma",1:11)          ,"d1"),
#     value = dcs1_normald$otimizados$par[-c(1:4)],
#     lower = c(-Inf       ,-Inf     ,rep(-Inf,11)                  ,-Inf),
#     upper = c(Inf        ,Inf      ,rep(Inf,11)                   , Inf)
#   ),
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# parametros1.2_normald
# 
# dcs1.2_normald <- dcs_fk_estimation(ipc, initial = parametros1.2_normald, type = "BSM1_normal", outlier = T, otimo = T, parinitial = T)
# dcs1.2_normald$otimizados$par/sqrt(diag(MASS::ginv(dcs1.2_normald$hessian)))

k <- data.frame(name = parametros1_normald$par$name, 
           #lower = parametros1_normald$par$lower,
           #upper = parametros1_normald$par$upper,
           #initial = round(parametros1_normald$par$value,4), 
           otimo = round(dcs1_normald$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs1_normald$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k
ts.plot(ipc,dcs1_normald$out[,"mu"], col = 1:2)

ts.plot(dcs1_normald$out[,"epsilon"], col = 1)
round(dcs1_normald$out[,"epsilon"],2)

diag_dcs1_normald <- diag.dcs(out = dcs1_normald, type = "norm")
diag_dcs1_normald$stats

# BSM padrão : mu (sem beta) + gamma (normal) ------------------------------------

parametros2_normal <- list(
  par = data.frame(
    name =  c("k1","ks","f2","mu[1|0]",paste0("gamma",1:11)          ),
    value = c(0.1 ,0.5 ,5   ,0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0   ,0.05,-Inf,-Inf     ,rep(-Inf,11)                  ),
    upper = c(Inf ,Inf ,Inf ,Inf      ,rep(Inf,11)                   )
  ),
  Dummy = NA
)
parametros2_normal

dcs2_normal <- dcs_fk_estimation(ipc, initial = parametros2_normal, type = "BSM2_normal", outlier = F, otimo = T, parinitial = F)
dcs2_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs2_normal$hessian)))

# parametros2.2_normal <- list(
#   par = data.frame(
#     name =  c("k1","ks","f2"),
#     value = c(0.1 ,0.5 ,5   ),
#     lower = c(0   ,0.05,-Inf),
#     upper = c(Inf ,Inf ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("mu[1|0]",paste0("gamma",1:11)          ),
#     value = dcs2_normal$otimizados$par[-c(1:3)],
#     lower = c(-Inf     ,rep(-Inf,11)                  ),
#     upper = c(Inf      ,rep(Inf,11)                   )
#   ),
#   Dummy = NA
# )
# parametros2.2_normal
# 
# dcs2.2_normal <- dcs_fk_estimation(ipc, initial = parametros2.2_normal, type = "BSM2_normal", outlier = F, otimo = T, parinitial = T)
# dcs2.2_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs2.2_normal$hessian)))

k <- data.frame(name = parametros2_normal$par$name, 
           #lower = parametros2_normal$par$lower,
           #upper = parametros2_normal$par$upper,
           #initial = round(parametros2_normal$par$value,4), 
           otimo = round(dcs2_normal$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs2_normal$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k
ts.plot(ipc,dcs2_normal$out[,"mu"], col = 1:2)

ts.plot(dcs2_normal$out[,"epsilon"], col = 1)
round(dcs2_normal$out[,"epsilon"],2)

diag_dcs2_normal <- diag.dcs(out = dcs2_normal, type = "norm")
diag_dcs2_normal$stats

# adicionar dummy

parametros2_normald <- list(
  par = data.frame(
    name =  c("k1","ks","f2","mu[1|0]",paste0("gamma",1:11)          ,"d1"),
    value = c(0.1 ,0.5 ,5   ,0        ,as.vector(initial_gamma)[1:11],0   ),
    lower = c(0   ,0.05,-Inf,-Inf     ,rep(-Inf,11)                  ,-Inf),
    upper = c(Inf ,Inf ,Inf ,Inf      ,rep(Inf,11)                   , Inf)
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros2_normald

dcs2_normald <- dcs_fk_estimation(ipc, initial = parametros2_normald, type = "BSM2_normal", outlier = T, otimo = T)
dcs2_normald$otimizados$par/sqrt(diag(MASS::ginv(dcs2_normald$hessian)))


k <- data.frame(name = parametros2_normald$par$name, 
           #lower = parametros2_normald$par$lower,
           #upper = parametros2_normald$par$upper,
           #initial = round(parametros2_normald$par$value,4), 
           otimo = round(dcs2_normald$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs2_normald$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k
ts.plot(ipc,dcs2_normald$out[,"mu"], col = 1:2)

ts.plot(dcs2_normald$out[,"epsilon"], col = 1)
round(dcs2_normald$out[,"epsilon"],2)

diag_dcs2_normald <- diag.dcs(out = dcs2_normald, type = "norm")
diag_dcs2_normald$stats


# BSM padrão : mu(beta[t]) + gamma ------------------------------------

parametros1 <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","df","beta[1|0]","mu[1|0]",paste0("gamma",1:11)),
    value = c(0.3 ,0.5 ,0.5 ,5   ,6   ,0          ,0        , as.vector(initial_gamma)[1:11]),
    lower = c(0.0 ,0.0 ,0.05,-Inf,4   ,-Inf       ,-Inf     ,rep(-Inf,11)),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf ,Inf        ,Inf      ,rep(Inf,11))
  )
)
parametros1

dcs1 <- dcs_fk_estimation(ipc, initial = parametros1, type = "BSM1", outlier = F, otimo = T, parinitial = F)
dcs1$otimizados$par/sqrt(diag(MASS::ginv(dcs1$hessian)))

# parametros1.2 <- list(
#   par = data.frame(
#     name =  c("k1","k2","ks","f2","df"),
#     value = dcs1$otimizados$par[1:5],
#     lower = c(0.0 ,0.0 ,0.05,-Inf,4   ),
#     upper = c(Inf ,Inf ,Inf ,Inf ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("beta[1|0]","mu[1|0]",paste0("gamma",1:11)),
#     value = c(dcs1$otimizados$par[-c(1:5)]),
#     lower = c(-Inf       ,-Inf     ,rep(-Inf,11)),
#     upper = c(Inf        ,Inf      ,rep(Inf,11))
#   )
# )
# parametros1.2
# 
# dcs1.2 <- dcs_fk_estimation(ipc, initial = parametros1.2, type = "BSM1", outlier = F, otimo = T, parinitial = T)
# cbind(dcs1$otimizados$par[1:5],dcs1.2$otimizados$par)
# dcs1.2$otimizados$par/sqrt(diag(MASS::ginv(dcs1.2$hessian)))
# dcs1$otimizados$par[1:5]/sqrt(diag(MASS::ginv(dcs1$hessian[1:5,1:5])))

k <- data.frame(name = parametros1$par$name, 
           #lower = parametros1$par$lower,
           #upper = parametros1$par$upper,
           #initial = round(parametros1$par$value,4), 
           otimo = round(dcs1$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs1$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k
ts.plot(ipc,dcs1$out[,"mu"],dcs1.2$out[,"mu"], col = c(1,2,4))
ts.plot(dcs1$out[,"mu"],dcs1.2$out[,"mu"], col = c(1,2,4))

ts.plot(dcs1$out[,"epsilon"], col = 1)
round(dcs1$out[,"epsilon"],2)
dcs1$out[,"beta"]

diag_dcs1 <- diag.dcs(out = dcs1, type = "t")
diag_dcs1$stats

diag_dcs1.2 <- diag.dcs(out = dcs1.2, type = "t")
diag_dcs1.2$stats

# adicionar dummy
parametros1_d <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","df","beta[1|0]","mu[1|0]",paste0("gamma",1:11)          ,"d1"),
    value = c(dcs1$otimizados$par,0), #c(0.3 ,0.5 ,0.5 ,5   ,6   ,0          ,0        ,as.vector(initial_gamma)[1:11],0   ),
    lower = c(0.0 ,0.0 ,0.05,-Inf,4   ,-Inf       ,-Inf     ,rep(-Inf,11)                  ,-Inf),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf ,Inf        ,Inf      ,rep(Inf,11)                   , Inf)
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros1_d

dcs1_d <- dcs_fk_estimation(ipc, initial = parametros1_d, type = "BSM1", outlier = T, otimo = T, parinitial = F)
dcs1_d$otimizados$par/sqrt(diag(MASS::ginv(dcs1_d$hessian)))

# parametros1.2_d <- list(
#   par = data.frame(
#     name =  c("k1","k2","ks","f2","df"),
#     value = dcs1_d$otimizados$par[1:5],
#     lower = c(0.0 ,0.0 ,0.05,-Inf,4   ),
#     upper = c(Inf ,Inf ,Inf ,Inf ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("beta[1|0]","mu[1|0]",paste0("gamma",1:11)          ,"d1"),
#     value = dcs1_d$otimizados$par[-c(1:5)],
#     lower = c(-Inf       ,-Inf     ,rep(-Inf,11)                  ,-Inf),
#     upper = c(Inf        ,Inf      ,rep(Inf,11)                   , Inf)
#   ),
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# parametros1.2_d
# 
# dcs1.2_d <- dcs_fk_estimation(ipc, initial = parametros1.2_d, type = "BSM1", outlier = T, otimo = T, parinitial = T)
# dcs1.2_d$otimizados$par[c(1,3,4,5,19)]/sqrt(diag(MASS::ginv(dcs1_d$hessian[c(1,3,4,5,19),c(1,3,4,5,19)])))
# dcs1.2_d$otimizados$par/sqrt(diag(MASS::ginv(dcs1.2_d$hessian)))
# dcs1_d$otimizados$par/sqrt(diag(MASS::ginv(dcs1_d$hessian)))
# cbind(dcs1_d$otimizados$par[1:5],dcs1.2_d$otimizados$par)

k <- data.frame(name = parametros1_d$par$name, 
           #lower = parametros1_d$par$lower,
           #upper = parametros1_d$par$upper,
           #initial = round(parametros1_d$par$value,4), 
           otimo = round(dcs1_d$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs1_d$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k

ts.plot(ipc,dcs1_d$out[,"mu"], col = 1:2)
ts.plot(dcs1_d$out[,"epsilon"], col = 1)
round(dcs1_d$out[,"epsilon"],2)

diag_dcs1_d <- diag.dcs(out = dcs1_d, type = "t")
diag_dcs1_d$stats

# BSM padrão : mu (sem beta) + gamma ------------------------------------

parametros2 <- list(
  par = data.frame(
    name =  c("k1","ks","f2","df","mu[1|0]",paste0("gamma",1:11)),
    value = c(0.5 ,0.5 ,5   ,6   ,0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0   ,0.05,-Inf,4   ,-Inf     ,rep(-Inf,11)),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf      ,rep(Inf,11))
  ),
  Dummy = NA
)
parametros2

dcs2 <- dcs_fk_estimation(ipc, initial = parametros2, type = "BSM2", outlier = F, otimo = T)
dcs2$otimizados$par/sqrt(diag(MASS::ginv(dcs2$hessian)))

# parametros2.2 <- list(
#   par = data.frame(
#     name =  c("k1","ks","f2","df"),
#     value = c(0.5 ,0.5 ,5   ,6   ),
#     lower = c(0   ,0.05,-Inf,4   ),
#     upper = c(Inf ,Inf ,Inf ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("mu[1|0]",paste0("gamma",1:11)),
#     value = dcs2$otimizados$par[-c(1:4)],
#     lower = c(-Inf,rep(-Inf,11)),
#     upper = c(Inf ,rep(Inf,11))
#   ),
#   Dummy = NA
# )
# parametros2.2
# 
# dcs2.2 <- dcs_fk_estimation(ipc, initial = parametros2.2, type = "BSM2", outlier = F, otimo = T, parinitial = T)
# cbind(dcs2$otimizados$par[1:4], dcs2.2$otimizados$par)
# dcs2$otimizados$par/sqrt(diag(MASS::ginv(dcs2$hessian)))
# dcs2.2$otimizados$par/sqrt(diag(MASS::ginv(dcs2.2$hessian)))

k <- data.frame(name = parametros2$par$name, 
           #lower = parametros2$par$lower,
           #upper = parametros2$par$upper,
           #initial = round(parametros2$par$value,4), 
           otimo = round(dcs2$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs2$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k
           
ts.plot(ipc,dcs2$out[,"mu"], col = 1:2)
ts.plot(dcs2$out[,"mu"],nucleo_tf[,3], col = 1:2)
ts.plot(dcs2$out[,"gamma"], col = 1)
ts.plot(dcs2$out[,"epsilon"], col = 1)
ts.plot(dcs2$out[,"u"], col = 1)
round(dcs2$out[,"epsilon"],2)

diag_dcs2 <- diag.dcs(out = dcs2, type = "t")
diag_dcs2$stats

# adicionar dummy

parametros2_d <- list(
  par = data.frame(
    name =  c("k1","ks","f2","df","mu[1|0]",paste0("gamma",1:11),"d1"),
    value = c(0.1 ,0.1 ,-2  ,11  ,0.5      ,as.vector(initial_gamma)[1:11],2   ),
    lower = c(0   ,0.05,-Inf,4   ,-Inf     ,rep(-Inf,11),-Inf),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf      ,rep(Inf,11), Inf)
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros2_d

dcs2_d <- dcs_fk_estimation(ipc, initial = parametros2_d, type = "BSM2", outlier = T, otimo = T)
dcs2_d$otimizados$par/sqrt(diag(MASS::ginv(dcs2_d$hessian)))

# parametros2.2_d <- list(
#   par = data.frame(
#     name =  c("k1","ks","f2","df"),
#     value = c(0.1 ,0.1 ,-2  ,11  ),
#     lower = c(0   ,-Inf,-Inf,4   ),
#     upper = c(Inf ,Inf ,Inf ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("mu[1|0]",paste0("gamma",1:11),"d1"),
#     value = c(dcs2_d$otimizados$par[-c(1:4)]),
#     lower = c(-Inf     ,rep(-Inf,11),-Inf),
#     upper = c(Inf      ,rep(Inf,11), Inf)
#   ),
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# parametros2.2_d
# 
# dcs2.2_d <- dcs_fk_estimation(ipc, initial = parametros2.2_d, type = "BSM2", outlier = T, otimo = T, parinitial = T)
# dcs2.2_d
# dcs2.2_d$otimizados$par/sqrt(diag(MASS::ginv(dcs2.2_d$hessian)))

k <- data.frame(name = parametros2_d$par$name, 
           #lower = parametros2_d$par$lower,
           #upper = parametros2_d$par$upper,
           #initial = round(parametros2_d$par$value,4), 
           otimo = round(dcs2_d$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs2_d$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k

ts.plot(ipc,dcs2_d$out[,"mu"], col = 1:2)

ts.plot(dcs2_d$out[,"epsilon"], col = 1)
round(dcs2_d$out[,"epsilon"],2)

diag_dcs2_d <- diag.dcs(out = dcs2_d, type = "t")
diag_dcs2_d$stats


# BSM padrão : mu (sem beta) + gamma + psi (normal) ------------------------------------

parametros3_normal <- list(
  par = data.frame(
    name =  c("k1","ks","f2","phi","k3","psi[1|0]", "mu[1|0]",paste0("gamma",1:11)          ),
    value = c(0.1 ,0.5 ,5   ,0.1  ,0.5 ,1         , 0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0.0 ,0.05,-Inf,-1   ,0   ,-Inf      , -Inf     ,rep(-Inf,11)                  ),
    upper = c(Inf ,Inf ,Inf ,1    ,Inf ,Inf       , Inf      ,rep(Inf,11)                   )
  ),
  Dummy = NA
)
parametros3_normal

dcs3_normal <- dcs_fk_estimation(ipc, initial = parametros3_normal, type = "BSM3_normal", outlier = F, otimo = T, parinitial = F)
dcs3_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs3_normal$hessian)))

# parametros3.2_normal <- list(
#   par = data.frame(
#     name =  c("k1","ks","f2","phi","k3"),
#     value = c(0.1 ,0.5 ,5   ,0.1  ,0.5 ),
#     lower = c(0.0 ,0.05,-Inf,-1   ,0   ),
#     upper = c(Inf ,Inf ,Inf ,1    ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("psi[1|0]", "mu[1|0]",paste0("gamma",1:11)          ),
#     value = dcs3_normal$otimizados$par[-c(1:5)],
#     lower = c(-Inf      , -Inf     ,rep(-Inf,11)                  ),
#     upper = c(Inf       , Inf      ,rep(Inf,11)                   )
#   ),
#   Dummy = NA
# )
# parametros3.2_normal
# 
# dcs3.2_normal <- dcs_fk_estimation(ipc, initial = parametros3.2_normal, type = "BSM3_normal", outlier = F, otimo = T, parinitial = T)
# dcs3.2_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs3.2_normal$hessian)))
# cbind(dcs3_normal$otimizados$par[1:5],dcs3.2_normal$otimizados$par)

k <- data.frame(name = parametros3_normal$par$name, 
           #lower = parametros3_normal$par$lower,
           #upper = parametros3_normal$par$upper,
           #initial = round(parametros3_normal$par$value,4), 
           otimo = round(dcs3_normal$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs3_normal$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k

ts.plot(ipc,dcs3_normal$out[,"mu"], col = 1:2)

ts.plot(dcs3_normal$out[,"epsilon"], col = 1)
round(dcs3_normal$out[,"epsilon"],2)
dcs3_normal$out[,"beta"]

diag_dcs3_normal <- diag.dcs(out = dcs3_normal, type = "norm")
diag_dcs3_normal$stats

# adicionar dummy

parametros3_normald <- list(
  par = data.frame(
    name =  c("k1","ks","f2","phi","k3","psi[1|0]","mu[1|0]",paste0("gamma",1:11)          ,"d1"),
    value = c(0.1 ,0.5 ,5   ,0.1  ,0.5 ,1         ,0        ,as.vector(initial_gamma)[1:11],0   ),
    lower = c(0.0 ,0.05,-Inf,-1   ,0   ,-Inf      ,-Inf     ,rep(-Inf,11)                  ,-Inf),
    upper = c(Inf ,Inf ,Inf ,1    ,Inf ,Inf       ,Inf      ,rep(Inf,11)                   ,Inf )
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11))
  )
)
parametros3_normald

#parametros3_normald$par$value <- dcs3_normald$otimizados$par
dcs3_normald <- dcs_fk_estimation(ipc, initial = parametros3_normald, type = "BSM3_normal", outlier = T, otimo = T, parinitial = F)
dcs3_normald$otimizados$par/sqrt(diag(MASS::ginv(dcs3_normald$hessian)))

# parametros3.2_normald <- list(
#   par = data.frame(
#     name =  c("k1","ks","f2","phi","k3"),
#     value = c(0.1 ,0.5 ,5   ,0.1  ,0.5 ),
#     lower = c(0.0 ,0.05,-Inf,-1   ,0   ),
#     upper = c(Inf ,Inf ,Inf ,1    ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("psi[1|0]","mu[1|0]",paste0("gamma",1:11)          ,"d1"),
#     value = dcs3_normald$otimizados$par[-c(1:5)],
#     lower = c(-Inf      ,-Inf     ,rep(-Inf,11)                  ,-Inf),
#     upper = c(Inf       ,Inf      ,rep(Inf,11)                   ,Inf )
#   ),
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# parametros3.2_normald
# 
# dcs3.2_normald <- dcs_fk_estimation(ipc, initial = parametros3.2_normald, type = "BSM3_normal", outlier = T, otimo = T, parinitial = T)
# dcs3.2_normald$otimizados$par/sqrt(diag(MASS::ginv(dcs3.2_normald$hessian)))
# cbind(dcs3_normald$otimizados$par[1:5],dcs3.2_normald$otimizados$par)

k <- data.frame(name = parametros3_normald$par$name, 
           #lower = parametros3_normald$par$lower,
           #upper = parametros3_normald$par$upper,
           #initial = round(parametros3_normald$par$value,4), 
           otimo = round(dcs3_normald$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs3_normald$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k
ts.plot(ipc,dcs3_normald$out[,"mu"], col = 1:2)
ts.plot(dcs3_normald$out[,"mu"],dcs3.2_normald$out[,"mu"], col = 1:2)

ts.plot(dcs3_normald$out[,"epsilon"], col = 1)
round(dcs3_normald$out[,"epsilon"],2)
dcs3_normal$out[,"beta"]

diag_dcs3_normald <- diag.dcs(out = dcs3_normald, type = "norm")
diag_dcs3_normald$stats

diag_dcs3.2_normald <- diag.dcs(out = dcs3.2_normald, type = "norm")
diag_dcs3.2_normald$stats


# exportar resíduos

# res <- na.omit(data.frame(data = as.Date(dcs3_normal$out[,"epsilon"]), normal = dcs3_normal$out[,"epsilon"], `normal com dummies` = dcs3_normald$out[,"epsilon"]))
# res$`diferença em módulo` <- (res[,3] - res[,2])
# res$`variação percentual` <- (res[,3] - res[,2])/res[,3]*100
# write.csv2(res, "resíduos DCS-Normal.csv")
# ts.plot(ts(res[,4], start = c(1999,1), freq = 12))


# BSM padrão : mu (sem beta) + gamma + psi ------------------------------------
parametros3 <- list(
  par = data.frame(
    name =  c("k1","ks","f2","df","phi","k3","psi[1|0]","mu[1|0]",paste0("gamma",1:11)          ),
    value = c(0.1 ,0.5 ,5   ,6   ,0.1  ,0.5 ,1         ,0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0.0 ,0.05,-Inf,4   ,-1   ,-Inf,-Inf      ,0        ,rep(-Inf,11)                  ),
    upper = c(Inf ,Inf ,Inf ,Inf ,1    ,Inf ,Inf       ,Inf      ,rep(Inf,11)                   )
  ),
  Dummy = NA
)
parametros3

dcs3 <- dcs_fk_estimation(ipc, initial = parametros3, type = "BSM3", outlier = F, otimo = T, parinitial = F)
dcs3$otimizados$par/sqrt(diag(MASS::ginv(dcs3$hessian)))


# parametros3.2 <- list(
#   par = data.frame(
#     name =  c("k1","ks","f2","df","phi","k3"),
#     value = c(0.1 ,0.5 ,5   ,6   ,0.1  ,0.5 ),
#     lower = c(0.0 ,0.05,-Inf,4   ,-1   ,-Inf),
#     upper = c(Inf ,Inf ,Inf ,Inf ,1    ,Inf )
#   ),
#   par2 = data.frame(
#     name =  c("psi[1|0]","mu[1|0]",paste0("gamma",1:11)          ),
#     value = dcs3$otimizados$par[-c(1:6)],
#     lower = c(-Inf      ,0        ,rep(-Inf,11)                  ),
#     upper = c(Inf       ,Inf      ,rep(Inf,11)                   )
#   ),
#   Dummy = NA
# )
# parametros3.2
# 
# dcs3.2 <- dcs_fk_estimation(ipc, initial = parametros3.2, type = "BSM3", outlier = F, otimo = T, parinitial = T)
# dcs3.2$otimizados$par/sqrt(diag(MASS::ginv(dcs3.2$hessian)))

k <- data.frame(name = parametros3$par$name, 
           #lower = parametros3$par$lower,
           #upper = parametros3$par$upper,
           #initial = round(parametros3$par$value,4), 
           otimo = round(dcs3$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs3$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k
ts.plot(dcs3$out[,"gamma"], col = 1:2)
ts.plot(dcs3$out[,"psi"], col = 1:2)
ts.plot(dcs3$out[,"mu"], col = 1:2)
ts.plot(dcs3$out[,"mu"], nucleo_tf[,3], col = 1:2)
ts.plot(dcs3$out[,"epsilon"], col = 1:2)
round(dcs3$out[,"epsilon"],2)

diag_dcs3 <- diag.dcs(out = dcs3, type = "t")
diag_dcs3$stats

# add dummy

parametros3_d <- list(
  par = data.frame(
    name =  c("k1","ks","f2","df","phi","k3","psi[1|0]","mu[1|0]",paste0("gamma",1:11)       ,"d1"),
    value = c(0.3 ,0.5 ,5   ,8   ,0.6  ,0.5 ,-1   ,0.5        ,as.vector(initial_gamma)[1:11],0   ),
    lower = c(0.0 ,0.05,-Inf,4   ,-1   ,0   ,-Inf ,-Inf     ,rep(-Inf,11)                    ,-Inf),
    upper = c(Inf ,Inf ,Inf ,Inf ,1    ,Inf ,Inf  ,Inf      ,rep(Inf,11)                     ,Inf )
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros3_d

dcs3_d <- dcs_fk_estimation(ipc, initial = parametros3_d, type = "BSM3", outlier = T, otimo = T)
dcs3_d$otimizados$par/sqrt(abs(diag(MASS::ginv(dcs3_d$hessian))))


k <- data.frame(name = parametros3_d$par$name, 
           #lower = parametros3_d$par$lower,
           #upper = parametros3_d$par$upper,
           #initial = round(parametros3_d$par$value,4), 
           otimo = round(dcs3_d$otimizados$par,4),
           desvio = round(sqrt(diag(MASS::ginv(dcs3_d$hessian))),4))
k$stat <- round(k$otimo/k$desvio,4)
k

ts.plot(dcs3_d$out[,"gamma"], col = 1:2)
ts.plot(dcs3_d$out[,"psi"], col = 1:2)
ts.plot(dcs3_d$out[,"mu"], col = 1:2)
ts.plot(dcs3_d$out[,"mu"], nucleo_tf[,3], col = 1:2)
ts.plot(dcs3_d$out[,"epsilon"], col = 1:2)
round(dcs3_d$out[,"epsilon"],2)

diag_dcs3_d <- diag.dcs(out = dcs3_d, type = "t")
diag_dcs3_d$stats

saveRDS(dcs3_d$out[,"mu"], "data/nucleo_dcs.rds")


parametros_otimos <- c(0.1195, 0.0000, 0.0500,-1.3745,23.0410, 0.0000, 0.6920, 0.5998,-0.0366, 0.1132, 0.1348,-0.0234,-0.2988,-0.1661,-0.2255,-0.2229,-0.1013, 0.0470,-0.6622, 0.5447, 0.5387, 1.7590)

# COMPARAÇÃO ----------------------------------------
ts.plot(dcs1_d$out[,"mu"],dcs2_d$out[,"mu"],dcs1_normald$out[,"mu"],dcs2_normald$out[,"mu"])

ts.plot(dcs2$out[,"mu"],dcs2_d$out[,"mu"], col = 1:2)
ts.plot(dcs2_normal$out[,"mu"],dcs2_normald$out[,"mu"], col = 1:2)

ts.plot(dcs1_d$out[,"mu"],dcs2_d$out[,"mu"],dcs3_d$out[,"mu"])


data.frame(name = parametros1_normald$par$name, 
           m1 = c(round(dcs1_normal$otimizados$par,4),NA),
           m2 = round(dcs1_normald$otimizados$par,4)
           )

# modelos DCS-Normal
data.frame(name = parametros1_normal$par$name, 
           lower = parametros1_normal$par$lower,
           upper = parametros1_normal$par$upper,
           initial = round(parametros1_normal$par$value,4), 
           otimo = round(dcs1_normal$otimizados$par,4))
data.frame(name = parametros1_normald$par$name, 
           lower = parametros1_normald$par$lower,
           upper = parametros1_normald$par$upper,
           initial = round(parametros1_normald$par$value,4), 
           otimo = round(dcs1_normald$otimizados$par,4))
data.frame(name = parametros2_normal$par$name, 
           lower = parametros2_normal$par$lower,
           upper = parametros2_normal$par$upper,
           initial = round(parametros2_normal$par$value,4), 
           otimo = round(dcs2_normal$otimizados$par,4))
data.frame(name = parametros2_normald$par$name, 
           lower = parametros2_normald$par$lower,
           upper = parametros2_normald$par$upper,
           initial = round(parametros2_normald$par$value,4), 
           otimo = round(dcs2_normald$otimizados$par,4))
data.frame(name = parametros3_normal$par$name, 
           lower = parametros3_normal$par$lower,
           upper = parametros3_normal$par$upper,
           initial = round(parametros3_normal$par$value,4), 
           otimo = round(dcs3_normal$otimizados$par,4))
data.frame(name = parametros3_normald$par$name, 
           lower = parametros3_normald$par$lower,
           upper = parametros3_normald$par$upper,
           initial = round(parametros3_normald$par$value,4), 
           otimo = round(dcs3_normald$otimizados$par,4))

# modelos DCS-t

data.frame(name = parametros1$par$name, 
           lower = parametros1$par$lower,
           upper = parametros1$par$upper,
           initial = round(parametros1$par$value,4), 
           otimo = round(dcs1$otimizados$par,4))
data.frame(name = parametros1_d$par$name, 
           lower = parametros1_d$par$lower,
           upper = parametros1_d$par$upper,
           initial = round(parametros1_d$par$value,4), 
           otimo = round(dcs1_d$otimizados$par,4))
data.frame(name = parametros2$par$name, 
           lower = parametros2$par$lower,
           upper = parametros2$par$upper,
           initial = round(parametros2$par$value,4), 
           otimo = round(dcs2$otimizados$par,4))
data.frame(name = parametros2_d$par$name, 
           lower = parametros2_d$par$lower,
           upper = parametros2_d$par$upper,
           initial = round(parametros2_d$par$value,4), 
           otimo = round(dcs2_d$otimizados$par,4))
data.frame(name = parametros3$par$name, 
           lower = parametros3$par$lower,
           upper = parametros3$par$upper,
           initial = round(parametros3$par$value,4), 
           otimo = round(dcs3$otimizados$par,4))
data.frame(name = parametros3_d$par$name, 
           lower = parametros3_d$par$lower,
           upper = parametros3_d$par$upper,
           initial = round(parametros3_d$par$value,4), 
           otimo = round(dcs3_d$otimizados$par,4))

ts.plot(dcs2_normald$out[,"mu"], nucleo_tf[,3], col = 1:2)
ts.plot(dcs2_d$out[,"mu"], nucleo_tf[,3], col = 1:2)

# saveRDS(dcs3_d$out[,"mu"], "data/nucleo_dcs.rds")
# saveRDS(dcs3_normald$out[,"mu"], "data/nucleo_dcs_normal.rds")


# gráficos --------------------

par(mar = c(2,4,1,2), mfrow = c(1,1))
# IPC-Br vs. núcleo dcs
plot(ipc, main = "", lwd = 1, lty = 4, ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(window(dcs3_normald$out[,"mu"], end = end(ipc), freq = 12), lwd = 2, lty = 1, col = "#1874CD")
abline(h = seq(0,3,1), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
#abline(h = 0, col = "#CC3232", lty = 2)
legend(2005,3, legend = c("IPC-Br","Núcleo-DCS"), lwd = c(1,2,2), lty = c(4,1), y.intersp = 1.5,
       col = c(1,"#1874CD","#CD0000"), cex = 1.3,bg = "white", box.col = "white",box.lwd = 0)

# as tendências dos três modelos DCS

plot(window(dcs1_normald$out[,"mu"], end = end(ipc), freq = 12), main = "", lwd = 1, lty = 4, # ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(window(dcs2_normald$out[,"mu"], end = end(ipc), freq = 12), lwd = 2, lty = 5, col = "orangered")
lines(window(dcs3_normald$out[,"mu"], end = end(ipc), freq = 12), lwd = 2, lty = 1, col = "#1874CD")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.8, legend = c("DCS-N 1","DCS-N 2","DCS-N 3"), lwd = c(1,2,2), lty = c(4,5,1), y.intersp = 1.5,
       col = c(1,"orangered","#1874CD"), cex = 1.2,bg = "white", box.col = "white",box.lwd = 0)


plot(dcs1_d$out[,"mu"], main = "", lwd = 1, lty = 4, # ylim = c(-0.5,3.5),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(dcs2_d$out[,"mu"], lwd = 2, lty = 5, col = "orangered")
lines(dcs3_d$out[,"mu"], lwd = 2, lty = 1, col = "#1874CD")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.8, legend = c("DCS-Normal 1","DCS-Normal 2","DCS-Normal 3"), lwd = c(1,2,2), lty = c(4,5,1), y.intersp = 1.5,
       col = c(1,"orangered","#1874CD"), cex = 1.2,bg = "white", box.col = "white",box.lwd = 0)


par(mar = c(2,4,1,2), mfrow = c(2,2))
layout(mat = matrix(c(1,2,3,3), byrow = T, ncol = 2))
plot(dcs1_d$out[,"mu"], main = "", lwd = 2, lty = 1, ylim = c(-0.2,1.8),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(dcs1_normald$out[,"mu"], lwd = 2, lty = 5, col = "orangered")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.7, legend = c("DCS-Normal 1","DCS-t 1"), lwd = c(2,2), lty = c(5,1), y.intersp = 1.5,
       col = c("orangered",1), cex = 1,bg = "white", box.col = "white",box.lwd = 0)

plot(dcs2_d$out[,"mu"], main = "", lwd = 2, lty = 1, ylim = c(-0.2,1.8),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(dcs2_normald$out[,"mu"], lwd = 2, lty = 5, col = "orangered")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.7, legend = c("DCS-Normal 2","DCS-t 2"), lwd = c(2,2), lty = c(5,1), y.intersp = 1.5,
       col = c("orangered",1), cex = 1,bg = "white", box.col = "white",box.lwd = 0)

plot(dcs3_d$out[,"mu"], main = "", lwd = 2, lty = 1, ylim = c(-0.2,1.8),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(dcs3_normald$out[,"mu"], lwd = 2, lty = 5, col = "orangered")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.7, legend = c("DCS-Normal 3","DCS-t 3"), lwd = c(2,2), lty = c(5,1), y.intersp = 1.5,
       col = c("orangered",1), cex = 1,bg = "white", box.col = "white",box.lwd = 0)

# Exercícios professor -------------------------------------------------------

# 1.estima o modelo normal sem dummies, guarda a  média condicional (previsão) e o resíduo simples ( observado - previsão)
media_cond <- ts(rowSums(dcs3_normal$out[,c("mu","gamma","psi")]), start = start(dcs3_normal$out), freq = 12)
res_simples <- ipc - media_cond
# 2. estima o modelo normal com dummies e guarda a previsão
media_condd <- ts(rowSums(dcs3_normald$out[,c("mu","gamma","psi","dummy")]), start = start(dcs3_normal$out), freq = 12)
# 3. pega o residuo simples da regressão em 1 e faz uma regressão por MQO usando como variável explicativa as dummies dos outliers
data <- na.omit(cbind(res_simples, parametros3_normald$Dummy))
colnames(data) <- c("residuos","d1","d2")
m <- lm(residuos ~ 0 + d1 + d2, data = data)
summary(m)

# 4. pega os parâmetros estimados em 3 e as dummies e forma a seguinte "previsão corrigida"
# previsão corrigida = previsão sem dummies (obtida em 1) + parametros*dummies (obtidas em 3)
prev_corrigida <- media_cond + m$coefficients[1]*parametros3_normald$Dummy[,1] + m$coefficients[2]*parametros3_normald$Dummy[,2] 

# 5. compara a previsão corrigida com a previsão do modelo normal com dummies (obtido em 2): grafico no tempo e diagrama de dispersão.
ts.plot(prev_corrigida, media_condd, col = 1:2, main = "Previsão")
legend("topright", legend = c("DCS com dummies","Previsão corrigida"), col = 2:1, lty = 1, cex = 0.9, bty = "n")

plot(prev_corrigida, media_condd, main = "Gráfico de dispersão: Previsão", xlab = "corrigida", ylab = "DCS com dummies")

# faz a regressão simples de uma na outra : mostra parâmetros e R2
data2 <- na.omit(cbind(media_condd, prev_corrigida))
m2 <- lm(media_condd ~ prev_corrigida, data = data2)
summary(m2)

# # BSM padrão (sem dummy sem psi) -------------------------------------------
# 
# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta", paste0("gamma",1:11)),
#     value = c(0.1,0.1,-1,12, 0.5,0.1,as.vector(initial_gamma)[1:11]),
#     lower = c(0.1,0,-Inf,5,-Inf,-Inf, rep(-Inf,11)),
#     upper = c(0.2,Inf,Inf,Inf,Inf,Inf,rep(Inf,11))
#   ),
#   gamma = NA,
#   Dummy = NA
# )
# 
# dcs_padrao <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2_beta", outlier = F)
# 
# ts.plot(ipc,dcs_padrao$out[,"mu"], col = 1:2)
# 
# ts.plot(dcs_padrao$out[,"epsilon"], col = 1)
# round(dcs_padrao$out[,"epsilon"],2)
# # dois resíduos mt grandes: julho de 2000 e novembro de 2002
# 
# diag_dcs_padrao <- diag.dcs(out = dcs_padrao, type = "t")
# diag_dcs_padrao$stats
# 
# # BSM padrão (com dummy sem psi) -------------------------------------------
# 
# parametros_dummy1 <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta", paste0("gamma",1:11), "d1", "d2"),
#     value = c(0.1,0.1,-1,12, 0.5,-0.1, as.vector(initial_gamma)[1:11], 1,1),
#     lower = c(0,0,-Inf,4,-Inf,-Inf, rep(-Inf,11), -Inf, -Inf),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,rep(Inf,11), Inf, Inf)
#   ),
#   gamma = NA,
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2000,7)),
#                 d2 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# 
# dcs_padrao_dummy <- dcs_fk_estimation(ipc, initial = parametros_dummy1, type = "BSM2_beta", outlier = T)
# 
# ts.plot(ipc,dcs_padrao_dummy$out[,"mu"], col = 1:2)
# 
# ts.plot(dcs_padrao_dummy$out[,"epsilon"], col = 1)
# round(dcs_padrao_dummy$out[,"epsilon"],2)
# 
# diag_dcs_padrao_dummy <- diag.dcs(out = dcs_padrao_dummy, type = "t")
# 
# # BSM padrão (sem dummy com psi) -------------------------------------------
# 
# parametros_psi <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11)),
#     value = c(0.3,0.1,-1,18, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11]),
#     lower = c(0.05,0,-Inf,5,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11)),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11))
#   ),
#   gamma = NA,
#   Dummy = NA
# )
# 
# dcs_padrao_psi <- dcs_fk_estimation(ipc, initial = parametros_psi, type = "BSM2_beta_psi", outlier = F)
# 
# ts.plot(ipc,dcs_padrao_psi$out[,"mu"], col = 1:2)
# 
# ts.plot(dcs_padrao_psi$out[,"epsilon"], col = 1)
# round(dcs_padrao_psi$out[,"epsilon"],2)
# # dois resíduos mt grandes: julho de 2000 e novembro de 2002
# 
# diag_dcs_padrao_psi <- diag.dcs(out = dcs_padrao_psi, type = "t")
# diag_dcs_padrao_psi$stats
# 
# # BSM padrão (com dummy com psi) -------------------------------------------
# 
# parametros_psi_dummy1 <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11), "d1", "d2"),
#     value = c(0.1,0.1,-1,12, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11], 1,1),
#     lower = c(0,0,-Inf,4,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11), -Inf, -Inf),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11), Inf, Inf)
#   ),
#   gamma = NA,
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2000,7)),
#                 d2 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# 
# dcs_padrao_psi_dummy <- dcs_fk_estimation(ipc, initial = parametros_psi_dummy1, type = "BSM2_beta_psi", outlier = T)
# 
# ts.plot(ipc,dcs_padrao_psi_dummy$out[,"mu"], col = 1:2)
# 
# ts.plot(dcs_padrao_psi_dummy$out[,"epsilon"], col = 1)
# round(dcs_padrao_psi_dummy$out[,"epsilon"],2)
# 
# diag_dcs_padrao_psi_dummy <- diag.dcs(out = dcs_padrao_psi_dummy, type = "t")
# diag_dcs_padrao_psi_dummy$stats
# 
# # BSM padrão (normal) -------------------------------------------
# 
# parametros_normal <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","mu[0]","beta", paste0("gamma",1:11)),
#     value = c(0.1,0.1,-1, 0.5,-0.1, as.vector(initial_gamma)[1:11]),
#     lower = c(0,0,-Inf,-Inf,-Inf, rep(-Inf,11)),
#     upper = c(1,Inf,Inf,Inf,Inf,rep(Inf,11))
#   ),
#   gamma = NA,
#   Dummy = NA
# )
# 
# parametros_psi_normal <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11)),
#     value = c(0.5,0.1,-1, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11]),
#     lower = c(0,0,-Inf,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11)),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11))
#   ),
#   gamma = NA,
#   Dummy = NA
# )
# 
# parametros_normal_dummy <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","mu[0]","beta", paste0("gamma",1:11), "d1", "d2"),
#     value = c(0.1,0.1,-1, 0.5,-0.1, as.vector(initial_gamma)[1:11], 1,1),
#     lower = c(0,0,-Inf,-Inf,-Inf, rep(-Inf,11),-Inf,-Inf),
#     upper = c(1,Inf,Inf,Inf,Inf,rep(Inf,11),Inf,Inf)
#   ),
#   gamma = NA,
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2000,7)),
#                 d2 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# 
# parametros_psi_normal_dummy <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11), "d1", "d2"),
#     value = c(0.5,0.1,-1, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11],1,1),
#     lower = c(0,0,-Inf,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11),-Inf,-Inf),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11),Inf,Inf)
#   ),
#   gamma = NA,
#   Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2000,7)),
#                 d2 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
# )
# 
# dcs_padrao_normal <- dcs_fk_estimation(ipc, initial = parametros_normal, type = "BSM2_beta_norm", outlier = F)
# dcs_padrao_psi_normal <- dcs_fk_estimation(ipc, initial = parametros_psi_normal, type = "BSM2_beta_psi_norm", outlier = F)
# dcs_padrao_normal_dummy <- dcs_fk_estimation(ipc, initial = parametros_normal_dummy, type = "BSM2_beta_norm", outlier = T)
# dcs_padrao_psi_normal_dummy <- dcs_fk_estimation(ipc, initial = parametros_psi_normal_dummy, type = "BSM2_beta_psi_norm", outlier = T)
# 
# ts.plot(ipc,dcs_padrao_normal$out[,"mu"], col = 1:2)
# ts.plot(ipc,dcs_padrao_normal_dummy$out[,"mu"], col = 1:2)
# ts.plot(ipc,dcs_padrao_psi_normal$out[,"mu"], col = 1:2)
# ts.plot(ipc,dcs_padrao_psi_normal_dummy$out[,"mu"], col = 1:2)
# 
# ts.plot(dcs_padrao_normal$out[,"mu"], col = 1:2)
# ts.plot(dcs_padrao_psi_normal$out[,"mu"], col = 1:2)
# ts.plot(dcs_padrao_normal_dummy$out[,"mu"], col = 1:2)
# ts.plot(dcs_padrao_psi_normal_dummy$out[,"mu"], col = 1:2)
# 
# ts.plot(dcs_padrao_normal$out[,"gamma"], col = 1:2)
# ts.plot(dcs_padrao_psi_normal$out[,"gamma"], col = 1:2)
# ts.plot(dcs_padrao_normal_dummy$out[,"gamma"], col = 1:2)
# ts.plot(dcs_padrao_psi_normal_dummy$out[,"gamma"], col = 1:2)
# 
# ts.plot(dcs_padrao_psi_normal$out[,"psi"], col = 1:2)
# ts.plot(dcs_padrao_psi_normal_dummy$out[,"psi"], col = 1:2)
# 
# ts.plot(dcs_padrao_psi_normal$out[,"epsilon"], col = 1)
# round(dcs_padrao_psi_normal$out[,"epsilon"],2)
# ts.plot(dcs_padrao_normal$out[,"epsilon"], col = 1)
# round(dcs_padrao_normal$out[,"epsilon"],2)
# ts.plot(dcs_padrao_psi_normal_dummy$out[,"epsilon"], col = 1)
# round(dcs_padrao_psi_normal_dummy$out[,"epsilon"],2)
# ts.plot(dcs_padrao_normal_dummy$out[,"epsilon"], col = 1)
# round(dcs_padrao_normal_dummy$out[,"epsilon"],2)
# 
# diag_dcs_padrao_normal <- diag.dcs(out = dcs_padrao_normal, type = "norm")
# diag_dcs_padrao_normal$stats
# diag_dcs_padrao_psi_normal <- diag.dcs(out = dcs_padrao_psi_normal, type = "norm")
# diag_dcs_padrao_psi_normal$stats
# diag_dcs_padrao_normal_dummy <- diag.dcs(out = dcs_padrao_normal_dummy, type = "norm")
# diag_dcs_padrao_normal_dummy$stats
# diag_dcs_padrao_psi_normal_dummy <- diag.dcs(out = dcs_padrao_psi_normal_dummy, type = "norm")
# diag_dcs_padrao_psi_normal_dummy$stats
# 
# 
# x <- rbind(diag_dcs_padrao_normal$stats,
#       diag_dcs_padrao_psi_normal$stats,
#       diag_dcs_padrao_normal_dummy$stats,
#       diag_dcs_padrao_psi_normal_dummy$stats)
# 
# rownames(x) <- c("sem psi sem dummy","com psi sem dummy","sem psi com dummy","com psi com dummy")
# x
# ts.plot(dcs_padrao_psi_normal_dummy$out[,"mu"],
#         dcs_padrao_psi_dummy$out[,"mu"], col = 1:2, lty = c(1,1,1), lwd = c(1,1,2))
# 
# legend("top", legend = c("Normal", "t-Student"), lty = 1, col = 1:2, bty = "n", cex = 0.8)
# 
# ts.plot(dcs_padrao_normal_dummy$out[,"mu"],
#         dcs_padrao_dummy$out[,"mu"], col = 1:2, lty = c(1,1,1), lwd = c(1,1,2))
# legend("top", legend = c("Normal", "t-Student"), lty = 1, col = 1:2, bty = "n", cex = 0.8)
# 
# 
# ts.plot(dcs_padrao_psi_normal_dummy$out[,"psi"],
#         dcs_padrao_psi_dummy$out[,"psi"], col = 1:2, lty = c(1,1,1), lwd = c(1,1,2))
# legend("top", legend = c("Normal", "t-Student"), lty = 1, col = 1:2, bty = "n", cex = 0.8)
# 
# 
# k <- round(cbind(normal = c(dcs_padrao_psi_normal_dummy$otimizados$par[1:3],NA,dcs_padrao_psi_normal_dummy$otimizados$par[-c(1:3)]),
#                  t = dcs_padrao_psi_dummy$otimizados$par),4)
# rownames(k) <- parametros_psi_dummy1$par$name
# k
# 
# # IMAGENS (do artigo) ---------------------------
# 
# # BSM PADRAO (sem psi sem dummy)
# 
# # fig 6
# par(mfrow = c(2,2))
# ts.plot(ipc,dcs_padrao$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
# ts.plot(dcs_padrao$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
# ts.plot(dcs_padrao$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
# abline(h=0, lty = 3, col = 2)
# ts.plot(dcs_padrao$out[,c("u")],dcs_padrao$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
# 
# # fig 7
# par(mfrow = c(2,2))
# ts.plot(dcs_padrao$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
# ts.plot(dcs_padrao$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
# acf(dcs_padrao$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
# acf(dcs_padrao$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
# 
# # BSM PADRAO (sem psi com dummy)
# 
# # fig 6
# par(mfrow = c(2,2))
# ts.plot(ipc,dcs_padrao_dummy$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
# ts.plot(dcs_padrao_dummy$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
# ts.plot(dcs_padrao_dummy$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
# abline(h=0, lty = 3, col = 2)
# ts.plot(dcs_padrao_dummy$out[,c("u")],dcs_padrao_dummy$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
# 
# # fig 7
# par(mfrow = c(2,2))
# ts.plot(dcs_padrao_dummy$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
# ts.plot(dcs_padrao_dummy$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
# acf(dcs_padrao_dummy$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
# acf(dcs_padrao_dummy$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
# 
# 
# # BSM PADRAO (com psi sem dummy)
# 
# # fig 6
# par(mfrow = c(2,2))
# ts.plot(ipc,dcs_padrao_psi$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
# ts.plot(dcs_padrao_psi$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
# ts.plot(dcs_padrao_psi$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
# abline(h=0, lty = 3, col = 2)
# ts.plot(dcs_padrao_psi$out[,c("u")],dcs_padrao_psi$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
# 
# # fig 7
# par(mfrow = c(2,2))
# ts.plot(dcs_padrao_psi$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
# ts.plot(dcs_padrao_psi$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
# acf(dcs_padrao_psi$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
# acf(dcs_padrao_psi$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
# 
# 
# # BSM PADRAO (com psi com dummy)
# 
# # fig 6
# par(mfrow = c(2,2))
# ts.plot(ipc,dcs_padrao_psi_dummy$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
# ts.plot(dcs_padrao_psi_dummy$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
# ts.plot(dcs_padrao_psi_dummy$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
# abline(h=0, lty = 3, col = 2)
# ts.plot(dcs_padrao_psi_dummy$out[,c("u")],dcs_padrao_psi_dummy$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
# 
# # fig 7
# par(mfrow = c(2,2))
# ts.plot(dcs_padrao_psi_dummy$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
# ts.plot(dcs_padrao_psi_dummy$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
# acf(dcs_padrao_psi_dummy$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
# acf(dcs_padrao_psi_dummy$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
# 
# 
# # EXPORTAR -------------------------------------------
# 
# todos_dcs <- list(
#   ipc = list(
#     dcs_padrao = dcs_padrao,
#     dcs_padrao_dummy = dcs_padrao_dummy,
#     dcs_padrao_psi = dcs_padrao_psi,
#     dcs_padrao_psi_dummy = dcs_padrao_psi_dummy
#   ),
#   ma = list(
#     dcs_padrao = dcs_padrao_ma,
#     dcs_padrao_dummy = dcs_padrao_dummy_ma,
#     dcs_padrao_psi = dcs_padrao_psi_ma,
#     dcs_padrao_psi_dummy = dcs_padrao_psi_dummy_ma
#   )
# )
# saveRDS(todos_dcs, "shiny_simulacao/data/todos_dcs.rds")
# saveRDS(todos_dcs, "dados/todos_dcs.rds")
# 
# todos_diags <- list(
#   ipc = list(
#     diag_dcs_padrao = diag_dcs_padrao,
#     diag_dcs_padrao_dummy = diag_dcs_padrao_dummy,
#     diag_dcs_padrao_psi = diag_dcs_padrao_psi,
#     diag_dcs_padrao_psi_dummy = diag_dcs_padrao_psi_dummy
#   ),
#   ma = list(
#     diag_dcs_padrao = diag_dcs_padrao_ma,
#     diag_dcs_padrao_dummy = diag_dcs_padrao_dummy_ma,
#     diag_dcs_padrao_psi = diag_dcs_padrao_psi_ma,
#     diag_dcs_padrao_psi_dummy = diag_dcs_padrao_psi_dummy_ma
#   )
# )
# saveRDS(todos_diags, "shiny_simulacao/data/todos_diags.rds")
# saveRDS(todos_diags, "dados/todos_diags.rds")
# 
# 
# tendencias <- cbind(dcs_padrao = dcs_padrao$out[,"mu"],
#                          dcs_padrao_dummy = dcs_padrao_dummy$out[,"mu"],
#                          dcs_padrao_psi = dcs_padrao_psi$out[,"mu"],
#                          dcs_padrao_psi_dummy = dcs_padrao_psi_dummy$out[,"mu"],
#                          dcs_padrao_ma = dcs_padrao_ma$out[,"mu"],
#                          dcs_padrao_dummy_ma = dcs_padrao_dummy_ma$out[,"mu"],
#                          dcs_padrao_psi_ma = dcs_padrao_psi_ma$out[,"mu"],
#                          dcs_padrao_psi_dummy_ma = dcs_padrao_psi_dummy_ma$out[,"mu"]
# )
# saveRDS(tendencias, "shiny_simulacao/data/tendencias.rds")
# saveRDS(tendencias, "dados/tendencias.rds")
# write.csv2(data.frame(data = as.Date(tendencias), round(((tendencias/100+1)^12-1)*100,2)), "ultimas_tendencias.csv", row.names = F)
# 
# # BSM com psi e dummy --------------------------------------------------------
# 
# dummy <- cbind(BETS.dummy(start = start(ipc), end =  end(ipc), date = c(2002,11), freq = 12),
#                BETS.dummy(start = start(ipc), end =  end(ipc), date = c(2000,7), freq = 12))
# 
# parametros_beta_psi <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11), "d1","d2"),
#     value = c(0.1,0.1,-1,12, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11], 2,2),
#     lower = c(0,0,-Inf,4,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11), -Inf,-Inf),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11), Inf,Inf)
#   ),
#   gamma = NA,
#   Dummy = dummy
# )
# 
# parametros_beta_psi <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11), "d1"),
#     value = c(0.1,0.1,-1,12, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11], 2),
#     lower = c(0,0,-Inf,4,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11), -Inf),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11), Inf)
#   ),
#   gamma = NA,
#   Dummy = dummy[,1]
# )
# 
# 
# dcs_betapsidummy <- dcs_fk_estimation(ipc, initial = parametros_beta_psi, type = "BSM2_beta_psi", outlier = T)
# 
# ts.plot(ipc,dcs_betapsidummy$out[,"mu"], col = 1:2)
# ts.plot(ipc,dcs_betapsidummy$out[,"psi_mu"], col = 1:2)
# ts.plot(ipc,dcs_betapsidummy$out[,"psi"], col = 1:2)
# 
# # fig 6
# par(mfrow = c(2,2))
# ts.plot(ipc,dcs_betapsidummy$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
# ts.plot(dcs_betapsidummy$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
# ts.plot(dcs_betapsidummy$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
# abline(h=0, lty = 3, col = 2)
# ts.plot(dcs_betapsidummy$out[,c("u")],dcs_betapsidummy$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
# 
# # fig 7
# par(mfrow = c(2,2))
# ts.plot(dcs_betapsidummy$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
# ts.plot(dcs_betapsidummy$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
# acf(dcs_betapsidummy$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
# acf(dcs_betapsidummy$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
# 
# write.csv2(data.frame(data = as.Date(dcs_betapsidummy$out[,"mu"]), round(((dcs_betapsidummy$out[,"mu"]/100+1)^12-1)*100,2)), "tendencia_betapsi_dummy2.csv", row.names = F)
# 
# # diagnóstico
# diag <- diag.dcs(y = ipc, out = dcs_betapsidummy, psi = T, dummy = T)
# diag$stats
# 
# hist(diag$resid.q, main = "residuo quantilico")
# hist(dcs_betapsidummy$out[,"epsilon"], main = "epsilon")
# 
# 
# # BSM com psi --------------------------------------------------------
# 
# parametros_beta_psi <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta","psi","phi","k3", paste0("gamma",1:11)),
#     value = c(0.1,0.1,-1,12, 0.5,-0.1,0.1,0.1,0.1, as.vector(initial_gamma)[1:11]),
#     lower = c(0,0,-Inf,4,-Inf,-Inf,-Inf,-1,0, rep(-Inf,11)),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,rep(Inf,11))
#   ),
#   gamma = NA
# )
# 
# parametros_beta <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta",paste0("gamma",1:11)),
#     value = c(0.11,0.1,-1,12, 0.5,-0.1, as.vector(initial_gamma)[1:11]),
#     lower = c(0,0,-Inf,4,-Inf,-Inf, rep(-Inf,11)),
#     upper = c(0.11,Inf,Inf,Inf,Inf,Inf,rep(Inf,11))
#   ),
#   gamma = NA
# )
# dcs_betapsi <- dcs_fk_estimation(ipc, initial = parametros_beta_psi, type = "BSM2_beta_psi")
# dcs_beta <- dcs_fk_estimation(ipc, initial = parametros_beta, type = "BSM2_beta")
# dcs_betapsi_ma <- dcs_fk_estimation(ipc_ma1, initial = parametros_beta_psi, type = "BSM2_beta_psi")
# dcs_beta_ma <- dcs_fk_estimation(ipc_ma1, initial = parametros_beta, type = "BSM2_beta")
# 
# ts.plot(ipc,dcs_betapsi$out[,"mu"], col = 1:2)
# ts.plot(ipc,dcs_betapsi$out[,"psi_mu"], col = 1:2)
# ts.plot(ipc,dcs_betapsi$out[,"psi"], col = 1:2)
# 
# ts.plot(dcs_betapsi$out[,"mu"], dcs_beta$out[,"mu"], col = 1:2)
# ts.plot(dcs_betapsi_ma$out[,"mu"], dcs_beta_ma$out[,"mu"], col = 1:2)
# ts.plot(dcs_betapsi$out[,"mu"], dcs_beta$out[,"mu"], dcs_beta_ma$out[,"mu"], col = 1:2)
# ts.plot(dcs_beta_ma$out[,"mu"], ipc_ma3, col = 1:2)
# ts.plot(ipc, dcs_beta_ma$out[,"mu"], ipc_ma3, col = c(1,2,4))
# 
# ts.plot(ipc,dcs_beta$out[,"mu"], col = 1:2)
# # pseudo.y <- (1 - dcs_res$out[,"b"])*ipc + dcs_res$out[,"b"]*(dcs_res$out[,"mu"] + dcs_res$out[,"gamma"])
# # # mu smooth
# # bsm_res <- bsm(pseudo.y, type = "BSM2_beta", iter = 10000)
# # ts.plot(bsm_res[,"mu"],ipc, col = 2:1)
# # 
# # ts.plot(bsm_res[,"mu"], dcs_res$out[,"mu"], col = 2:1)
# # 
# # saveRDS(list(core_res = dcs_res$out[,"mu"], core_ress = bsm_res[,"mu"]), "./dados/nucleo.rds")
# 
# 
# # fig 6
# par(mfrow = c(2,2))
# ts.plot(ipc,dcs_betapsi$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
# ts.plot(dcs_betapsi$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
# ts.plot(dcs_betapsi$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
# abline(h=0, lty = 3, col = 2)
# ts.plot(dcs_betapsi$out[,c("u")],dcs_betapsi$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
# 
# # fig 7
# par(mfrow = c(2,2))
# ts.plot(dcs_betapsi$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
# ts.plot(dcs_betapsi$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
# acf(dcs_betapsi$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
# acf(dcs_betapsi$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
# 
# 
# # fig 6
# par(mfrow = c(2,2))
# ts.plot(ipc,dcs_beta$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
# ts.plot(dcs_beta$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
# ts.plot(dcs_beta$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
# abline(h=0, lty = 3, col = 2)
# ts.plot(dcs_beta$out[,c("u")],dcs_beta$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
# 
# # fig 7
# par(mfrow = c(2,2))
# ts.plot(dcs_beta$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
# ts.plot(dcs_beta$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
# acf(dcs_beta$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
# acf(dcs_beta$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
# 
# 
# # fig 8
# # ts.plot(fk2[,"mu"],bsm2$out[,"mu"], col = 1, lty = c(1,3))
# # ts.plot(fk2_beta[,"mu"],bsm2_beta$out[,"mu"], col = 1, lty = c(1,3))
# # ts.plot(fk2[,"mu"],fk2_beta[,"mu"], col = 1, lty = c(1,3))
# # ts.plot(fk2[,"mu"],ipc, col = 1, lty = c(1,3))
# # ts.plot(fk2_beta[,"mu"],ipc, col = 1, lty = c(1,3), lwd = c(2,1))
# 
# write.csv2(data.frame(data = as.Date(dcs_betapsi$out[,"mu"]), round(((dcs_betapsi$out[,"mu"]/100+1)^12-1)*100,2)), "tendencia_betapsi.csv", row.names = F)
# 
# # diagnóstico
# diag <- diag.dcs(y = ipc, out = dcs_betapsi, psi = T)
# 
# hist(diag$resid.q, main = "residuo quantilico")
# hist(dcs_res$out[,"epsilon"], main = "residuo de pearson")
# 
# # BSM para IPC ----------------------------------------------------
# # valores iniciais
# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta"),
#     value = c(0.02,0.5,0.5,10, 0.49,-0.01),
#     lower = c(0,0.001,-Inf,4,-Inf,-Inf),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf)
#   ),
#   gamma = as.vector(initial_gamma)[1:11]
# )
# 
# # IRRESTRITO
# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta", paste0("gamma",1:11)),
#     value = c(0.1,0.1,-1,12, 0.5,-0.1, as.vector(initial_gamma)[1:11]),
#     lower = c(0,0,-Inf,4,-Inf,-Inf, rep(-Inf,11)),
#     upper = c(1,Inf,Inf,Inf,Inf,Inf,rep(Inf,11))
#   ),
#   gamma = NA
# )
# parametros
# dcs_full <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2_beta")
# ts.plot(ipc,dcs_full$out[,"mu"], col = 1:2)
# 
# pseudo.y <- (1 - dcs_full$out[,"b"])*ipc + dcs_full$out[,"b"]*(dcs_full$out[,"mu"] + dcs_full$out[,"gamma"])
# 
# # mu smooth
# bsm_full <- bsm(pseudo.y, type = "BSM2_beta", iter = 10000)
# ts.plot(bsm_full[,"mu"],ipc, col = 2:1)
# 
# ts.plot(bsm_full[,"mu"], dcs_full$out[,"mu"], col = 2:1)
# 
# # RESTRITO
# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]","beta", paste0("gamma",1:11)),
#     value = c(0.1,0.1,-1,12, 0.5,-0.1, as.vector(initial_gamma)[1:11]),
#     lower = c(0,0,-Inf,4,-Inf,-Inf, rep(-Inf,11)),
#     upper = c(0.1,Inf,Inf,Inf,Inf,Inf,rep(Inf,11))
#   ),
#   gamma = NA
# )
# parametros
# dcs_res <- dcs_fk_estimation(ipc, initial = parametros, type = "BSM2_beta")
# ts.plot(ipc,dcs_res$out[,"mu"], col = 1:2)
# ts.plot(dcs_res$out[,"gamma"], col = 1:2)
# pseudo.y <- (1 - dcs_res$out[,"b"])*ipc + dcs_res$out[,"b"]*(dcs_res$out[,"mu"] + dcs_res$out[,"gamma"])
# 
# x <- data.frame(parametros$par, otimo = round(dcs_res$otimizados$par,4))
# x[2,3:5] <- NA 
# x
# # mu smooth
# bsm_res <- bsm(pseudo.y, type = "BSM2_beta", iter = 10000)
# ts.plot(bsm_res[,"mu"],ipc, col = 2:1)
# 
# ts.plot(bsm_res[,"mu"], dcs_res$out[,"mu"], col = 2:1)
# 
# saveRDS(list(core_res = dcs_res$out[,"mu"], core_ress = bsm_res[,"mu"]), "./dados/nucleo.rds")
# 
# 
# # fig 6
# par(mfrow = c(2,2), mar = c(3,3,2,3))
# ts.plot(ipc,dcs_res$out[,"mu"], col = 1, lwd = c(1,2), main = "Y e Tendência", ylab = "")
# ts.plot(dcs_res$out[,"mu"], col = 1, lwd = c(1), main = "Tendência", ylab = "")
# ts.plot(dcs_res$out[,"gamma"], col = 1, lwd = c(1,2), main = "Sazonalidade", ylab = "")
# abline(h=0, lty = 3, col = 2)
# ts.plot(dcs_res$out[,c("u")],dcs_res$out[,c("nu")], col = 1, lty = c(1,3), main = "u e nu = exp(f2)*epsilon.", ylab = "")
# 
# # fig 7
# par(mfrow = c(2,2), mar = c(3,3,2,3))
# ts.plot(dcs_res$out[,c("nu")], col = 1, lty = c(1,3), main = "nu = exp(f2)*epsilon.", ylab = "", ylim = c(-1,2.1))
# ts.plot(dcs_res$out[,c("u")], col = 1, lty = c(1,3), main = "u", ylab = "", ylim = c(-1,2.1))
# acf(dcs_res$out[,c("nu")], 48, drop.lag.0 = T, main = "acf nu")
# acf(dcs_res$out[,c("u")], 48, drop.lag.0 = T, main = "acf u")
# 
# # fig 8
# # ts.plot(fk2[,"mu"],bsm2$out[,"mu"], col = 1, lty = c(1,3))
# # ts.plot(fk2_beta[,"mu"],bsm2_beta$out[,"mu"], col = 1, lty = c(1,3))
# # ts.plot(fk2[,"mu"],fk2_beta[,"mu"], col = 1, lty = c(1,3))
# # ts.plot(fk2[,"mu"],ipc, col = 1, lty = c(1,3))
# # ts.plot(fk2_beta[,"mu"],ipc, col = 1, lty = c(1,3), lwd = c(2,1))
# 
# write.csv2(data.frame(data = as.Date(bsm2_beta$out[,"mu"]), round(((bsm2_beta$out[,"mu"]/100+1)^12-1)*100,2)), "tendencia5.csv", row.names = F)
# 
# # diagnóstico
# diag <- diag.dcs(y = dcs_res$out[,"epsilon"], df = dcs_res$otimizados$par[4])
# 
# hist(diag$resid.q, main = "modelo sem beta")
# 
# diag$stats
# 
# saveRDS(bsm2, "shiny_simulacao/data/bsm_ipc.rds")
# saveRDS(pseudo.y, "shiny_simulacao/data/pseudoy_ipc.rds")
# saveRDS(fk2, "shiny_simulacao/data/smooth_ipc.rds")
# saveRDS(diag, "shiny_simulacao/data/diag_ipc.rds")
# #saveRDS(psd, "shiny_simulacao/data/smooth_ipc_psd.rds")
# 
# # nucleo final
# nucleo <- bsm2$out[,"mu"]
# nucleo <- fk2[,"mu"]
# # exercício
# 
# out <- NULL
# vero <- NULL
# k1 <- NULL
# ks <- seq(0,1,by = 0.01)
# mu <- NULL
# beta <- NULL
# dfs <- 4:30
# for(i in 1:length(ks)){
#   #out[[i]] <- dcs_fk_estimation_exercise(ipc, initial = parametros, type = "BSM2_beta", k1 = k1[i])
#   #out[[i]] <- dcs_fk_estimation_exercise(ipc, initial = parametros, type = "BSM2_beta", df = dfs[i])
#   out[[i]] <- dcs_fk_estimation_exercise(ipc, initial = parametros, type = "BSM2_beta", ks = ks[i])
#   vero[i] <- out[[i]]$loglik
#   
# }
# k <- cbind(ks = ks, loglik = -vero)
# plot(k[,1:2], xlab = "graus de liberdade", ylab = "loglik")
# round(k,3)
# 
# library(zoo)
# exportar <- data.frame(data = as.Date(cbind(bsm2$out, fk2)),cbind(bsm2$out, fk2))
# colnames(exportar) <- c("data", colnames(bsm2$out), paste0(colnames(fk2),"_smoother"))
# write.csv2(exportar,"shiny_simulacao/www/ultimos_resultados.csv", row.names = F)
# 

# lixo ------------------------------------------------------
# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]",paste0("gamma",1:11)),
#     value = c(0.5,-1,-1.22,11,mean(ipc0),c(0.63634560, -0.11379156, 0.08584090, 0.04461767, 0.27043329, -0.11460482, -0.09033392, -0.07185322, -0.51438998, -0.29311867, -0.02459362)),
#     lower = c(0,-Inf,-Inf,2, rep(-Inf,12)),
#     upper = c(1,Inf,Inf,Inf, rep(Inf,12))
#   ),
#   mu = 0.92148341, 
#   gamma = c(0.04130164, 0.10012847, 0.01462955,-0.20743980,-0.39845356,-0.15230706, 0.67272960,-0.02919001,-0.41970979, 0.06692283, 0.27037323)
# )

# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]",paste0("gamma",1:11)),
#     value = c(0.8,-1.1,-1.6,6,mean(ipc),initial_gamma[1:11]),
#     lower = c(0,-Inf,-Inf,2, rep(-Inf,12)),
#     upper = c(1,Inf,Inf,Inf, rep(Inf,12))
#   ),
#   mu = initial_mu,
#   gamma = as.vector(initial_gamma[1:11])
# )
# 
# parametros <- list(
#   par = data.frame(
#     name = c("k1","ks","f2","df","mu[0]",paste0("gamma",1:11)),
#     value = c(0.5,-1,-1.22,11,mean(ipc0),c(0.63634560, -0.11379156, 0.08584090, 0.04461767, 0.27043329, -0.11460482, -0.09033392, -0.07185322, -0.51438998, -0.29311867, -0.02459362)),
#     lower = c(0,-Inf,-Inf,2, rep(-Inf,12)),
#     upper = c(1,Inf,Inf,Inf, rep(Inf,12))
#   ),
#   mu = 0.92148341, 
#   gamma = c(0.04130164, 0.10012847, 0.01462955,-0.20743980,-0.39845356,-0.15230706, 0.67272960,-0.02919001,-0.41970979, 0.06692283, 0.27037323)
# )

# # USGDP
# 
# usgdp <- window(ts(read.csv("dados/GDPC1.csv")[,2], start = c(1947,1), freq = 4), end = c(2012,1), freq = 4)
# y <- diff(log(usgdp))
# m <- dcs_fk_estimation(y, initial = c(0.2,0.2,0.2,4,5), type = "mean")
# round(m$otimizados$par,4)
# ts.plot(y,m$out[,c(1)], col = 1:2, lwd = 1:2)
# 
# 
# # AWHMAN
# aw <- window(ts(read.csv("dados/AWHMAN.csv")[,2], start = c(1939,1), freq = 12), start = c(1992,2), end = c(2010,5), freq = 12)
# 
# m <- dcs_fk_estimation(aw, initial = c(0.2,0.2,4,5,1), type = "cn")
# round(m$otimizados$par,4)
# 
# k <- window(cbind(aw,lag(m$out[,c(1)],1)), end = c(1997,1), freq = 12)
# ts.plot(k, col = 1:2, lwd = 1:2)
# 
# # IPC
# ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
# m1 <- dcs_fk_estimation(ipc, initial = c(0.2,0.2,0.2,4,5), type = "mean")
# m2 <- dcs_fk_estimation(ipc, initial = c(0.2,0.2,4,5,1), type = "cn")
# m3 <- dcs_fk_estimation(ipc, type = "cn2")
# 
# round(m1$otimizados$par,4)
# round(m2$otimizados$par,4)
# ts.plot(ipc,m1$out[,c(1)], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m2$out[,c(1)], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m3$out[,"mu"] + m3$out[,"beta"] + m3$out[,"gamma"], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m3$out[,"mu"], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m3$out[,"mu"] + m3$out[,"beta"], col = 1:2, lwd = 1:2)
# ts.plot(ipc,m3$out[,"beta"], col = 1:2, lwd = 1:2)
# 
# plot(m3$out)
# 
# # initial dummies
# 
# y = ipc0
# dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 11),
#                BETS::BETS.dummy(start = start(y), end = end(y), month = 12)
# )
# 
# y_filter <- y[1:13]
# dummy_filter <- dummy[1:13,]
# t <- 1:13
# dados <- data.frame(cbind(y_filter,t,dummy_filter))
# dados <- data.frame(window(cbind(ipc0,dummy),start = start(ipc0), freq = 12))
# colnames(dados) <- c("y",paste0("D",1:12))
# md <- lm(y ~ ., data = dados)
# data.frame(as.vector(md$coefficients)[2:12])
# c(0.37195906,-0.20488304,-0.02277778,-0.09172515,-0.22804094,-0.45383041,-0.21172515,-0.37804094,-0.39611111,-0.25611111,0.01666667)
# c(0.43032967,-0.20752747, 0.03032967,-0.06038462,-0.13109890,-0.41538462,-0.39324176,-0.39967033,-0.35461538,-0.25538462,-0.09615385)
# c(0.50767857,-0.13357143, 0.07642857,-0.09607143,-0.19357143,-0.54357143,-0.54107143,-0.46982143,-0.29142857,-0.20142857,-0.04857143)
# c(0.42891026,-0.19647436, 0.04583333,-0.04032051,-0.14724359,-0.45878205,-0.42032051,-0.44262821,-0.33250000,-0.23250000,-0.08250000)
# # 
# data <- data.frame(param = c(paste0("D",1:11),            "df","w1","w2","beta","a11","a12","a13","a2","b2","f1","f2"),
#                    value = c(c(md$coefficients[-c(1,2,14)]), 3,   1,   1,  0.10,  0.5,  0.5,  0.5, 0.5, 0.9,   2, 0.5),
#                    lower = c(rep(-Inf,11),                   2,-Inf,-Inf, -Inf,-Inf,-Inf,-Inf,-Inf,-1,  -Inf, -Inf),
#                    upper = c(rep( Inf,11),                 Inf, Inf, Inf,  Inf, Inf, Inf, Inf, Inf, 1,   Inf,  Inf))   
# 
# data$value <- round(c(1.395931e+00, 7.717124e-01, 1.143014e+00, 1.017791e+00,-9.203240e+00, 7.213934e-01, 8.484142e-01, 7.339786e-01, 7.438413e-01, 
#  7.914208e-01, 9.271763e-01, 3.920452e+00,-8.925764e-01,-8.142339e-02,-3.349298e-03, 7.269610e-03, -7.303699e-05, 6.873406e-03,
#  2.295186e-01, 9.511262e-01, 6.880529e-01, 1.036890e+00),4)
# 
# data$value <- m$otimizados$par
# 
# round(m$otimizados$par[1:11], 3)
# round(md$coefficients[-c(1,2,14)], 3)
# 
# data_new <- data.frame(param = c(paste0("D",1:11),            "df","w1","w2","beta","a11","a12","a13","a2","b2","f1","f2"),
#                    value = c(m$otimizados$par))   
# 
# 
# data <- data.frame(param = c(paste0("D",1:11),"df","w1","w2","beta","a11","a12","a13","a2","b2","f1","f2"),
#                    value = c(rep(1,11),          4,   0, 0.1,  0.10,  0.3,  0.3,  0.5, 0.5, 0.1, 0.5, 0.5),
#                    lower = c(rep(-Inf,11),       2,-Inf,-Inf, -Inf,-Inf,-Inf,-Inf,-Inf,-1,  -Inf, -Inf),
#                    upper = c(rep( Inf,11),     Inf, Inf, Inf,  Inf, Inf, Inf, Inf, Inf, 1,   Inf,  Inf))   
# 
# 
# 
# a <- Sys.time()
# m <- dcs_fk_estimation(ipc, type = "cn3", initial = data, dummy = dummy)
# b <- Sys.time()
# b-a
# y=ipc
# ts.plot(y,m$out[,"mu"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência + Sazonalidade")
# ts.plot(y,m$out[,"f1"], col = 1:2, lwd = c(1,2), main = "IPC e Tendência")
# ts.plot(y,m$out[,"f2"], col = 1:2, lwd = c(1,2), main = "IPC e Sazonalidade")
# ts.plot(y,m$out[,"sigma"], col = 1:2, lwd = c(1,2), main = "IPC e Sigma")
# ts.plot(m[,c("f1","f2")], col = 1:2, lwd = c(1,2), main = "Tendência e Sazonalidade")
# 
# 
# # initial_mu <- median(ipc)
# # gammas <- data.frame((ipc0 - mean(ipc0))/sd(ipc0), cycle(ipc0))
# # gammas <- data.frame((ipc - mean(ipc))/sd(ipc), cycle(ipc))
# # gammas <- data.frame(ipc0, cycle(ipc0))
# # initial_gamma <- tapply(gammas[,1],gammas[,2], FUN = mean)
# # 
# # y = ipc0
# # dummy <- cbind(BETS::BETS.dummy(start = start(y), end = end(y), month = 1),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 2),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 3),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 4),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 5),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 6),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 7),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 8),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 9),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 10),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 11),
# #                BETS::BETS.dummy(start = start(y), end = end(y), month = 12)
# # )
# # dados <- cbind(y, dummy)
# # colnames(dados) <- c("y",paste0("D",1:12))
# # md <- lm(y ~ ., data = dados)
# # data.frame(as.vector(md$coefficients)[2:12])
# #initial_gamma <- c(0.37195906,-0.20488304,-0.02277778,-0.09172515,-0.22804094,-0.45383041,-0.21172515,-0.37804094,-0.39611111,-0.25611111, 0.01666667)
# #initial_gamma <- c(0.50767857,-0.13357143, 0.07642857,-0.09607143,-0.19357143,-0.54357143,-0.54107143,-0.46982143,-0.29142857,-0.20142857,-0.04857143)
# initial_gamma <- c(0.05516970,0.09563170,0.05699505,-0.18377041,-0.43230580,-0.17613783,0.69144612,-0.04360480,-0.44064156,0.07471822,0.27193868)
# initial_gamma <- c(0.153357232, 0.017716173, 0.132833222,-0.139159857,-0.404060430,-0.145550980, 0.660517961,-0.048933837,-0.449613456, 0.022034873, 0.223473635)
# #initial_gamma <- c(0.512692417,-0.041189734, 0.144115903, 0.088499102,-0.060261330,-0.266698219,-0.138449064,-0.195460536,-0.221822482,-0.069103134, 0.092489454)
# initial_gamma <- c(0.5193,-0.0469, 0.1511, 0.0878,-0.0606,-0.2625,-0.1520,-0.1967,-0.2275,-0.0687, 0.0995)
# initial_gamma <- c(0.48070729,-0.05859547, 0.15370156, 0.10032539,-0.03874261,-0.26932386,-0.12546543,-0.18050091,-0.20419881,-0.07149072, 0.06560687)

