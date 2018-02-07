
# pacotes, funções e leitura
library(dygraphs)
library(BETS)
library(TSA)
library(MARSS)
library(zoo)
library(MASS)
library(numDeriv)
source("functions/dcs_fk_estimation.R")
source("functions/dcs_fk_estimation_exercise.R")
source("functions/ERRO.R")
source("functions/bsm.R")
source("functions/diag.dcs.R", encoding = "utf8")

# leitura -----------
ipc <- window(readRDS("data/ipc.rds"), start = c(2003,7), freq = 12)
nucleo_tf <- readRDS("data/nucleo_tf.rds")

# obter parâmetros para entrar na função de otimização -----------
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
initial_gamma <- c(0.45928571,-0.17857143, 0.05928571,-0.03142857,-0.10214286,-0.38642857,-0.36428571,-0.37071429,-0.34571429,-0.22857143,-0.07857143)
initial_gamma <- c(0.46971429,-0.16814286, 0.06971429,-0.02100000,-0.09171429,-0.37600000,-0.34600000,-0.36600000,-0.30066667,-0.22800000,-0.08000000)

# BSM padrão : mu (beta[t]) + gamma (normal) ------------------------------------

parametros1_normal <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","beta[1|0]","mu[1|0]",paste0("gamma",1:11)          ),
    value = c(0.1 ,0.5 ,0.5 ,5   ,0          ,0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0   ,0   ,0   ,-Inf,-Inf       ,-Inf     ,rep(-Inf,11)                  ),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf        ,Inf      ,rep(Inf,11)                   )
  ),
  Dummy = NA
)
parametros1_normal

dcs1_normal <- dcs_fk_estimation(ipc, initial = parametros1_normal, type = "BSM1_normal", outlier = F, otimo = T, parinitial = F)
dcs1_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs1_normal$hessian)))

data.frame(name = parametros1_normal$par$name, 
           lower = parametros1_normal$par$lower,
           upper = parametros1_normal$par$upper,
           initial = round(parametros1_normal$par$value,4), 
           otimo = round(dcs1_normal$otimizados$par,4))
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
    lower = c(0   ,0   ,0   ,-Inf,-Inf  ,-Inf     ,rep(-Inf,11)                  ,-Inf),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf   ,Inf      ,rep(Inf,11)                   , Inf)
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2002,11)))
)
parametros1_normald

dcs1_normald <- dcs_fk_estimation(ipc, initial = parametros1_normald, type = "BSM1_normal", outlier = T, otimo = T)
dcs1_normald$otimizados$par/sqrt(diag(MASS::ginv(dcs1_normald$hessian)))

data.frame(name = parametros1_normald$par$name, 
           lower = parametros1_normald$par$lower,
           upper = parametros1_normald$par$upper,
           initial = round(parametros1_normald$par$value,4), 
           otimo = round(dcs1_normald$otimizados$par,4))
ts.plot(ipc,dcs1_normald$out[,"mu"], col = 1:2)

ts.plot(dcs1_normald$out[,"epsilon"], col = 1)
round(dcs1_normald$out[,"epsilon"],2)
dcs1_normal$out[,"beta"]

diag_dcs1_normald <- diag.dcs(out = dcs1_normald, type = "norm")
diag_dcs1_normald$stats

# BSM padrão : mu (sem beta) + gamma (normal) ------------------------------------

parametros2_normal <- list(
  par = data.frame(
    name =  c("k1","ks","f2","mu[1|0]",paste0("gamma",1:11)          ),
    value = c(0.1 ,0.5 ,5   ,0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0   ,0   ,-Inf,-Inf     ,rep(-Inf,11)                  ),
    upper = c(Inf ,Inf ,Inf ,Inf      ,rep(Inf,11)                   )
  ),
  Dummy = NA
)
parametros2_normal

dcs2_normal <- dcs_fk_estimation(ipc, initial = parametros2_normal, type = "BSM2_normal", outlier = F, otimo = T)
dcs2_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs2_normal$hessian)))


data.frame(name = parametros2_normal$par$name, 
           lower = parametros2_normal$par$lower,
           upper = parametros2_normal$par$upper,
           initial = round(parametros2_normal$par$value,4), 
           otimo = round(dcs2_normal$otimizados$par,4))
ts.plot(ipc,dcs2_normal$out[,"mu"], col = 1:2)

ts.plot(dcs2_normal$out[,"epsilon"], col = 1)
round(dcs2_normal$out[,"epsilon"],2)

diag_dcs2_normal <- diag.dcs(out = dcs2_normal, type = "norm")
diag_dcs2_normal$stats


# BSM padrão : mu(beta[t]) + gamma ------------------------------------

parametros1 <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","df","beta[1|0]","mu[1|0]",paste0("gamma",1:11)),
    value = c(0.1 ,0.5 ,0.4 ,-2  ,6  ,-0.2        ,0.5    , as.vector(initial_gamma)[1:11]),
    lower = c(0.0 ,0.0 ,0.0 ,-Inf,4   ,-Inf       ,-Inf     ,rep(-Inf,11)),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf ,Inf        ,Inf      ,rep(Inf,11))
  )
)
parametros1
dcs1 <- dcs_fk_estimation(ipc, initial = parametros1, type = "BSM1", outlier = F, otimo = T, parinitial = F)
dcs1$otimizados$par/sqrt(diag(MASS::ginv(dcs1$hessian)))

data.frame(name = parametros1$par$name, 
           lower = parametros1$par$lower,
           upper = parametros1$par$upper,
           initial = round(parametros1$par$value,4), 
           otimo = round(dcs1$otimizados$par,4))
ts.plot(ipc,dcs1$out[,"mu"], col = 1:2)
ts.plot(dcs1$out[,"epsilon"], col = 1)
round(dcs1$out[,"epsilon"],2)
dcs1$out[,"beta"]

diag_dcs1 <- diag.dcs(out = dcs1, type = "t")
diag_dcs1$stats

# BSM padrão : mu (sem beta) + gamma ------------------------------------

parametros2 <- list(
  par = data.frame(
    name =  c("k1","ks","f2","df","mu[1|0]",paste0("gamma",1:11)),
    value = c(0.1 ,0.2 ,5   ,6   ,0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0   ,0.0 ,-Inf,4   ,-Inf     ,rep(-Inf,11)),
    upper = c(Inf ,Inf ,Inf ,Inf ,Inf      ,rep(Inf,11))
  ),
  Dummy = NA
)
parametros2

dcs2 <- dcs_fk_estimation(ipc, initial = parametros2, type = "BSM2", outlier = F, otimo = T, parinitial = F)
dcs2$otimizados$par/sqrt(diag(MASS::ginv(dcs2$hessian)))

data.frame(name = parametros2$par$name, 
           lower = parametros2$par$lower,
           upper = parametros2$par$upper,
           initial = round(parametros2$par$value,4), 
           otimo = round(dcs2$otimizados$par,4))
ts.plot(ipc,dcs2$out[,"mu"], col = 1:2)
ts.plot(dcs2$out[,"mu"],nucleo_tf[,3], col = 1:2)
ts.plot(dcs2$out[,"gamma"], col = 1)
ts.plot(dcs2$out[,"epsilon"], col = 1)
ts.plot(dcs2$out[,"u"], col = 1)
round(dcs2$out[,"epsilon"],2)

diag_dcs2 <- diag.dcs(out = dcs2, type = "t")
diag_dcs2$stats

# BSM padrão : mu (sem beta) + gamma + psi (normal) ------------------------------------

parametros3_normal <- list(
  par = data.frame(
    name =  c("k1","k2","ks","f2","df","beta","mu[1|0]",paste0("gamma",1:11)          ,"psi","phi","k3"),
    value = c(0.5 ,0   ,0.5 ,1   ,0   ,0     ,0.5        ,as.vector(initial_gamma)[1:11],0.2    ,0.5  ,0.2),
    lower = c(0.0 ,0   ,0.05,-Inf,0   ,0     ,-Inf     ,rep(-Inf,11)                  ,-Inf ,-1   ,0),
    upper = c(Inf ,0   ,Inf ,Inf ,0   ,0     ,Inf      ,rep(Inf,11)                   ,Inf  ,1    ,Inf)
  ),
  Dummy = NA
)
parametros3_normal

dcs3_normal <- dcs_fk_estimation(ipc, initial = parametros3_normal, type = "BSM3_normal", outlier = F, otimo = T)
dcs3_normal$otimizados$par/sqrt(diag(MASS::ginv(dcs3_normal$hessian)))


data.frame(name = parametros3_normal$par$name, 
           lower = parametros3_normal$par$lower,
           upper = parametros3_normal$par$upper,
           initial = round(parametros3_normal$par$value,4), 
           otimo = round(dcs3_normal$otimizados$par,4))
ts.plot(ipc,dcs3_normal$out[,"mu"], col = 1:2)

ts.plot(dcs3_normal$out[,"epsilon"], col = 1)
round(dcs3_normal$out[,"epsilon"],2)

diag_dcs3_normal <- diag.dcs(out = dcs3_normal, type = "norm")
diag_dcs3_normal$stats

# adicionar dummy

parametros3_normald <- list(
  par = data.frame(
    name =  c("k1","ks","f2","phi","k3","psi[1|0]","mu[1|0]",paste0("gamma",1:11)          ,"d1"),
    value = c(0.5 ,0.5 ,5   ,0.1  ,0.5 ,1         ,0        ,as.vector(initial_gamma)[1:11],0   ),
    lower = c(0.0 ,0.05,-Inf,-1   ,0   ,-Inf      ,-Inf     ,rep(-Inf,11)                  ,-Inf),
    upper = c(Inf ,Inf ,Inf ,1    ,Inf ,Inf       ,Inf      ,rep(Inf,11)                   ,Inf )
  ),
  Dummy = cbind(d1 = BETS.dummy(start = start(ipc), end = end(ipc), freq = 12, date = c(2016,1)))
)
parametros3_normald

dcs3_normald <- dcs_fk_estimation(ipc, initial = parametros3_normald, type = "BSM3_normal", outlier = T, otimo = T)

data.frame(name = parametros3_normald$par$name, 
           lower = parametros3_normald$par$lower,
           upper = parametros3_normald$par$upper,
           initial = round(parametros3_normald$par$value,4), 
           otimo = round(dcs3_normald$otimizados$par,4))
ts.plot(ipc,dcs3_normald$out[,"mu"], col = 1:2)
ts.plot(ipc,dcs3_normald$out[,"psi"], col = 1:2)

ts.plot(dcs3_normald$out[,"epsilon"], col = 1)
round(dcs3_normald$out[,"epsilon"],2)
dcs3_normal$out[,"beta"]

diag_dcs3_normald <- diag.dcs(out = dcs3_normald, type = "norm")
diag_dcs3_normald$stats


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
    value = c(0.1 ,0.3 ,5   ,6   ,0.1  ,0.5 ,1         ,0        ,as.vector(initial_gamma)[1:11]),
    lower = c(0.0 ,0.05,-Inf,4   ,-1   ,0   ,-Inf      ,-Inf     ,rep(-Inf,11)                  ),
    upper = c(Inf ,Inf ,Inf ,Inf ,1    ,Inf ,Inf       ,Inf      ,rep(Inf,11)                   )
  ),
  Dummy = NA
)
parametros3

dcs3 <- dcs_fk_estimation(ipc, initial = parametros3, type = "BSM3", outlier = F, otimo = T)
dcs3$otimizados$par/sqrt(diag(MASS::ginv(dcs3$hessian)))


parametros3.2 <- list(
  par = data.frame(
    name =  c("k1","ks","f2","df","phi","k3"),
    value = c(0.1 ,0.3 ,5   ,6   ,0.1  ,0.5 ),
    lower = c(0.0 ,0.05,-Inf,4   ,-1   ,0   ),
    upper = c(Inf ,Inf ,Inf ,Inf ,1    ,Inf )
  ),
  par2 = data.frame(
    name =  c("psi[1|0]","mu[1|0]",paste0("gamma",1:11)          ),
    value = dcs3$otimizados$par[-c(1:6)],
    lower = c(-Inf      ,-Inf     ,rep(-Inf,11)                  ),
    upper = c(Inf       ,Inf      ,rep(Inf,11)                   )
  ),
  Dummy = NA
)
parametros3.2

dcs3.2 <- dcs_fk_estimation(ipc, initial = parametros3.2, type = "BSM3", outlier = F, otimo = T, parinitial = T)
dcs3.2$otimizados$par/sqrt(diag(MASS::ginv(dcs3.2$hessian)))
cbind(dcs3$otimizados$par[1:6],dcs3.2$otimizados$par)
data.frame(name = parametros3$par$name, 
           lower = parametros3$par$lower,
           upper = parametros3$par$upper,
           initial = round(parametros3$par$value,4), 
           otimo = round(dcs3$otimizados$par,4))
ts.plot(dcs3$out[,"gamma"], col = 1:2)
ts.plot(dcs3$out[,"psi"], col = 1:2)
ts.plot(dcs3$out[,"mu"], col = 1:2)
ts.plot(dcs3$out[,"mu"], nucleo_tf[,3], col = 1:2)
ts.plot(dcs3$out[,"epsilon"], col = 1:2)
round(dcs3$out[,"epsilon"],2)

diag_dcs3 <- diag.dcs(out = dcs3, type = "t")
diag_dcs3$stats



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
ts.plot(dcs3_normald$out[,"mu"], dcs3$out[,"mu"], nucleo_tf[,3], col = 1:2)
ts.plot(ipc, dcs3$out[,"mu"], col = 1:2)


saveRDS(dcs3_d$out[,"mu"], "data/nucleo_dcs.rds")
saveRDS(dcs3_normald$out[,"mu"], "data/nucleo_dcs_normal.rds")


# gráficos --------------------

par(mar = c(2,4,1,2), mfrow = c(1,1))
# IPC-Br vs. núcleo dcs
plot(ipc, main = "", lwd = 1, lty = 4, ylim = c(-0.5,2),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(dcs3_normald$out[,"mu"],lwd = 2, lty = 5, col = "#1874CD")
lines(dcs3$out[,"mu"], lwd = 2, lty = 1, col = "#CD0000")
abline(h = seq(-0.5,2,0.5),v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,2, legend = c("IPC-Br","DCS-N dummy", "DCS-t"), lwd = c(1,2,2), lty = c(4,5,1), y.intersp = 1.5,
       col = c(1,"#1874CD","#CD0000"), cex = 1.2,bg = "white", box.col = "white",box.lwd = 0)

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


par(mar = c(2,4,1,2), mfrow = c(1,3))

plot(dcs1_d$out[,"mu"], main = "", lwd = 1, lty = 1, ylim = c(-0.2,1.7),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(dcs1_normald$out[,"mu"], lwd = 1, lty = 4, col = "orangered")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.7, legend = c("DCS-Normal 1","DCS-t 1"), lwd = c(1,1), lty = c(1,4), y.intersp = 1.5,
       col = c("orangered",1), cex = 1,bg = "white", box.col = "white",box.lwd = 0)

plot(dcs2_d$out[,"mu"], main = "", lwd = 1, lty = 1, ylim = c(-0.2,1.7),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(dcs2_normald$out[,"mu"], lwd = 1, lty = 4, col = "orangered")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.7, legend = c("DCS-Normal 2","DCS-t 2"), lwd = c(1,1), lty = c(1,4), y.intersp = 1.5,
       col = c("orangered",1), cex = 1,bg = "white", box.col = "white",box.lwd = 0)

plot(dcs3$out[,"mu"], main = "", lwd = 1, lty = 1, ylim = c(-0.2,1.7),
     col = 1, ylab = "variação mensal percentual (%)", xlab = "")
lines(dcs3_normald$out[,"mu"], lwd = 1, lty = 4, col = "orangered")
lines(ipc, lwd = 1, lty = 3, col = "orangered")
abline(h = seq(-0.5,3.5,0.5), col = "#C9C9C9", lty = 3)
abline(v = 1999:2018, col = "#C9C9C9", lty = 3)
legend(2005,1.7, legend = c("DCS-Normal 3","DCS-t 3"), lwd = c(1,1), lty = c(1,4), y.intersp = 1.5,
       col = c("orangered",1), cex = 1,bg = "white", box.col = "white",box.lwd = 0)

