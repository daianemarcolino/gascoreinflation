# script desenvolvido por Daiane Marcolino de Mattos
# objetivo: simular GAS
# início da programação: 13/08/2017 às 18h26 - domingo

# carregar funções --------------------------------------
source("functions/simularGAS.R")

# simulação 1 -------------------------------------------
# y[t] = mu[t-1] + epsilon[t], epsilon[t] ~ N(0,sigma2)
# theta (parâmetros estáticos): c(sigma2, w, A0, B1)

# simular série
sim1 <- simularGAS(density = "normal", n = 300, mean = "variavel", sd = "fixa", seed = 1, theta = c(0.5, 0, 0.2, 0.4))
sim1 <- simularGAS(density = "normal", n = 300, mean = "variavel", sd = "fixa", seed = 1, theta = c(1, 2, 0.5, 1.5))
plot(sim1)

# estimar GAS
gas1 <- estimarGAS(y = sim1[,"y"], density = "normal", mean = "variavel", sd = "fixa", initial = c(2, 0.5, 0.5, 0.5))
gas1$optim
round(gas1$optim,1)

# Y vs. média estimada
ts.plot(sim1[,"y"], gas1$ft, col = c("#909090","red"), lwd = 1:2,
        main = "w = 0; A0 = 0.2; B1 = 0.4; sigma2 = 0.5")

# média estimada vs. média observada
ts.plot(sim1[,"mu"], gas1$ft, col = c("#909090","red"), lwd = 1,
        main = "w = 0; A0 = 0.2; B1 = 0.4; sigma2 = 0.5")


# simular série
sim2 <- simularGAS(density = "normal", n = 300, mean = "variavel", sd = "fixa", seed = 1, theta = c(0.5, 0, 0.2, 1))
plot(sim2)

# estimar GAS
gas2 <- estimarGAS(y = sim2[,"y"], density = "normal", mean = "variavel", sd = "fixa", initial = c(2, 0.5, 0.5, 0.5))
gas2$optim
round(gas2$optim,1)

ts.plot(sim2[,"y"], gas2$ft, col = c("#909090","red"), lwd = 1:2,
        main = "w = 0; A0 = 0.2; B1 = 1; sigma2 = 0.5")

ts.plot(sim1[,"mu"], gas1$ft, col = c("#909090","red"), lwd = 1,
        main = "w = 0; A0 = 0.2; B1 = 1; sigma2 = 0.5")


# simular série
sim3 <- simularGAS(density = "normal", n = 300, mean = "variavel", sd = "fixa", seed = 1, theta = c(0.5, 1.5, 0.5, 0.7))
plot(sim3)

# estimar GAS
gas3 <- estimarGAS(y = sim3[,"y"], density = "normal", mean = "variavel", sd = "fixa", initial = c(2, 0.5, 0.5, 0.5))
gas3$optim
round(gas3$optim,1)

ts.plot(sim3[,"y"], gas3$ft, col = c("#909090","red"), lwd = 1:2,
        main = "w = 1.5; A0 = 0.5; B1 = 0.7; sigma2 = 0.5")

ts.plot(sim3[,"mu"], gas3$ft, col = c("#909090","red"), lwd = 1,
        main = "w = 1.5; A0 = 0.5; B1 = 0.7; sigma2 = 0.5")


# simulação 2 -------------------------------------------
# y[t] = mu + sigma2[t-1]*epsilon[t], epsilon[t] ~ N(0,1)
# theta (parâmetros estáticos): c(mean y, w, A0, B1)

# simular série
sim1 <- simularGAS(density = "normal", n = 300, mean = "fixa", sd = "variavel", seed = 1, theta = c(1, 0.5, 0.2, 0.6))
plot(sim1)

# estimar GAS
gas1 <- estimarGAS(y = sim1[,"y"], density = "normal", mean = "fixa", sd = "variavel", initial = c(1, 1, 1, 1))
gas1$optim
round(gas1$optim,1)

# Y vs. sigma estimada
ts.plot(sim1[,"y"], gas1$ft, col = c("#909090","red"), lwd = 1:2,
        main = "mean y = 1; w = 0.5; A0 = 0.2; B1 = 0.6")

# sigma estimada vs. sigma observada
ts.plot(sim1[,"sigma2"], gas1$ft, col = c("#909090","red"), lwd = 1,
        main = "mean y = 1; w = 0.5; A0 = 0.2; B1 = 0.6")


# simulação 3 -------------------------------------------
# y[t] = mu[t-1] + sigma2[t-1]*epsilon[t], epsilon[t] ~ N(0,1)
# theta (parâmetros estáticos): c(w, A0, B1, w, A0, B1)

# simular série
sim1 <- simularGAS(density = "normal", n = 300, mean = "variavel", sd = "variavel", seed = 1, theta = c(1, 0.5, 0.2,
                                                                                                        0.5, 0.3, 0.5))
plot(sim1)

# estimar GAS
gas1 <- estimarGAS(y = sim1[,"y"], density = "normal", mean = "variavel", sd = "variavel", initial = c(1, 1, 1,
                                                                                                       1, 1, 1))
gas1$optim
round(gas1$optim,1)


# sigma estimada vs. sigma observada
ts.plot(sim1[,"sigma2"], gas1$ft[,"sigma2"], col = c("#909090","red"), lwd = 1,
        main = "mean y = 1; w = 0.5; A0 = 0.2; B1 = 0.6")

# média estimada vs. média observada
ts.plot(sim1[,"sigma2"], gas1$ft, col = c("#909090","red"), lwd = 1,
        main = "mean y = 1; w = 0.5; A0 = 0.2; B1 = 0.6")

# > RASCUNHO ----------------------------------------------

# # SIMULACAO 1
# 
# # número de observações
# n <- 1000
# 
# # vetor theta (parâmetros estáticos): c(sigma2, w, A0, B1)
# theta <- c(0.5, 0, 0.2, 0.4)
# 
# # epsilon 
# set.seed(1308)
# epsilon <- ts(rnorm(n, mean = 0, sd = sqrt(theta[1])))
# plot(epsilon)
# 
# # valor inicial de mu = 0
# mu <- 0
# y <- NA
# for(t in 2:n){
#   y[t] <- mu[t-1] + epsilon[t]
#   mu[t] <- theta[2] + theta[3]*(y[t] -  mu[t-1]) + theta[4]*mu[t-1]
# }
# 
# y <- ts(y)
# mu <- ts(mu)
# plot(y)
# plot(mu)
# # processo simulado
# 
# otimizar <- function(y, par){
#   mu <- 0 
#   for(t in 2:length(y)){
#     mu[t] <- par[2] + par[3]*(y[t] -  mu[t-1]) + par[4]*mu[t-1]
#   }
#   mu <- ts(mu)
#   loglik <- sum(-0.5*log(2*pi*par[1]) - 0.5*(y - lag(mu,-1))^2/par[1], na.rm = T)
#   -loglik
# }
# 
# 
# otimizados <- optim(par = c(2, 0.5, 0.5, 0.5), fn = otimizar, y = sim$y, method = "Nelder-Mead", hessian = F)
# otimizados
# param <- otimizados$par
# round(param,4)




# GAS package --------------------
library(GAS)
GASSpec <- UniGASSpec(Dist = "norm", ScalingType = "Inv",
                      GASPar = list(location = TRUE, scale = FALSE))
GASSpec
Fit <- UniGASFit(GASSpec,  k)
Fit
