# script desenvolvido por Daiane Marcolino de Mattos
# objetivo: simular GAS
# início da programação: 13/08/2017 às 18h26 - domingo

# carregar funções --------------------------------------
source("functions/simularGAS.R")

# simulação 1 -------------------------------------------
# y[t] = mu[t-1] + epsilon[t], epsilon[t] ~ N(0,sigma2)
# theta (parâmetros estáticos): c(sigma2, w, A0, B1)

# simular série
sim1 <- simularGAS(density = "normal", n = 300, link = F, mean = "variavel", sd = "fixa", seed = 1, theta = c(0.5, 1, 0.5, 0.8))
plot(sim1)


# simulação 2 -------------------------------------------
# y[t] = mu + sigma2[t-1]*epsilon[t], epsilon[t] ~ N(0,1)
# theta (parâmetros estáticos): c(mean y, w, A0, B1)

# simular série
sim2 <- simularGAS(density = "normal", n = 300, link = F, mean = "fixa", sd = "variavel", 
                   seed = 1, theta = c(-1, 2, 0.5, 0.7))
plot(sim2)

# simular série
sim2 <- simularGAS(density = "normal", n = 300, link = T, mean = "fixa", sd = "variavel", 
                   seed = 1, theta = c(-1, 2, 0.5, 0.7))
plot(sim2)


# simulação 3 -------------------------------------------
# y[t] = mu[t-1] + sigma2[t-1]*epsilon[t], epsilon[t] ~ N(0,1)
# theta (parâmetros estáticos): c(w, A0, B1, w, A0, B1)

# simular série
sim1 <- simularGAS(density = "normal", n = 300, link = F, mean = "variavel", sd = "variavel", 
                   seed = 1, theta = c(1, 0.5, 0.2, 0.5, 0.3, 0.5))
plot(sim1)

# simular série
sim1 <- simularGAS(density = "normal", n = 300, link = T, mean = "variavel", sd = "variavel", 
                   seed = 1, theta = c(1, 0.5, 0.2, 0.5, 0.3, 0.5))
plot(sim1)


# GAS package --------------------
library(GAS)
GASSpec <- UniGASSpec(Dist = "norm", ScalingType = "Inv",
                      GASPar = list(location = TRUE, scale = FALSE))
GASSpec
Fit <- UniGASFit(GASSpec,  k)
Fit
