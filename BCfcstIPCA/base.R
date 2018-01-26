library(BETS)
library(seasonal)
source("functions/ifa_functions.R")

ipca_livres <- BETS.get("11428")
ipca_adms <- BETS.get("4449")
cambio_usd <- BETS.get("3695")
selic <- BETS.get("4189")
igpm <- BETS.get("189")
ipca <- BETS.get("433")
ipca12 <- acum12(ipca)


juros_reais_igpm <- ((selic/100 + 1)/(igpm/100 + 1) - 1)*100
juros_reais_ipca <- ((selic/100 + 1)/(ipca12/100 + 1) - 1)*100

m4 <-  BETS.get("1843")
b <-  BETS.get("1788")
m4b <- m4/b

pim_sajuste <- ts(BETS.sidra.get(x = 3653, from = c("200201"), to = paste0(format(Sys.Date(),"%Y"),format(Sys.Date(),"%m")), variable = 3135, sections = 129314, cl = 544)[[1]]$Valor,
                  start = c(2002,01), freq = 12)

pim_cajuste <- ts(BETS.sidra.get(x = 3653, from = c("200201"), to = paste0(format(Sys.Date(),"%Y"),format(Sys.Date(),"%m")), variable = 3134, sections = 129314, cl = 544)[[1]]$Valor,
                  start = c(2002,01), freq = 12)

m1 <- BETS.get("1827")

# base modelos VAR

base_var1 <- window(cbind(ipca_livres, ipca_adms, cambio_usd, juros_reais_igpm), start = c(2001,1), end = c(2017,12), freq = 12)
base_var2 <- window(cbind(ipca_livres, m4b, cambio_usd, juros_reais_ipca), start = c(2001,1), end = c(2017,12), freq = 12)
base_var3 <- window(cbind(ipca_livres, pim = pim_cajuste, juros_reais_igpm), start = c(2001,1), end = c(2017,12), freq = 12)
base_var4 <- window(cbind(ipca_livres, m1, juros_reais_igpm), start = c(2001,1), end = c(2017,12), freq = 12)

saveRDS(list(base_var1 = base_var1, base_var2 = base_var2, base_var3 = base_var3, base_var4 = base_var4), "BCfcstIPCA/base_var.rds")

# base modelo estrutural de pequeno porte

