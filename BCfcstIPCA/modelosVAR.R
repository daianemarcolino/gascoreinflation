
nucleo_tf <- readRDS("data/nucleo_tf.rds")[,3]
nucleo_dcs <- window(readRDS("data/nucleo_dcs.rds"), end = c(2017,12), freq = 12)
base_var <- readRDS("BCfcstIPCA/base_var.rds")

base_var1 <- window(base_var$base_var1, start = c(2002,1), end = c(2017,11), freq = 12)
base_var2 <- window(base_var$base_var2, start = c(2002,1), end = c(2017,11), freq = 12)
base_var3 <- window(base_var$base_var3, start = c(2002,1), end = c(2017,11), freq = 12)
base_var4 <- window(base_var$base_var4, start = c(2002,1), end = c(2017,11), freq = 12)

base_var2[,"m4b"] <- final(seas(base_var2[,"m4b"]))

plot(base_var1)
base_var1[,3] <- c(NA, diff(base_var1[,3]))
base_var1[,4] <- c(NA, diff(base_var1[,4]))
plot(base_var2)
base_var2[,2] <- c(NA, diff(base_var2[,2]))
base_var2[,3] <- c(NA, diff(base_var2[,3]))
base_var2[,4] <- c(NA, diff(base_var2[,4]))
plot(base_var3)
base_var3[,2] <- c(NA, diff(base_var3[,2]))
base_var3[,3] <- c(NA, diff(base_var3[,3]))
plot(base_var4)
base_var4[,2] <- c(NA, diff(base_var4[,2]))
base_var4[,3] <- c(NA, diff(base_var4[,3]))

base_var1 <- window(base_var1, start = c(2005,1), end = c(2017,11), freq = 12)
base_var2 <- window(base_var2, start = c(2005,1), end = c(2017,11), freq = 12)
base_var3 <- window(base_var3, start = c(2005,1), end = c(2017,11), freq = 12)
base_var4 <- window(base_var4, start = c(2005,1), end = c(2017,11), freq = 12)

plot(base_var4)
qs(seas(base_var2[,4]))


# modelo VAR 1
VARselect(base_var1)
VARselect(base_var2)
VARselect(base_var3)
VARselect(base_var4)

VAR1 <- VAR(base_var1, p = 2)
summary(VAR1)

# modelo VAR 2
VAR2 <- VAR(base_var2, p = 1)
summary(VAR2)

# modelo VAR 3
VAR3 <- VAR(base_var3, p = 3)
summary(VAR3)

# modelo VAR 4
VAR4 <- VAR(base_var4, p = 1)
summary(VAR4)


h = 24
M <- matrix(NA, ncol = 4, nrow = h)
for(i in h:1){
  var1 <- VAR(base_var1[1:(nrow(base_var1)-i),], p = 2)
  var1_fcst <- predict(var1, n.ahead = 1, ci = 0.95)
  var2 <- VAR(base_var2[1:(nrow(base_var2)-i),], p = 1)
  var2_fcst <- predict(var2, n.ahead = 1, ci = 0.95)
  var3 <- VAR(base_var3[1:(nrow(base_var3)-i),], p = 3)
  var3_fcst <- predict(var3, n.ahead = 1, ci = 0.95)
  var4 <- VAR(base_var4[1:(nrow(base_var4)-i),], p = 1)
  var4_fcst <- predict(var4, n.ahead = 1, ci = 0.95)
  
  M[h-i+1,1] <- round(var1_fcst$fcst$ipca_livres[,1],4)
  M[h-i+1,2] <- round(var2_fcst$fcst$ipca_livres[,1],4)
  M[h-i+1,3] <- round(var3_fcst$fcst$ipca_livres[,1],4)
  M[h-i+1,4] <- round(var4_fcst$fcst$ipca_livres[,1],4)
}
M <- ts(M, end = end(base_var1), freq = 12)

sqrt(colMeans(M - ipca_livres)^2)

ts.plot(M)

ts.plot(window(base_var1[,"ipca_livres"], start = c(2015,12), freq = 12), M, lty = c(3,1,1,1,1), col = c(2,1,1,1,1),
        lwd = c(1,2,3,4,5))

ts.plot(window(base_var1[,"ipca_livres"], start = c(2015,12), freq = 12), 
        M[,1], lty = c(3,1,1,1,1), col = c(2,1,1,1,1), lwd = c(1,2,3,4,5))
