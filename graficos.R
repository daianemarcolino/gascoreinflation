# gráficos dos últimos resultados

smooth <- ts(readRDS("ultimos_resultados.rds")[,-1], start = c(1999,1), freq = 12)
head(smooth)

dygraph(smooth[,c("ipc","mu")]) %>%
  dySeries("ipc", strokePattern = "dotted", color = "black") %>%
  dySeries("mu", strokeWidth = 2, color = "orangered")

dygraph(smooth[,c("mu")]) %>%
  dySeries("V1", strokeWidth = 2, color = "orangered")

dygraph(smooth[,c("ipc","gamma")]) %>%
  dySeries("ipc", strokePattern = "dotted", color = "black") %>%
  dySeries("gamma", strokeWidth = 2, color = "dodgerblue")

dygraph(smooth[,c("gamma")]) %>%
  dySeries("V1", strokeWidth = 2, color = "dodgerblue")

# SMOOTHER
dygraph(smooth[,c("ipc","mu_smooth")]) %>%
  dySeries("ipc", strokePattern = "dotted", color = "black") %>%
  dySeries("mu_smooth", strokeWidth = 2, color = "orangered")

dygraph(smooth[,c("mu_smooth")]) %>%
  dySeries("V1", strokeWidth = 2, color = "orangered")

dygraph(smooth[,c("ipc","gamma_smooth")]) %>%
  dySeries("ipc", strokePattern = "dotted", color = "black") %>%
  dySeries("gamma_smooth", strokeWidth = 2, color = "dodgerblue")

dygraph(smooth[,c("gamma_smooth")]) %>%
  dySeries("V1", strokeWidth = 2, color = "dodgerblue")

dygraph(smooth[,c("epsilon","score")]) %>%
  dySeries("epsilon", strokeWidth = 1, color = "#F64D54") %>%
  dySeries("score", strokeWidth = 2, color = "#91219E")

hist(smooth[,"epsilon"], border = "#F64D54", col = "#FFC0CB", breaks = seq(min(smooth[,"epsilon"]), max(smooth[,"epsilon"]), length.out = 12), xlim = c(-5,13), main = "epsilon")
hist(smooth[,"score"], border = "#91219E", col = "#DCA2CD", breaks = seq(min(smooth[,"score"]), max(smooth[,"score"]), length.out = 12), xlim = c(-0.4,0.4), main = "score")
