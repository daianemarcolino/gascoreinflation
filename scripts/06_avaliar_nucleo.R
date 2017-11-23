nucleos <- readRDS("dados/nucleo.rds")

ipc <- ts(read.csv2("dados/VARIAÇÃO IPC-DI.csv")[,2], start = c(1999,1), freq = 12)
nucleo <- ((nucleo/100+1)^12-1)*100

ipc <- ((ipc/100+1)^12-1)*100
y=ipc
core=nucleos$core_res

y=window(ipc, start = c(2008,5), freq = 12)
core=window(nucleo, start = c(2008,5), freq = 12)

ts.plot(y,core)
write.csv2(data.frame(data = as.Date(core), round(((core/100+1)^12-1)*100,2)), "tendencia_final.csv", row.names = F)

