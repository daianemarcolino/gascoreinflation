library(shiny)
library(shinythemes)
library(shinyjs)
library(dygraphs)
library(TSA)
library(DT)
source("data/simularGAS.R")
source("data/dcs_fk_estimation.R")
graphs <- readRDS("data/dygraphs_nasdaq.rds")
resul <- readRDS("data/resultados_betategarch_estimation.rds")
resul_novo <- readRDS("data/novos_modelos.rds")


# turistas
bsm <- readRDS("data/bsm_turistas.rds")
turistas <- readRDS("data/turistas.rds")
smooth_turistas <- readRDS("data/smooth_turistas.rds")
diag_turistas <-  readRDS("data/diag_turistas.rds")

# ipc
ipc <- readRDS("data/ipc.rds")
ipc_ma1 <- window(readRDS("data/ipc_medias_aparadas.rds"), start = start(ipc), end = end(ipc), freq = 12)
ipc_ma2 <- window(readRDS("data/ipc_medias_aparadas_sa.rds"), start = start(ipc), end = end(ipc), freq = 12)
ipc_ma3 <- window(readRDS("data/ipc_medias_aparadas_sa_mm3.rds"), start = start(ipc), end = end(ipc), freq = 12)
todos_dcs <- readRDS("data/todos_dcs.rds")
todos_diags <- readRDS("data/todos_diags.rds")
tendencias <- readRDS("data/tendencias.rds")

# bsm2 <- readRDS("data/bsm_ipc.rds")
# smooth_ipc <- readRDS("data/smooth_ipc.rds")
# diag_ipc <-  readRDS("data/diag_ipc.rds")
# psd <- readRDS("data/smooth_ipc_psd.rds")
# diag_ipc <-  readRDS("data/diag_ipc.rds")

