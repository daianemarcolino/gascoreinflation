library(shiny)
library(shinythemes)
library(shinyjs)
library(dygraphs)
library(TSA)
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
bsm2 <- readRDS("data/bsm_ipc.rds")
ipc <- readRDS("data/ipc.rds")
smooth_ipc <- readRDS("data/smooth_ipc.rds")
diag_ipc <-  readRDS("data/diag_ipc.rds")
psd <- readRDS("data/smooth_ipc_psd.rds")
