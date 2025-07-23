setwd("~/Documents/McGill/PhD/Chapter 3/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(shinystan)
library(rstan)
library(posterior)
library(loo)
library(bayesplot)

source("LoadDataandFunctions.R")

#Predictions for Chinook
ChinookList <- getdatalist("Chinook")
ChinookPreds <- readRDS("../Results/HPCfits/CSE4chain/CSEFRV_Chinook.RDS")

#Predictions for Pink and Sockeye
# PinkList <- getdatalist("Pink")
# PinkPreds <- readRDS("../Results/HPCfits/Third8ChainRun/CatchSupplyEffortValue2T_Pink.RDS")
SockeyeList <- getdatalist("Sockeye")
SockeyePreds <- readRDS("../Results/HPCfits/CSE4chain/CSEF_Sockeye.RDS") ## also CSEFRV and CSEFRVStp & CSEFVH

#Predictions for Chum and Coho
ChumList <- getdatalist("Chum")
ChumPreds <- readRDS("../Results/HPCfits/CSE4chain/CSEFV_Chum.RDS")
CohoList <- getdatalist("Coho")
CohoPreds <- readRDS("../Results/HPCfits/CSE4chain/CSEFRV_Coho.RDS")

# #Load predictions only Stan files
# ValueRec_pred <- cmdstan_model(stan_file = "StanModels/Preds_CatchSupplyEffortValueRec.stan", pedantic=T)
# Value_pred <- cmdstan_model(stan_file = "StanModels/Preds_CatchSupplyEffortValue.stan", pedantic=T)
# 
# #Run predictions
# ChinookPreds <- ValueRec_pred$generate_quantities(ChinookList, fitted_params = ChinookMod)
# CohoPreds <- ValueRec_pred$generate_quantities(CohoList, fitted_params = CohoMod)
# ChumPreds <- ValueRec_pred$generate_quantities(ChumList, fitted_params = ChumMod)
# PinkPreds <- ValueRec_pred$generate_quantities(PinkList, fitted_params = PinkMod)
# SockeyePreds <- ValueRec_pred$generate_quantities(SockeyeList, fitted_params = SockeyeMod)
# 
# ChinookPreds <- Value_pred$generate_quantities(ChinookList, fitted_params = ChinookMod)
# CohoPreds <- Value_pred$generate_quantities(CohoList, fitted_params = CohoMod)
# ChumPreds <- Value_pred$generate_quantities(ChumList, fitted_params = ChumMod)
# PinkPreds <- Value_pred$generate_quantities(PinkList, fitted_params = PinkMod)
# SockeyePreds <- Value_pred$generate_quantities(SockeyeList, fitted_params = SockeyeMod)


#Plot model outputs - set T/F depending on presence of rec catch and value
plotTSpreds(ChumPreds,"Chum",F,T)
plotTSpreds(ChinookPreds,"Chinook",T,T)
plotTSpreds(CohoPreds,"Coho",T,T)
plotTSpreds(SockeyePreds,"Sockeye",F,F)

#Plot model outputs on the same scale - set T/F depending on presence of rec catch and value
plotTSpreds2(ChumPreds,"Chum",F,T)
plotTSpreds2(ChinookPreds,"Chinook",T,T)
plotTSpreds2(CohoPreds,"Coho",T,T)
plotTSpreds2(SockeyePreds,"Sockeye",F,F)

#Plot harvest rates time series
plothrate2(ChinookPreds,"Chinook")
plothrate(ChumPreds,"Chum")
plothrate2(CohoPreds,"Coho")
plothrate(SockeyePreds,"Sockeye")

#Chinook param plots
plot_param(ChinookPreds,"Chinook", "alpha", 3, "blue")
plot_param(ChinookPreds,"Chinook","beta", 3, "green")
plot_param(ChinookPreds,"Chinook","log_q", NULL, "red")
plot_param(ChinookPreds,"Chinook","log_qrec", NULL, "pink")
plot_param(ChinookPreds,"Chinook","Fslope", NULL, "yellow")
plot_param(ChinookPreds,"Chinook","Fexp", NULL, "orange")
plot_param(ChinookPreds,"Chinook","Cslope", NULL, "brightblue")
plot_param(ChinookPreds,"Chinook","Cexp", NULL, "teal")

#Chum param plots
plot_param(ChumPreds,"Chum", "alpha", NULL, "blue")
plot_param(ChumPreds,"Chum","beta", NULL, "green")
plot_param(ChumPreds,"Chum","log_q", NULL, "red")
plot_param(ChumPreds,"Chum","Fslope", NULL, "yellow")
plot_param(ChumPreds,"Chum","Fexp", NULL, "orange")
plot_param(ChumPreds,"Chum","Cslope", NULL, "brightblue")
plot_param(ChumPreds,"Chum","Cexp", NULL, "teal")

#Coho param plots
plot_param(CohoPreds,"Coho", "alpha", NULL, "blue")
plot_param(CohoPreds,"Coho","beta", NULL, "green")
plot_param(CohoPreds,"Coho","log_q", NULL, "red")
plot_param(CohoPreds,"Coho","log_qrec", NULL, "pink")
plot_param(CohoPreds,"Coho","Fslope", NULL, "yellow")
plot_param(CohoPreds,"Coho","Fexp", NULL, "orange")
plot_param(CohoPreds,"Coho","Cslope", NULL, "brightblue")
plot_param(CohoPreds,"Coho","Cexp", NULL, "teal")

#Sockeye param plots
plot_param(SockeyePreds,"Sockeye", "alpha", NULL, "blue")
plot_param(SockeyePreds,"Sockeye","beta", NULL, "green")
plot_param(SockeyePreds,"Sockeye","log_q", NULL, "red")
plot_param(SockeyePreds,"Sockeye","Fslope", NULL, "yellow")
plot_param(SockeyePreds,"Sockeye","Fexp", NULL, "orange")
