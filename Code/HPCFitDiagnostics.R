setwd("~/Documents/McGill/PhD/Chapter 3/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(shinystan)
library(rstan)
library(posterior)
library(loo)

source("LoadDataandFunctions.R")

#List the models to be used
modelnames <- c(
  "CSE",
  "CSEF", 
  "CSEFV", 
  "CSEFRV",
  "CSEFRVStp"
  )

# List of species to iterate over
species_list <- c("Chinook","Chum","Coho","Pink","Sockeye")

#Produce all PPC plots and save all R2 values for all models, variables and species
#also save summaries
#hash out any section to reduce run time
for (spp in species_list) {
  print(paste("Starting species", spp))
  
  print("Loading fits")
  fit1 <- readRDS(paste0("../Results/HPCfits/CSE4chain/CSE_",spp,".RDS"))
  fit2 <- readRDS(paste0("../Results/HPCfits/CSE4chain/CSEF_",spp,".RDS"))
  fit3 <- readRDS(paste0("../Results/HPCfits/CSE4chain/CSEFV_",spp,".RDS"))
  fit4 <- readRDS(paste0("../Results/HPCfits/CSE4chain/CSEFRV_",spp,".RDS"))
  fit5 <- readRDS(paste0("../Results/HPCfits/CSE4chain/CSEFRVStp_",spp,".RDS"))

  getR2s(fit1,spp,"CSE_")
  getR2s(fit2,spp,"CSEF_")
  getR2s(fit3,spp,"CSEFV_")
  getR2s(fit4,spp,"CSEFRV_")
  getR2s(fit5,spp,"CSEFRVStp_")
  
  #set region to F to avoid plotting regional PPCs
  modellist <- list(fit1,fit2,fit3,fit4,fit5)
  SaveAllPlots(modellist,spp,modelnames,region=T)
}

#Conduct LOO validation
spp <- "Coho" #change species here
modelfit <- readRDS(paste0("../Results/HPCfits/CSE4chain/CSEFV_",spp,".RDS")) #change model here
computeLOO(spp,modelfit,"log_lik_spawners","log_Sobs")
computeLOO(spp,modelfit,"log_lik_catch","log_Cobs")
computeLOO(spp,modelfit,"log_lik_effort","log_Eobs")
computeLOO(spp,modelfit,"log_lik_fleet","log_Fobs")
computeLOO(spp,modelfit,"log_lik_weight","log_Wobs")
computeLOO(spp,modelfit,"log_lik_value","log_Vobs")
computeLOO2(spp,modelfit,"log_lik_reccatch","log_RecCobs") #recheck this for computeLOO normal
computeLOO2(spp,modelfit,"log_lik_receffort","log_RecEobs") #recheck this for computeLOO normal
computeLOO(spp,modelfit,"log_lik_hatcheries","log_Hobs")
computeLOO(spp,modelfit,"log_lik_license","log_Lobs")
computeLOO(spp,modelfit,"log_lik_stamps","log_Lrecobs")


#Get R2 summaries
summarize_R2("../Results/Diagnostics/","Chinook/")
summarize_R2("../Results/Diagnostics/","Coho/")
summarize_R2("../Results/Diagnostics/","Chum/")
summarize_R2("../Results/Diagnostics/","Pink/")
summarize_R2("../Results/Diagnostics/","Sockeye/")

#Get fit summaries
summarize_fit("../Results/Diagnostics/","Chinook/")
summarize_fit("../Results/Diagnostics/","Coho/")
summarize_fit("../Results/Diagnostics/","Chum/")
summarize_fit("../Results/Diagnostics/","Pink/")
summarize_fit("../Results/Diagnostics/","Sockeye/")
