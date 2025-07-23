library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)

# Ensure the "Results" folder exists (relative to the job's working directory)
if (!dir.exists("../Results")) {
  dir.create("../Results", recursive = TRUE)
}

# Load necessary data
SES <- read_csv("RegionalData.csv") |> 
  mutate(total_catch = ifelse(is.na(abundance) & total_catch == 0, NA, total_catch)) |>
  mutate(harvest = ifelse(is.na(abundance) & harvest == 0, NA, harvest)) |>
  mutate(landed = ifelse(is.na(abundance) & landed == 0, NA, landed)) |>
  mutate(ccatch = ifelse(total_catch == 0, 1, ifelse(total_catch == 1, 0.99, total_catch)),
         rcatch = ifelse(total_rec.catch == 0, 1, total_rec.catch)) |>
  mutate(log_ccatch = log(ccatch),
         log_rcatch = log(rcatch)) |> 
  mutate(spp.n = case_when(
    species == "Chinook"    ~ 1,
    species == "Chum"    ~ 2,
    species == "Coho"   ~ 3,
    species == "Pink"   ~ 4,
    species == "Sockeye" ~ 5,
    TRUE ~ NA_integer_  # For unmatched cases
  )) |> 
  mutate(region.n = case_when(
    region == "CC"    ~ 1,
    region == "Fraser"    ~ 2,
    region == "HG"   ~ 3,
    region == "Nass"   ~ 4,
    region == "Skeena" ~ 5,
    region == "VIMI" ~ 6,
    TRUE ~ NA_integer_  # For unmatched cases
  )) |> mutate(year.n = year - 1995) |>
  mutate(
    harvest.1000 = harvest*1000,
    charvest = ifelse(harvest.1000 == 0, 1, harvest.1000),
    log_harvest = log(charvest),
    wholesale.1000 = landed*1000,
    cwholesale = ifelse(wholesale.1000 == 0, 1, wholesale.1000),
    log_value = log(cwholesale)) |>
  mutate(recatch = ifelse(total_rec.catch == 0, 1, total_rec.catch),
         log_reccatch = log(recatch))

reclic <- read_csv("reclicences.csv") |> select(year, stamps)
SES <- SES |> left_join(reclic, by = c("year"))

gens <- read_csv("Genlengths.csv")
sst.gen <- read_csv("SST_genlengths.csv")

# List of species to iterate over
species_list <- c("Pink")

# List of Stan models to run
model_files <- c(
  "StanModels/CSE.stan",
  "StanModels/CSEF.stan",
  "StanModels/CSEFV.stan",
  #"StanModels/CSEFVH.stan",
  #"StanModels/CSEFLV.stan",
  "StanModels/CSEFRV.stan",
  "StanModels/CSEFRVStp.stan"
  #"StanModels/CSEFLRV.stan",
  #"StanModels/CSEFLRVStp.stan",
  #"StanModels/CSEFLVH.stan",
  #"StanModels/CSEFLRVStpH.stan"
)

# Import functions
source("RunFunctionsRaw.R")

# --- Main part of the script to run the models ---
num_chains_global <- 4
num_iterations_global <- 20000
warmup_length_global <- 10000

for (species in species_list) {
  data_for_species <- prepare_data(species, SES)
  
  for (model_file in model_files) {
    run_model_and_plot(
      species_name = species,
      model_file = model_file,
      data_prepared = data_for_species,
      num_chains = num_chains_global,
      num_iterations = num_iterations_global,
      warmup_length = warmup_length_global
    )
  }
}
# Plots are not coming out as they should for variables with missing data. Probably needs redoing after the fit 
# rather than fixing in the function

cat("--- All models and species processed. Results and plots saved in the 'Results' folder. ---\n")