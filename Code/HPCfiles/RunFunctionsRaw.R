# Function to prepare data for a given species
prepare_data <- function(species_name, all_data) {
  dat <- all_data |> filter(species == species_name)
  dat$temperature[is.na(dat$temperature)] <- 0 # Handle NA temperature
  
  list(
    N_total = nrow(dat),
    N_loc = length(unique(dat$region.n)),
    N_year_max = length(unique(dat$year.n)),
    
    log_Eobs = log(dat$total_effort),
    location_id = dat$region.n,
    year_index = dat$year.n,
    
    N_Cobs_total = sum(!is.na(dat$log_ccatch)),
    log_Cobs = dat$log_ccatch[!is.na(dat$log_ccatch)],
    Cobs_location_id = dat$region.n[!is.na(dat$log_ccatch)],
    Cobs_year_index = dat$year.n[!is.na(dat$log_ccatch)],
    
    N_Fobs_total = sum(!is.na(dat$total_fleet)),
    log_Fobs = log(dat$total_fleet[!is.na(dat$total_fleet)]),
    Fobs_location_id = dat$region.n[!is.na(dat$total_fleet)],
    Fobs_year_index = dat$year.n[!is.na(dat$total_fleet)],
    
    N_Sobs_total = sum(!is.na(dat$abundance)),
    log_Sobs = log(dat$abundance[!is.na(dat$abundance)]),
    Sobs_location_id = dat$region.n[!is.na(dat$abundance)],
    Sobs_year_index = dat$year.n[!is.na(dat$abundance)],
    
    # Hatcheries data
    N_Hobs_total = sum(!is.na(dat$release)),
    log_Hobs = log(dat$release[!is.na(dat$release)] / 1000),
    Hobs_location_id = dat$region.n[!is.na(dat$release)],
    Hobs_year_index = dat$year.n[!is.na(dat$release)],
    
    # Weight data
    N_Wobs_total = length(dat$log_harvest[!is.na(dat$log_harvest) & (dat$log_harvest != 0)]),
    log_Wobs = (dat$log_harvest[!is.na(dat$log_harvest) & (dat$log_harvest != 0)]),
    Wobs_location_id = dat$region.n[!is.na(dat$log_harvest) & (dat$log_harvest != 0)],
    Wobs_year_index = dat$year.n[!is.na(dat$log_harvest) & (dat$log_harvest != 0)],
    
    # Value data
    N_Vobs_total = length(dat$log_value[!is.na(dat$log_value) & (dat$log_value != 0)]),
    log_Vobs = (dat$log_value[!is.na(dat$log_value) & (dat$log_value != 0)]),
    Vobs_location_id = dat$region.n[!is.na(dat$log_value) & (dat$log_value != 0)],
    Vobs_year_index = dat$year.n[!is.na(dat$log_value) & (dat$log_value != 0)],
    
    # Recreational catch data
    N_Recobs = sum(!is.na(dat$log_reccatch)),
    log_RecCobs = (dat$log_reccatch[!is.na(dat$log_reccatch)]),
    log_RecEobs = log(dat$total_rec.effort[!is.na(dat$total_rec.effort)] / 100),
    Recobs_location_id = dat$region.n[!is.na(dat$log_reccatch)],
    Recobs_year_index = dat$year.n[!is.na(dat$log_reccatch)],
    
    #Salmon stamp data
    N_Lrecobs_total = length(unique(log(dat$stamps[!is.na(dat$stamps)]))),
    log_Lrecobs = unique(log(dat$stamps[!is.na(dat$stamps)])),
    Lrecobs_year_index = unique(dat$year.n[!is.na(dat$stamps)]),
    
    # License data
    N_Lobs_total = sum(!is.na(dat$licenses)),
    log_Lobs = log(dat$licenses[!is.na(dat$licenses)]),
    Lobs_location_id = dat$region.n[!is.na(dat$licenses)],
    Lobs_year_index = dat$year.n[!is.na(dat$licenses)],
    
    # Farms data
    N_farm_obs = sum(!is.na(dat$active.farms)),
    farms_obs = dat$active.farms[!is.na(dat$active.farms)],
    farm_obs_location_id = dat$region.n[!is.na(dat$active.farms)],
    farm_obs_year_index = dat$year.n[!is.na(dat$active.farms)],
    N_farm_missing = sum(is.na(dat$active.farms)),
    farm_missing_location_id = dat$region.n[is.na(dat$active.farms)],
    farm_missing_year_index = dat$year.n[is.na(dat$active.farms)],
    
    # Averages for centered regression
    mean_log_Cobs_loc = (tapply(dat$log_ccatch, dat$region.n, mean, na.rm = T)),
    mean_log_Wobs_loc = (tapply(dat$log_harvest, dat$region.n, mean, na.rm = T)),
    mean_log_Fobs_loc = (tapply(log(dat$total_fleet), dat$region.n, mean, na.rm = T)),
    mean_log_Lobs_loc = (tapply(log(dat$licenses), dat$region.n, mean, na.rm = T)),
    
    SST = dat$temperature[dat$region.n == 1],
    gentemp = sst.gen[[species_name]], ### Species-specific SST gen length
    genlength = gens[[species_name]]    ### Species-specific generation length
  )
}

## Function to run a Stan model, save results, and generate PPC plots
## Function to run a Stan model, save results, and generate PPC plots
run_model_and_plot <- function(species_name, model_file, data_prepared, num_chains, num_iterations, warmup_length) {
  cat(paste0("--- Running model: ", basename(model_file), " for ", species_name, " ---\n"))
  cat(paste0("--- Chains: ", num_chains, ", Iterations: ", num_iterations, ", Warmup: ", warmup_length, " ---\n"))
  
  # Compile the Stan model
  mod <- tryCatch(
    cmdstan_model(stan_file = model_file, pedantic = TRUE),
    error = function(e) {
      cat(paste0("Error compiling model ", basename(model_file), ": ", e$message, "\n"))
      return(NULL)
    }
  )
  
  if (is.null(mod)) {
    return(NULL)
  }
  
  # Define the base name for the output files, specific to the model and species
  output_base <- paste0(gsub("\\.stan$", "", basename(model_file)), "_", species_name)
  output_dir <- "../Results/" # Consistent output directory
  
  # Sample from the posterior
  Samp <- tryCatch(
    mod$sample(
      data = data_prepared,
      chains = num_chains,
      parallel_chains = num_chains,
      iter_warmup = warmup_length,
      iter_sampling = num_iterations - warmup_length,
      adapt_delta = 0.9,
      seed = 123 # For reproducibility
    ),
    error = function(e) {
      cat(paste0("Error sampling from model ", basename(model_file), " for ", species_name, ": ", e$message, "\n"))
      return(NULL)
    }
  )
  
  if (is.null(Samp)) {
    return(NULL)
  }
  
  # Save the output files and the object, with error handling
  # tryCatch({
  #   Samp$save_output_files(dir = output_dir, basename = output_base)
  #   cat(paste0("Sampling results saved to: ", file.path(output_dir, paste0(output_base, "-chain-*.csv")), "\n"))
  # }, error = function(e) {
  #   cat(paste0("Error saving output files for ", species_name, ": ", e$message, "\n"))
  # })
  
  tryCatch({
    Samp$save_object(file = file.path(output_dir, paste0(output_base, ".RDS")))
    cat(paste0("Sampling object saved to: ", file.path(output_dir, paste0(output_base, ".RDS")), "\n"))
  }, error = function(e) {
    cat(paste0("Error saving sampling object for ", species_name, ": ", e$message, "\n"))
  })
  
  tryCatch({
    write_csv(Samp$summary(), file = file.path(output_dir, paste0(output_base, "_summary.csv")))
    cat(paste0("Sampling summary saved to: ", file.path(output_dir, paste0(output_base, "_summary.csv")), "\n"))
  }, error = function(e) {
    cat(paste0("Error saving sampling summary for ", species_name, ": ", e$message, "\n"))
  })
  
  return(invisible(NULL)) # Function completes here
}