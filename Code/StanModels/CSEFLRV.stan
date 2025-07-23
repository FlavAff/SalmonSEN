//
// This Stan program defines a partially hierarchical model for salmon abundance and catch
// across multiple locations.
//

data {
  int<lower=1> N_total;                // Total number of observations across all locations
  int<lower=1> N_loc;                  // Number of locations
  int<lower=1> N_year_max;             // Maximum number of years across all locations

  vector[N_total] log_Eobs;             // Observed effort
  array[N_total] int<lower=1, upper=N_loc> location_id; // Identifier for each observation's location
  array[N_total] int<lower=1, upper=N_year_max> year_index; // Year index for each observation

  // Catch data
  int<lower=0> N_Cobs_total;           // Total number of observed catch values
  vector[N_Cobs_total] log_Cobs;     // Observed log catch (flattened, only where available)
  array[N_Cobs_total] int<lower=1, upper=N_loc> Cobs_location_id; // Location ID for each catch observation
  array[N_Cobs_total] int<lower=1, upper=N_year_max> Cobs_year_index; // Year index for each catch observation

  // Fleet size data
  int<lower=0> N_Fobs_total;           // Total number of observed fleet size values
  vector[N_Fobs_total] log_Fobs;     // Observed log fleet size (flattened, only where available)
  array[N_Fobs_total] int<lower=1, upper=N_loc> Fobs_location_id; // Location ID for each fleet size observation
  array[N_Fobs_total] int<lower=1, upper=N_year_max> Fobs_year_index; // Year index for each fleet size observation

  // Spawners data
  int<lower=0> N_Sobs_total;           // Total number of observed spawner values
  vector[N_Sobs_total] log_Sobs;       // Observed spawner values (flattened)
  array[N_Sobs_total] int<lower=1, upper=N_loc> Sobs_location_id; // Location ID for each spawner observation
  array[N_Sobs_total] int<lower=1, upper=N_year_max> Sobs_year_index; // Year index for each spawner observation
  
  // Commercial license data
  int<lower=0> N_Lobs_total;            // Total number of observed license counts
  vector[N_Lobs_total] log_Lobs; // Observed log license counts (flattened, only where available)
  array[N_Lobs_total] int<lower=1, upper=N_loc> Lobs_location_id; // Location ID for each license observation
  array[N_Lobs_total] int<lower=1, upper=N_year_max> Lobs_year_index; // Year index for each license observation

  // Recreational catch data
  int<lower=0> N_Recobs;            // Total number of observed recreational catch & effort
  vector[N_Recobs] log_RecCobs;       // Observed recreational catch values (flattened)
  vector[N_Recobs] log_RecEobs;       // Observed recreational catch values (flattened)
  array[N_Recobs] int<lower=1, upper=N_loc> Recobs_location_id; // Location ID for each recreational observation
  array[N_Recobs] int<lower=1, upper=N_year_max> Recobs_year_index; // Year index for each recreational observation

  // Value data
  int<lower=0> N_Vobs_total;           // Total number of observed wholesale values
  vector[N_Vobs_total] log_Vobs;       // Observed wholesale values (flattened)
  array[N_Vobs_total] int<lower=1, upper=N_loc> Vobs_location_id; // Location ID for each wholesale value observation
  array[N_Vobs_total] int<lower=1, upper=N_year_max> Vobs_year_index; // Year index for each wholesale value observation

  // Genlength data
  array[N_loc] int<lower=1> genlength;      // Generation length for each location
  
  // Means for centering
  vector[N_loc] mean_log_Lobs_loc;
}

parameters {
  // Supply hyperparameters
  // Hyperparameters for alpha
  real mu_alpha;
  real<lower=0> sigma_alpha;
  // Hyperparameters for beta
  real mu_beta;
  real<lower=0> sigma_beta;
  // Catch hyperparameters
  // Hyperparameters for log_q
  real mu_log_q;
  real<lower=0> sigma_log_q;
  // Hyperparameters for theta (logit scale)
  real mu_theta_logit;
  real<lower=0> sigma_theta_logit;
  // Hyperparameter for exponential effect of fleet
  real mu_Fexp;
  real<lower=0> sigma_Fexp;
  // Hyperparameter for slope effect of fleet
  real mu_Fslope;
  real<lower=0> sigma_Fslope;
  // Weight & value hyperparameters
  // Hyperparameter for exponential effect of catch
  real mu_Cexp;
  real<lower=0> sigma_Cexp;
  // Hyperparameter for slope effect of catch
  real mu_Cslope;
  real<lower=0> sigma_Cslope;
  // Recreational catch hyperparameters
  // Hyperparameters for log_qrec 
  real mu_log_qrec;
  real<lower=0> sigma_log_qrec;
  // Hyperparameters for thetarec (logit scale)
  real mu_thetarec_logit;
  real<lower=0> sigma_thetarec_logit;
  // Licenses hyperparameters
  // Fleet size Hyperparameters for License Activation Rate Model
  real mu_activation_intercept;
  real<lower=0> sigma_activation_intercept;
  real mu_activation_slope;
  real<lower=0> sigma_activation_slope;
  
  // Missing data hyperparameters where no regressors are available
  // Hyperparameters for log_L
  real mu_log_L;
  real<lower=0> sigma_log_L;
  // Hyperparameters for log_Erec
  real mu_log_Erec;
  real<lower=0> sigma_log_Erec;

  // Location-specific parameters (drawn from the population-level distribution)
  vector[N_loc] alpha_raw;
  vector[N_loc] beta_raw;
  vector[N_loc] log_q_raw;
  vector[N_loc] theta_logit_raw;
  vector[N_loc] Fexp_raw;
  vector[N_loc] Fslope_raw;
  vector[N_loc] Cexp_raw;
  vector[N_loc] Cslope_raw;
  vector[N_loc] log_qrec_raw;
  vector[N_loc] thetarec_logit_raw;
  vector[N_loc] activation_intercept_raw;
  vector[N_loc] activation_slope_raw;

  // Parameters for latent objects without regressors
  array[N_loc] vector[N_year_max] log_L;
  array[N_loc] vector[N_year_max] log_Erec;
  
  // Location-specific parameters (not hierarchical)
  vector<lower=0>[N_loc] sigma_r;      // Recruitment process error for each location
  vector<lower=0>[N_loc] sigma_s;      // Spawner observation error for each location
  vector<lower=0>[N_loc] sigma_c;      // Catch observation error for each location
  vector<lower=0>[N_loc] sigma_e;      // Effort observation error for each location
  vector<lower=0>[N_loc] sigma_f;              // Observation error for fleet
  vector<lower=0>[N_loc] sigma_l;              // Observation error for licenses
  vector<lower=0>[N_loc] sigma_crec;       // Recreational catch observation error for each location
  vector<lower=0>[N_loc] sigma_erec;       // Recreational effort observation error for each location
  vector<lower=0>[N_loc] sigma_v;             // Observation error for value
  array[N_loc] vector[N_year_max] log_R_Perror; // Latent recruitment with error for each location
  array[N_loc] vector[N_year_max] start_log_R; // Starting years for each location
  vector<lower=0>[N_loc] df_c;  // Degrees of freedom for Student's t-distribution
}


transformed parameters {
  vector[N_loc] alpha;
  vector<lower=0>[N_loc] beta;
  vector[N_loc] log_q;
  vector<lower=0, upper=1>[N_loc] theta; // Hurdle parameter for each location
  vector<lower=0>[N_loc] Fexp;
  vector[N_loc] Fslope;
  vector<lower=0>[N_loc] Cexp;
  vector[N_loc] Cslope;
  
  array[N_loc] vector[N_year_max] logit_h;
  array[N_loc] vector[N_year_max] log_R;
  array[N_loc] vector[N_year_max] log_S;
  array[N_loc] vector[N_year_max] log_C;
  array[N_loc] vector[N_year_max] log_E;
  array[N_loc] vector[N_year_max] log_F;
  array[N_loc] vector[N_year_max] log_V;
  array[N_loc] vector[N_year_max] log_Crec;
  vector[N_loc] log_qrec;
  vector<lower=0, upper=1>[N_loc] thetarec; 
  array[N_loc] vector[N_year_max] logit_hrec;
  
  vector[N_loc] activation_intercept;
  vector[N_loc] activation_slope;
  array[N_loc] vector[N_year_max] logit_p_activation;
  
  // Estimate the real parameters by decomposing variance from site and overall averages using a centering and scaling approach
  alpha = mu_alpha + sigma_alpha * alpha_raw;
  beta = exp(mu_beta + sigma_beta * beta_raw);
  log_q = mu_log_q + sigma_log_q * log_q_raw;
  theta = inv_logit(mu_theta_logit + sigma_theta_logit * theta_logit_raw);
  Fexp = exp(mu_Fexp + sigma_Fexp * Fexp_raw);
  Fslope = mu_Fslope + sigma_Fslope * Fslope_raw;
  Cexp = exp(mu_Cexp + sigma_Cexp * Cexp_raw);
  Cslope = mu_Cslope + sigma_Cslope * Cslope_raw;
  log_qrec = mu_log_qrec + sigma_log_qrec * log_qrec_raw;
  thetarec = inv_logit(mu_thetarec_logit + sigma_thetarec_logit * thetarec_logit_raw);
  activation_intercept = mu_activation_intercept + sigma_activation_intercept * activation_intercept_raw;
  activation_slope = mu_activation_slope + sigma_activation_slope * activation_slope_raw;
  
  //Catch and abundance models
  for (l in 1:N_loc) {
    logit_h[l] = rep_vector(0.0, N_year_max);
    log_E[l] = rep_vector(0.0, N_year_max);
    log_R[l] = rep_vector(0.0, N_year_max);
    log_S[l] = rep_vector(0.0, N_year_max);
    log_C[l] = rep_vector(0.0, N_year_max);
    log_V[l] = rep_vector(0.0, N_year_max);
    log_F[l] = rep_vector(0.0, N_year_max);
    logit_hrec[l] = rep_vector(0.0, N_year_max);
    log_Crec[l] = rep_vector(0.0, N_year_max);
    logit_p_activation[l] = rep_vector(0.0, N_year_max);
    
    // Calculate logit_h based on effort
    for (y in 1:N_year_max) {
      logit_p_activation[l][y] = activation_intercept[l] + activation_slope[l] * (log_L[l][y] - mean_log_Lobs_loc[l]);
      log_F[l][y] = log_inv_logit(logit_p_activation[l][y]) + log_L[l][y];
      log_E[l][y] = Fslope[l] + Fexp[l] * log_F[l][y];
      logit_h[l][y] = log_E[l][y] + log_q[l]; // log_E[l][y] is now a parameter
      logit_hrec[l][y] = log_Erec[l][y] + log_qrec[l];
    }
    
    // Initialize the first 'genlength[l]' values of log_R for location 'l'
    for (i in 1:genlength[l]) {
      vector[3] unnormalized_log_probs_init;
      vector[3] log_probs_vector_init;
      
      // Define log-odds for commercial catch, recreational catch, and spawners
      unnormalized_log_probs_init[1] = logit_h[l][i];    // Logit for commercial catch proportion
      unnormalized_log_probs_init[2] = logit_hrec[l][i]; // Logit for recreational catch proportion
      unnormalized_log_probs_init[3] = 0;                // Logit for spawner proportion (reference category)
      
      // Transform to log probabilities that sum to log(1) = 0 (i.e., probabilities sum to 1)
      log_probs_vector_init = log_softmax(unnormalized_log_probs_init);
      
      // Calculate log catch and log spawners
      log_R[l][i] = start_log_R[l][i];
      log_C[l][i]    = log_R[l][i] + log_probs_vector_init[1];
      log_Crec[l][i] = log_R[l][i] + log_probs_vector_init[2];
      log_S[l][i]    = log_R[l][i] + log_probs_vector_init[3];
    }
    
    // Calculate the remaining values of log_R, log_C, and log_S for location 'l'
    for (t in (genlength[l] + 1):N_year_max) {
      vector[3] unnormalized_log_probs_dyn;
      vector[3] log_probs_vector_dyn;
      
      // Calculate recruits based on previous year's spawners
      // This Ricker formulation is standard.
      log_R[l][t] = log_S[l][t - genlength[l]] + alpha[l] - beta[l] * exp(log_S[l][t - genlength[l]]) + log_R_Perror[l][t];
      
      // Define log-odds for commercial catch, recreational catch, and spawners
      unnormalized_log_probs_dyn[1] = logit_h[l][t];    // Logit for commercial catch proportion
      unnormalized_log_probs_dyn[2] = logit_hrec[l][t]; // Logit for recreational catch proportion
      unnormalized_log_probs_dyn[3] = 0;                // Logit for spawner proportion (reference category)
      
      // Transform to log probabilities
      log_probs_vector_dyn = log_softmax(unnormalized_log_probs_dyn);
      
      // Calculate log catch and log spawners
      log_C[l][t]    = log_R[l][t] + log_probs_vector_dyn[1];
      log_Crec[l][t] = log_R[l][t] + log_probs_vector_dyn[2];
      log_S[l][t]    = log_R[l][t] + log_probs_vector_dyn[3];
    }
    
    // Now calculate harvest and value from populated other parameters
    for (y in 1:N_year_max) {
      log_V[l][y] = Cslope[l] + Cexp[l] * log_C[l][y];
    }
  }
}

model {
  // Hyperpriors for hierarchical parameters
  mu_alpha ~ normal(1, 1);
  sigma_alpha ~ exponential(1);
  mu_beta ~ normal(0, 1);
  sigma_beta ~ exponential(1);
  mu_log_q ~ normal(0, 1);
  sigma_log_q ~ exponential(1);
  mu_theta_logit ~ normal(0, 1);
  sigma_theta_logit ~ exponential(1);
  mu_Fexp ~ normal(0, 0.15);
  sigma_Fexp ~ exponential(1);
  mu_Fslope ~ normal(4, 1);
  sigma_Fslope ~ exponential(1);
  mu_Cexp ~ normal(0, 1);
  sigma_Cexp ~ exponential(1);
  mu_Cslope ~ normal(0, 1);
  sigma_Cslope ~ exponential(1);
  mu_log_qrec ~ normal(0, 1);
  sigma_log_qrec ~ exponential(1);
  mu_thetarec_logit ~ normal(0, 1);
  sigma_thetarec_logit ~ exponential(1);
  mu_activation_intercept ~ normal(-1, 1);    // Baseline logit-activation (e.g. -1 -> ~27% activation)
  sigma_activation_intercept ~ exponential(1);
  mu_activation_slope ~ normal(0, 0.1);     // How activation changes with log_L
  sigma_activation_slope ~ exponential(8);
  
  // Hyperpriors for missing data with no regressors
  mu_log_L ~ normal(6, 1);
  sigma_log_L ~ exponential(1);
  mu_log_Erec ~ normal(6, 1);
  sigma_log_Erec ~ exponential(1);

  // Priors for location-specific parameters (drawn from hyperparameters)
  alpha_raw ~ std_normal();
  beta_raw ~ std_normal();
  log_q_raw ~ std_normal();
  theta_logit_raw ~ std_normal();
  Fexp_raw ~ std_normal();
  Fslope_raw ~ std_normal();
  Cexp_raw ~ std_normal();
  Cslope_raw ~ std_normal();
  log_qrec_raw ~ std_normal();
  thetarec_logit_raw ~ std_normal();
  activation_intercept_raw ~ std_normal();
  activation_slope_raw ~ std_normal();

  // Priors for all latent values - missing data with no regressors
  for (l in 1:N_loc) {
    log_L[l] ~ normal(mu_log_L, sigma_log_L);
  }
  for (l in 1:N_loc) {
    log_Erec[l] ~ normal(mu_log_Erec, sigma_log_Erec);
  }

  // Priors for non-hierarchical location-specific parameters
  sigma_r ~ exponential(4);
  sigma_s ~ exponential(4);
  sigma_c ~ exponential(4);
  sigma_e ~ exponential(4);
  sigma_f ~ exponential(4);
  sigma_crec ~ exponential(4);
  sigma_erec ~ exponential(4);
  sigma_v ~ exponential(4);
  sigma_l ~ exponential(4);
  df_c ~ gamma(2, 0.1);  // A weakly informative prior for df_c
  for (l in 1:N_loc) {
      log_R_Perror[l][1:N_year_max] ~ normal(0, sigma_r[l]);
      start_log_R[l][1:genlength[l]] ~ normal(4, 2);    // Starting values for each location
    for (i in (genlength[l] + 1):N_year_max) {
      start_log_R[l][i] ~ normal(4, 2); // Or some other prior for the unused parts
    }
  }
  
  // Likelihood for observed fleet size
  for (n_f in 1:N_Fobs_total) {
    log_Fobs[n_f] ~ normal(log_F[Fobs_location_id[n_f]][Fobs_year_index[n_f]], sigma_f[Fobs_location_id[n_f]]);
  }
  
  // Likelihood for observed license counts
  for (n_f in 1:N_Lobs_total) {
    log_Lobs[n_f] ~ normal(log_L[Lobs_location_id[n_f]][Lobs_year_index[n_f]], sigma_l[Lobs_location_id[n_f]]);
  }

  // Observation model for value
  for (n in 1:N_Vobs_total) {
     if (log_Vobs[n] > 0){ // Assuming 0 indicates missing or true zero not handled by log-normal
        log_Vobs[n] ~ normal(log_V[Vobs_location_id[n]][Vobs_year_index[n]], sigma_v[Vobs_location_id[n]]);
     }
  }

  // Observation model for escapement
  for (n in 1:N_Sobs_total) {
    log_Sobs[n] ~ normal(log_S[Sobs_location_id[n]][Sobs_year_index[n]], sigma_s[Sobs_location_id[n]]);
  }

  // Hurdle observation model for catch
  for (n in 1:N_Cobs_total) {
    int loc = Cobs_location_id[n];
    int year = Cobs_year_index[n];
    if (log_Cobs[n] == 0) { // Assuming 0 indicates a true zero catch, not missing data
      target += log(theta[loc]);
    } else {
      target += log1m(theta[loc]) + student_t_lpdf(log_Cobs[n] | df_c[loc], log_C[loc][year], sigma_c[loc]);
    }
  }
  
  // Observation model for effort
  for (n in 1:N_total) {
    log_Eobs[n] ~ normal(log_E[location_id[n]][year_index[n]], sigma_e[location_id[n]]);
  }
  
  // Hurdle observation model for recreational catch
  for (n in 1:N_Recobs) {
    int loc = Recobs_location_id[n];
    int year = Recobs_year_index[n];
    if (log_RecCobs[n] == 0) {
      target += log(thetarec[loc]);
    } else {
      target += log1m(thetarec[loc]) + normal_lpdf(log_RecCobs[n] | log_Crec[loc][year], sigma_crec[loc]);
    }
    target += normal_lpdf(log_RecEobs[n] | log_Erec[loc][year], sigma_erec[loc]);
  }
}

generated quantities {
  // Produce some model predictions
  array[N_loc] vector[N_year_max] log_C_rep;
  array[N_loc] vector[N_year_max] log_S_rep;
  array[N_loc] vector[N_year_max] log_E_rep;
  array[N_loc] vector[N_year_max] log_F_rep;
  array[N_loc] vector[N_year_max] log_V_rep;
  array[N_loc] vector[N_year_max] log_L_rep;
  array[N_loc] vector[N_year_max] log_Crec_rep;
  array[N_loc] vector[N_year_max] log_Erec_rep;
  // Calculate pointwise log likelihoods
  array[N_loc] vector[N_year_max] log_lik_catch;
  array[N_loc] vector[N_year_max] log_lik_effort;
  array[N_loc] vector[N_year_max] log_lik_spawners;
  array[N_loc] vector[N_year_max] log_lik_fleet;
  array[N_loc] vector[N_year_max] log_lik_value;
  array[N_loc] vector[N_year_max] log_lik_reccatch;
  array[N_loc] vector[N_year_max] log_lik_receffort;
  array[N_loc] vector[N_year_max] log_lik_licenses;

  for (l in 1:N_loc) {
    log_C_rep[l] = rep_vector(negative_infinity(), N_year_max);
    log_S_rep[l] = rep_vector(0.0, N_year_max);
    log_E_rep[l] = rep_vector(0.0, N_year_max);
    log_F_rep[l] = rep_vector(0.0, N_year_max);
    log_V_rep[l] = rep_vector(0.0, N_year_max);
    log_Crec_rep[l] = rep_vector(negative_infinity(), N_year_max);
    log_Erec_rep[l] = rep_vector(negative_infinity(), N_year_max);

    // Initialize log_lik arrays for the current location
    log_lik_catch[l] = rep_vector(0.0, N_year_max); // Or negative_infinity()
    log_lik_effort[l] = rep_vector(0.0, N_year_max); // Or negative_infinity()
    log_lik_spawners[l] = rep_vector(0.0, N_year_max); // Or negative_infinity()
    log_lik_fleet[l] = rep_vector(0.0, N_year_max); // Or negative_infinity()
    log_lik_value[l] = rep_vector(0.0, N_year_max);
    log_lik_reccatch[l] = rep_vector(0.0, N_year_max); // Or negative_infinity()
    log_lik_receffort[l] = rep_vector(0.0, N_year_max); // Or negative_infinity()
    log_lik_licenses[l] = rep_vector(0.0, N_year_max); // Or negative_infinity()

    // Predictive distribution for catch, effort, fleet, value and harvest for location 'l'
    for (y in 1:N_year_max) {
        log_C_rep[l][y] = student_t_rng(df_c[l], log_C[l][y], sigma_c[l]);
        log_S_rep[l][y] = normal_rng(log_S[l][y], sigma_s[l]);
        log_E_rep[l][y] = normal_rng(log_E[l][y], sigma_e[l]);
        log_F_rep[l][y] = normal_rng(log_F[l][y], sigma_f[l]);
        log_V_rep[l][y] = normal_rng(log_V[l][y], sigma_v[l]);
        log_L_rep[l][y] = normal_rng(log_L[l][y], sigma_l[l]);
        log_Crec_rep[l][y] = normal_rng(log_Crec[l][y], sigma_crec[l]);
        log_Erec_rep[l][y] = normal_rng(log_Erec[l][y], sigma_erec[l]);
    }
    // Log likelihood for observed effort
    for (n in 1:N_total) {
      if (location_id[n] == l && year_index[n] <= N_year_max) {
        log_lik_effort[l][year_index[n]] = normal_lpdf(log_Eobs[n] | log_E[l][year_index[n]], sigma_e[l]);
      }
    }
    
    // Log likelihood for observed catch
    for (n in 1:N_Cobs_total) {
      if (Cobs_location_id[n] == l && Cobs_year_index[n] <= N_year_max) {
        if (log_Cobs[n] == 0) {
          log_lik_catch[l][Cobs_year_index[n]] = log(theta[l]);
        } else {
          log_lik_catch[l][Cobs_year_index[n]] = log1m(theta[l]) + student_t_lpdf(log_Cobs[n] | df_c[l], log_C[l][Cobs_year_index[n]], sigma_c[l]);
        }
      }
    }
    
    // Log likelihood for observed spawners
    for (n_s in 1:N_Sobs_total) {
      if (Sobs_location_id[n_s] == l && Sobs_year_index[n_s] <= N_year_max) {
        log_lik_spawners[l][Sobs_year_index[n_s]] = normal_lpdf(log_Sobs[n_s] | log_S[l][Sobs_year_index[n_s]], sigma_s[l]);
      }
    }
    // Log likelihood for observed fleet size
    for (n_f in 1:N_Fobs_total) {
      if (Fobs_location_id[n_f] == l && Fobs_year_index[n_f] <= N_year_max) {
        log_lik_fleet[l][Fobs_year_index[n_f]] = normal_lpdf(log_Fobs[n_f] | log_F[l][Fobs_year_index[n_f]], sigma_f[l]);
      }
    }
    // Log likelihood for observed recreational catch
    for (n_l in 1:N_Lobs_total) {
      if (Lobs_location_id[n_l] == l && Lobs_year_index[n_l] <= N_year_max) {
        log_lik_licenses[l][Lobs_year_index[n_l]] = normal_lpdf(log_Lobs[n_l] | log_L[l][Lobs_year_index[n_l]], sigma_l[l]);
      }
    }
    // Log likelihood for observed recreational effort
    for (n_rece in 1:N_Recobs) {
      if (Recobs_location_id[n_rece] == l && Recobs_year_index[n_rece] <= N_year_max) {
        log_lik_receffort[l][Recobs_year_index[n_rece]] = normal_lpdf(log_RecEobs[n_rece] | log_Erec[l][Recobs_year_index[n_rece]], sigma_erec[l]);
      }
    }
    // Log likelihood for observed recreational catch 
    for (n_rec in 1:N_Recobs) {
      if (Recobs_location_id[n_rec] == l && Recobs_year_index[n_rec] <= N_year_max) {
        if (log_RecCobs[n_rec] == 0) {
          log_lik_reccatch[l][Recobs_year_index[n_rec]] = log(thetarec[l]);
        } else {
          log_lik_reccatch[l][Recobs_year_index[n_rec]] = log1m(thetarec[l]) + normal_lpdf(log_RecCobs[n_rec] | log_Crec[l][Recobs_year_index[n_rec]], sigma_crec[l]);
        }
      }
    }
    // Log likelihood for observed value
    for (n_v in 1:N_Vobs_total) {
      if (log_Vobs[n_v] > 0) {
        if (Vobs_location_id[n_v] == l && Vobs_year_index[n_v] <= N_year_max) {
          log_lik_value[l][Vobs_year_index[n_v]] = normal_lpdf(log_Vobs[n_v] | log_V[l][Vobs_year_index[n_v]], sigma_v[l]);
        }
      }
    }
  }
}
