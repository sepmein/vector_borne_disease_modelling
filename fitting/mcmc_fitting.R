# ===============================================================================
# SEIR-SIR Vector-Borne Disease Model - MCMC Fitting and Analysis
# ===============================================================================
# Purpose: Bayesian inference for dengue transmission dynamics using Stan
# Model: 7-compartment SEIR-SIR (Human: S-E-I-R, Vector: S-I-R)
# Date: Created for Foshan dengue outbreak analysis
# ===============================================================================

# ===============================================================================
# 1. ENVIRONMENT SETUP AND CONFIGURATION
# ===============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(rstan)
  library(tidybayes)
  library(tidyverse)
  library(ggplot2)
  library(posterior)
  library(bayesplot)
  library(lubridate)
})

# Configure Stan settings for optimal performance
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set random seed for reproducibility
set.seed(42)

# ===============================================================================
# 2. DATA PREPARATION AND PREPROCESSING
# ===============================================================================

# Define observation period and case data
setup_data <- function() {
  cat("Loading and preprocessing epidemiological data...\n")
  
  # Time series data - Foshan dengue outbreak 2025
  dates <- as.Date(c(
    "2025-07-08", "2025-07-09", "2025-07-10", "2025-07-11",
    "2025-07-12", "2025-07-13", "2025-07-14", "2025-07-15",
    "2025-07-16", "2025-07-17", "2025-07-18", "2025-07-19",
    "2025-07-20", "2025-07-21", "2025-07-22"
  ))
  # Daily reported cases
  cases <- c(1, 1, 3, 8, 20, 48, 116, 280, 104, 296, 393, 473, 540, 373, 536)

  # Basic statistics
  cat(sprintf(
    "Observation period: %s to %s (%d days)\n",
    min(dates), max(dates), length(dates)
  ))
  cat(sprintf(
    "Total cases: %d (Peak: %d cases on %s)\n",
    sum(cases), max(cases), dates[which.max(cases)]
  ))

  # Convert to Stan format
  T <- length(cases)
  ts <- as.numeric(dates - dates[1]) + 1 # Time points starting from 1
  return(list(
    dates = dates,
    cases = cases,
    T = T,
    ts = ts
  ))
}

# ===============================================================================
# 3. MODEL PARAMETERS AND INITIAL CONDITIONS
# ===============================================================================

# Set up initial state and population parameters
setup_model_parameters <- function() {
  cat("Setting up model parameters and initial conditions...\n")
  
  # Population parameters
  N_h <- 969.89e4  # Total human population (Foshan city)
  
  # Initial conditions for 7-compartment model
  # Human compartments: S_h, E_h, I_h, R_h
  # Vector compartments: S_v, I_v, R_v
  init_state <- c(
    N_h - 1,    # S_h0: Initially susceptible humans (nearly all population)
    0,          # E_h0: Initially exposed humans
    1,          # I_h0: Initially infectious humans (index case)
    0,          # R_h0: Initially recovered humans
    1e5,        # S_v0: Initially susceptible vectors (adjustable)
    0,          # I_v0: Initially infectious vectors
    0           # R_v0: Initially recovered vectors
  )
  
  cat(sprintf("Human population: %.2e\n", N_h))
  cat(sprintf("Initial infectious humans: %d\n", init_state[3]))
  cat(sprintf("Initial susceptible vectors: %.0e\n", init_state[5]))
  
  return(list(
    N_h = N_h,
    init_state = init_state
  ))
}

# ===============================================================================
# 4. BAYESIAN INFERENCE WITH STAN
# ===============================================================================

# Run MCMC sampling using Stan
run_mcmc_inference <- function(stan_data, n_iter = 2000, n_chains = 12, 
                               adapt_delta = 0.9, max_treedepth = 15) {
  cat("Starting Bayesian inference with Stan...\n")
  cat(sprintf("Chains: %d, Iterations: %d per chain\n", n_chains, n_iter))
  cat(sprintf("Total posterior samples: %d (after warmup)\n", n_chains * n_iter/2))
  
  # Check if compiled model exists to save compilation time
  stan_file <- "fitting/seir_zxc.stan"
  
  if (!file.exists(stan_file)) {
    stop("Stan model file not found: ", stan_file)
  }
  
  # Run Stan with enhanced control parameters
  fit <- stan(
    file = stan_file,
    data = stan_data,
    iter = n_iter,
    chains = n_chains,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    verbose = TRUE
  )
  
  # Basic convergence diagnostics
  cat("\n=== MCMC Diagnostics ===\n")
  cat("Rhat summary (should be < 1.1):\n")
  rhat_vals <- rhat(fit)
  print(summary(rhat_vals[!is.na(rhat_vals)]))
  
  cat("\nEffective sample size summary:\n")
  eff_vals <- eff_sample(fit)
  print(summary(eff_vals[!is.na(eff_vals)]))
  
  return(fit)
}

# ===============================================================================
# 4.5. CHAIN CONVERGENCE ANALYSIS (NEW SECTION)
# ===============================================================================

# Detailed chain convergence analysis and visualization
analyze_chain_convergence <- function(fit) {
  cat("\n=== Detailed Chain Convergence Analysis ===\n")
  
  # Extract key parameters for convergence analysis
  key_params <- c("beta", "beta_hv", "sigma_h", "gamma_h", "vie", "beta_vh", "gamma_v")
  
  # Create trace plots to visualize chain mixing
  p_trace <- mcmc_trace(fit, pars = key_params[1:4])  # Show first 4 parameters
  
  # Create Rhat plot
  p_rhat <- mcmc_rhat(rhat(fit))
  
  # Create effective sample size plot  
  p_neff <- mcmc_neff(neff_ratio(fit))
  
  # Chain-specific statistics
  chain_summary <- fit %>%
    as_draws_df() %>%
    select(.chain, all_of(key_params)) %>%
    group_by(.chain) %>%
    summarise(
      across(all_of(key_params), list(mean = mean, sd = sd)),
      .groups = "drop"
    )
  
  cat("\n--- Chain-specific Statistics ---\n")
  print(chain_summary)
  
  # Check for chain convergence
  rhat_vals <- rhat(fit, pars = key_params)
  n_divergent <- sum(get_divergent_iterations(fit))
  
  cat(sprintf("\nConvergence Summary:\n"))
  cat(sprintf("- Max Rhat: %.4f (should be < 1.1)\n", max(rhat_vals, na.rm = TRUE)))
  cat(sprintf("- Divergent transitions: %d (should be 0)\n", n_divergent))
  
  if (max(rhat_vals, na.rm = TRUE) < 1.1 && n_divergent == 0) {
    cat(" All chains converged successfully!\n")
  } else {
    cat("  Convergence issues detected. Consider:\n")
    cat("   - Increasing adapt_delta (e.g., 0.95)\n")
    cat("   - Running more iterations\n")
    cat("   - Checking model specification\n")
  }
  
  return(list(
    trace_plot = p_trace,
    rhat_plot = p_rhat,
    neff_plot = p_neff,
    chain_summary = chain_summary
  ))
}

# ===============================================================================
# 5. RESULTS ANALYSIS AND VISUALIZATION
# ===============================================================================

# Extract and summarize parameter estimates
analyze_parameters <- function(fit) {
  cat("\n=== Parameter Estimates ===\n")
  
  # Key epidemiological parameters
  key_params <- c("beta", "beta_hv", "sigma_h", "gamma_h", "vie", "beta_vh", "gamma_v", "phi")
  
  # Print parameter summary
  print(fit, pars = key_params, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  # Extract draws for further analysis
  draws <- as_draws_df(fit)
  
  # Calculate derived quantities
  param_summary <- draws %>%
    select(all_of(key_params)) %>%
    summarise(
      across(everything(), list(
        mean = mean,
        median = median,
        sd = sd,
        q025 = ~quantile(.x, 0.025),
        q975 = ~quantile(.x, 0.975)
      ))
    )
  
  return(list(
    draws = draws,
    summary = param_summary
  ))
}

# Create comprehensive visualization
create_model_plots <- function(fit, data_list) {
  cat("\n=== Creating Visualizations ===\n")
  
  # Extract posterior predictions for incidence
  posterior_inc <- fit %>%
    as_draws_df() %>% 
    select(starts_with("incidence")) %>% 
    pivot_longer(everything(),
                 names_to = "day",
                 values_to = "value") %>%
    mutate(
      day = as.integer(str_extract(day, "\\d+")),
      date = data_list$dates[day]
    )
  # Calculate credible intervals
  summ_inc <- posterior_inc %>%
    group_by(date, day) %>%
    summarise(
      median = median(value),
      mean = mean(value),
      lower_50 = quantile(value, 0.25),
      upper_50 = quantile(value, 0.75),
      lower_95 = quantile(value, 0.025),
      upper_95 = quantile(value, 0.975),
      .groups = "drop"
    )
  
  # Main prediction plot
  p1 <- ggplot() +
    # 95% credible interval
    geom_ribbon(data = summ_inc,
                aes(x = date, ymin = lower_95, ymax = upper_95),
                fill = "steelblue", alpha = 0.2) +
    # 50% credible interval
    geom_ribbon(data = summ_inc,
                aes(x = date, ymin = lower_50, ymax = upper_50),
                fill = "steelblue", alpha = 0.4) +
    # Posterior median
    geom_line(data = summ_inc,
              aes(x = date, y = median),
              color = "steelblue", linewidth = 1.2) +
    # Observed data points
    geom_point(data = tibble(date = data_list$dates, cases = data_list$cases),
               aes(x = date, y = cases),
               color = "red", size = 2, alpha = 0.8) +
    scale_x_date(date_labels = "%m-%d", date_breaks = "2 days") +
    labs(
      title = "SEIR-SIR Model: Posterior Predictions vs Observed Cases",
      subtitle = "Red points: observed data, Blue: model predictions with credible intervals",
      x = "Date (2025)", 
      y = "Daily New Cases",
      caption = "50% (dark) and 95% (light) posterior credible intervals"
    ) +
    theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
  # Parameter posterior distributions
  draws <- as_draws_df(fit)
  key_params <- c("beta", "beta_hv", "sigma_h", "gamma_h", "vie", "beta_vh", "gamma_v")
  
  p2 <- draws %>%
    select(all_of(key_params)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 50, fill = "lightblue", alpha = 0.7, color = "black") +
    facet_wrap(~parameter, scales = "free", ncol = 3) +
    labs(
      title = "Posterior Distributions of Key Parameters",
      x = "Parameter Value",
      y = "Frequency"
    ) +
    theme_minimal()
  
  return(list(
    prediction_plot = p1,
    parameter_plot = p2,
    summary_data = summ_inc
  ))
}

# ===============================================================================
# 6. MODEL VALIDATION AND DIAGNOSTICS
# ===============================================================================

# Perform model checking and validation
model_diagnostics <- function(fit, data_list) {
  cat("\n=== Model Diagnostics ===\n")
  
  # Posterior predictive checks
  y_rep <- extract(fit, "incidence")[[1]]
  y_obs <- data_list$cases
  
  # Calculate Bayesian R-squared
  # (simplified version - could be enhanced)
  residuals <- sweep(y_rep, 2, y_obs, "-")
  ss_res <- rowSums(residuals^2)
  ss_tot <- sum((y_obs - mean(y_obs))^2)
  r_squared <- 1 - ss_res/ss_tot
  
  cat(sprintf("Bayesian R-squared: %.3f [%.3f, %.3f]\n", 
              median(r_squared), quantile(r_squared, 0.025), quantile(r_squared, 0.975)))
  
  # Root Mean Square Error
  rmse <- sqrt(rowMeans(residuals^2))
  cat(sprintf("RMSE: %.2f [%.2f, %.2f]\n", 
              median(rmse), quantile(rmse, 0.025), quantile(rmse, 0.975)))
  
  # Mean Absolute Error
  mae <- rowMeans(abs(residuals))
  cat(sprintf("MAE: %.2f [%.2f, %.2f]\n", 
              median(mae), quantile(mae, 0.025), quantile(mae, 0.975)))
  
  return(list(
    r_squared = r_squared,
    rmse = rmse,
    mae = mae
  ))
}

# ===============================================================================
# 7. MAIN EXECUTION WORKFLOW
# ===============================================================================

# Main analysis function
main_analysis <- function() {
  cat("===============================================================================\n")
  cat("SEIR-SIR Bayesian Analysis Pipeline\n")
  cat("===============================================================================\n")
  
  # Explain chains concept first
  explain_chains_concept()
  
  # Step 1: Data preparation
  data_list <- setup_data()
  
  # Step 2: Model setup
  model_params <- setup_model_parameters()
  
  # Step 3: Prepare Stan data
  stan_data <- list(
    T = data_list$T,
    ts = data_list$ts,
    cases = data_list$cases,
    N_h = model_params$N_h,
    init_state = model_params$init_state
  )
  
  # Step 4: Run MCMC inference
  fit <- run_mcmc_inference(stan_data)
  
  # Step 4.5: Analyze chain convergence
  convergence_analysis <- analyze_chain_convergence(fit)
  
  # Step 5: Analyze results
  analysis_results <- analyze_parameters(fit)
  
  # Step 6: Create visualizations
  plots <- create_model_plots(fit, data_list)
  
  # Step 7: Model diagnostics
  diagnostics <- model_diagnostics(fit, data_list)
  
  # Display main plot
  print(plots$prediction_plot)
  
  cat("\n===============================================================================\n")
  cat("Analysis completed successfully!\n")
  cat("===============================================================================\n")
  
  return(list(
    fit = fit,
    data = data_list,
    results = analysis_results,
    plots = plots,
    diagnostics = diagnostics,
    convergence = convergence_analysis
  ))
}

# ===============================================================================
# 8. EXECUTION
# ===============================================================================

# Run the complete analysis
# Uncomment the line below to execute
results <- main_analysis()

# For interactive use, you can run individual components:
# data_list <- setup_data()
# model_params <- setup_model_parameters()
# ...and so on


