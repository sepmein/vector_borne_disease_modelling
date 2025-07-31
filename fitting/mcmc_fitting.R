# ===============================================================================
# SEIR-SEI Vector-Borne Disease Model - MCMC Fitting and Analysis
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
  library(EpiEstim) # Added for Rt estimation
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

  # Read symptom-based cases data (biological incidence)
  symptons_file <- "data/daily_cases_foshan_symptons.csv"
  reported_file <- "data/daily_cases_foshan.csv"

  if (!file.exists(symptons_file)) {
    stop("Symptoms CSV file not found: ", symptons_file)
  }
  if (!file.exists(reported_file)) {
    stop("Reported CSV file not found: ", reported_file)
  }

  # Read both CSV files
  symptons_df <- read.csv(symptons_file, stringsAsFactors = FALSE)
  reported_df <- read.csv(reported_file, stringsAsFactors = FALSE)

  # Remove 2025-07-26 data from both datasets
  symptons_df <- symptons_df[symptons_df$date != "2025-07-26", ]
  reported_df <- reported_df[reported_df$date != "2025-07-26", ]

  # Convert date columns to Date format
  symptons_dates <- as.Date(symptons_df$date)
  reported_dates <- as.Date(reported_df$date)

  # Find overlapping date range
  start_date <- max(min(symptons_dates), min(reported_dates))
  end_date <- min(max(symptons_dates), max(reported_dates))

  cat(sprintf("Symptoms data range: %s to %s\n", min(symptons_dates), max(symptons_dates)))
  cat(sprintf("Reported data range: %s to %s\n", min(reported_dates), max(reported_dates)))
  cat(sprintf("Using overlapping period: %s to %s\n", start_date, end_date))

  # Filter both datasets to overlapping period
  symptons_filtered <- symptons_df[symptons_dates >= start_date & symptons_dates <= end_date, ]
  reported_filtered <- reported_df[reported_dates >= start_date & reported_dates <= end_date, ]

  # Create aligned date sequence
  dates <- seq(from = start_date, to = end_date, by = "day")

  # Align both datasets to the same date sequence
  symptons_aligned <- merge(
    data.frame(date = as.character(dates)),
    symptons_filtered,
    by = "date",
    all.x = TRUE
  )
  reported_aligned <- merge(
    data.frame(date = as.character(dates)),
    reported_filtered,
    by = "date",
    all.x = TRUE
  )

  # Fill missing values with 0
  symptons_aligned$daily_cases[is.na(symptons_aligned$daily_cases)] <- 0
  reported_aligned$daily_cases[is.na(reported_aligned$daily_cases)] <- 0

  # Extract aligned case vectors
  cases_symptons <- symptons_aligned$daily_cases
  cases_reported <- reported_aligned$daily_cases
  # Basic statistics
  cat(sprintf(
    "Final observation period: %s to %s (%d days)\n",
    min(dates), max(dates), length(dates)
  ))
  cat(sprintf(
    "Total symptom cases: %d (Peak: %d cases on %s)\n",
    sum(cases_symptons), max(cases_symptons), dates[which.max(cases_symptons)]
  ))
  cat(sprintf(
    "Total reported cases: %d (Peak: %d cases on %s)\n",
    sum(cases_reported), max(cases_reported), dates[which.max(cases_reported)]
  ))

  # Convert to Stan format
  T <- length(dates)
  ts <- as.numeric(dates - dates[1]) + 1 # Time points starting from 1

  return(list(
    dates = dates,
    cases = cases_symptons,           # Biological incidence (E→I transition)
    cases_reported = cases_reported,  # Reported cases
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
  N_h <- 969.89e4 # Total human population (Foshan city)

  # More realistic initial conditions for dengue outbreak
  # Based on early case detection and vector dynamics
  init_state <- c(
    N_h - 10, # S_h0: Initially susceptible humans (assume some early exposure)
    5, # E_h0: Initially exposed humans (realistic for dengue)
    5, # I_h0: Initially infectious humans (matches early cases)
    0, # R_h0: Initially recovered humans
    0, # Report_h0: Initially reported cumulative cases
    5000, # S_v0: Much smaller vector population (more realistic)
    10, # I_v0: Some initially infectious vectors
    0 # R_v0: Initially recovered vectors
  )

  cat(sprintf("Human population: %.2e\n", N_h))
  cat(sprintf("Initial infectious humans: %d\n", init_state[3]))
  cat(sprintf("Initial susceptible vectors: %.0e\n", init_state[6]))

  return(list(
    N_h = N_h,
    init_state = init_state
  ))
}

# ===============================================================================
# 4. BAYESIAN INFERENCE WITH STAN
# ===============================================================================

# Run MCMC sampling using Stan
run_mcmc_inference <- function(stan_data, n_iter = 2000, n_chains = 4,
                               adapt_delta = 0.99, max_treedepth = 15) {
  cat("Starting Bayesian inference with Stan...\n")
  cat(sprintf("Chains: %d, Iterations: %d per chain\n", n_chains, n_iter))
  cat(sprintf("Total posterior samples: %d (after warmup)\n", n_chains * n_iter / 2))

  # Check if compiled model exists to save compilation time
  stan_file <- "fitting/seir_zxc.stan"

  if (!file.exists(stan_file)) {
    stop("Stan model file not found: ", stan_file)
  }
  
  # Debug: Check Stan data structure before passing to Stan
  cat("\n=== DEBUG: Checking Stan Data Before Fitting ===\n")
  cat("Stan file:", stan_file, "\n")
  cat("Data names:", names(stan_data), "\n")
  cat("N_h value:", stan_data$N_h, "\n")
  cat("N_h class:", class(stan_data$N_h), "\n")

  # Create stable initial values for each chain
  create_stable_init <- function() {
    list(
      beta = 0.25 + rnorm(1, 0, 0.05), # Stable around 0.25
      beta_hv = 0.24 + rnorm(1, 0, 0.03), # Stable around 0.24
      sigma_h = 1 / 3 + rnorm(1, 0, 0.02), # Stable around 1/3
      gamma_h = 1 / 4 + rnorm(1, 0, 0.02), # Stable around 1/4
      vie = 0.25 + rnorm(1, 0, 0.05), # Stable around 0.25
      beta_vh = 0.24 + rnorm(1, 0, 0.03), # Stable around 0.24
      gamma_v = 1 / 3.5 + rnorm(1, 0, 0.02), # Stable around 1/3.5
      gamma_rh = 1 / 3.5 + rnorm(1, 0, 0.05), # Stable around 1/3.5
      awareness = 0.1 + abs(rnorm(1, 0, 0.02)), # Small positive
      mu_v = 1 / 10 + rnorm(1, 0, 0.01), # Stable around 1/10
      phi_incidence = 10 + abs(rnorm(1, 0, 2)), # Positive around 10
      phi_reported = 5 + abs(rnorm(1, 0, 1)) # Positive around 5
    )
  }

  init_list <- lapply(1:n_chains, function(x) create_stable_init())

  # Run Stan with enhanced control parameters and stable initialization
  fit <- stan(
    file = stan_file,
    data = stan_data,
    iter = n_iter,
    chains = n_chains,
    init = init_list,  # ✅ ADD STABLE INITIAL VALUES
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      stepsize = 0.01,           # ✅ Smaller steps for stability
      metric = "diag_e"          # ✅ Diagonal mass matrix
    ),
    verbose = TRUE
  )

  # Basic convergence diagnostics
  cat("\n=== MCMC Diagnostics ===\n")
  cat("Rhat summary (should be < 1.1):\n")
  rhat_vals <- rhat(fit)
  print(summary(rhat_vals[!is.na(rhat_vals)]))

  cat("\nEffective sample size summary:\n")
  eff_vals <- summary(fit)$summary[, "n_eff"]
  print(summary(eff_vals[!is.na(eff_vals)]))

  return(fit)
}

# ===============================================================================
# 4.5. CHAIN CONVERGENCE ANALYSIS
# ===============================================================================

# Detailed chain convergence analysis and visualization
analyze_chain_convergence <- function(fit) {
  cat("\n=== Detailed Chain Convergence Analysis ===\n")

  # Extract key parameters for convergence analysis
  key_params <- c("beta", "beta_hv", "sigma_h", "gamma_h", "vie", "beta_vh", "gamma_v", "gamma_rh", "awareness", "mu_v")

  # Create trace plots to visualize chain mixing
  p_trace <- mcmc_trace(fit, pars = key_params[1:4]) # Show first 4 parameters

  # Create Rhat plot
  p_rhat <- mcmc_rhat(rhat(fit))
# Create effective sample size plot
  p_neff <- mcmc_neff(summary(fit)$summary[, "n_eff"])

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
  key_params <- c("beta", "beta_hv", "sigma_h", "gamma_h", "vie", "beta_vh", "gamma_v", "gamma_rh", "awareness", "mu_v", "phi_incidence", "phi_reported")

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
        q025 = ~ quantile(.x, 0.025),
        q975 = ~ quantile(.x, 0.975)
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
      values_to = "value"
    ) %>%
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
    geom_ribbon(
      data = summ_inc,
      aes(x = date, ymin = lower_95, ymax = upper_95),
      fill = "steelblue", alpha = 0.2
    ) +
    # 50% credible interval
    geom_ribbon(
      data = summ_inc,
      aes(x = date, ymin = lower_50, ymax = upper_50),
      fill = "steelblue", alpha = 0.4
    ) +
    # Posterior median
    geom_line(
      data = summ_inc,
      aes(x = date, y = median),
      color = "steelblue", linewidth = 1.2
    ) +
    # Observed data points
    geom_point(
      data = tibble(date = data_list$dates, cases = data_list$cases),
      aes(x = date, y = cases),
      color = "red", size = 2, alpha = 0.8
    ) +
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
  key_params <- c("beta", "beta_hv", "sigma_h", "gamma_h", "vie", "beta_vh", "gamma_v", "gamma_rh", "awareness", "mu_v")

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

  # Future prediction extraction and plots
  future_pred <- NULL
  future_reported_pred <- NULL
  p_future <- NULL
  p_future_reported <- NULL
  if ("future_incidence" %in% names(rstan::extract(fit))) {
    # Extract biological incidence predictions
    future_inc <- rstan::extract(fit, pars = "future_incidence")$future_incidence
    future_pred <- as.data.frame(t(apply(future_inc, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
    colnames(future_pred) <- c("lower_95", "lower_50", "median", "upper_50", "upper_95")
    future_pred$day <- 1:nrow(future_pred)
    future_pred$date <- max(data_list$dates) + future_pred$day
    future_pred$type <- "Biological Incidence"

    # Plot biological incidence predictions
    p_future <- ggplot(future_pred, aes(x = date)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "green", alpha = 0.2) +
      geom_ribbon(aes(ymin = lower_50, ymax = upper_50), fill = "green", alpha = 0.4) +
      geom_line(aes(y = median), color = "green") +
      labs(title = "Future Biological Incidence Prediction (Next 7 Days)", 
           y = "Predicted Daily Biological Cases", x = "Date") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
  }
  
  if ("future_reported" %in% names(rstan::extract(fit))) {
    # Extract reported cases predictions
    future_rep <- rstan::extract(fit, pars = "future_reported")$future_reported
    future_reported_pred <- as.data.frame(t(apply(future_rep, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
    colnames(future_reported_pred) <- c("lower_95", "lower_50", "median", "upper_50", "upper_95")
    future_reported_pred$day <- 1:nrow(future_reported_pred)
    future_reported_pred$date <- max(data_list$dates) + future_reported_pred$day
    future_reported_pred$type <- "Reported Cases"

    # Plot reported cases predictions
    p_future_reported <- ggplot(future_reported_pred, aes(x = date)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "orange", alpha = 0.2) +
      geom_ribbon(aes(ymin = lower_50, ymax = upper_50), fill = "orange", alpha = 0.4) +
      geom_line(aes(y = median), color = "orange") +
      labs(title = "Future Reported Cases Prediction (Next 7 Days)", 
           y = "Predicted Daily Reported Cases", x = "Date") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
  }

  # Combine model fit and future prediction for plotting
  model_pred <- summ_inc %>%
    select(date, median, lower_50, upper_50, lower_95, upper_95) %>%
    mutate(type = "Model Fit")

  # Combine biological incidence predictions if available
  all_pred <- model_pred
  if (!is.null(future_pred)) {
    future_pred2 <- future_pred %>%
      select(date, median, lower_50, upper_50, lower_95, upper_95) %>%
      mutate(type = "Future Biological Incidence")
    all_pred <- bind_rows(all_pred, future_pred2)
  }

  return(list(
    prediction_plot = p1,
    parameter_plot = p2,
    summary_data = summ_inc,
    future_prediction_plot = p_future,
    future_reported_plot = p_future_reported,
    future_prediction_summary = future_pred,
    future_reported_summary = future_reported_pred,
    combined_prediction_plot = ggplot() +
      # 95% credible interval ribbons
      geom_ribbon(
        data = all_pred,
        aes(x = date, ymin = lower_95, ymax = upper_95, fill = type),
        alpha = 0.2
      ) +
      # 50% credible interval ribbons
      geom_ribbon(
        data = all_pred,
        aes(x = date, ymin = lower_50, ymax = upper_50, fill = type),
        alpha = 0.4
      ) +
      # Posterior median lines
      geom_line(
        data = all_pred,
        aes(x = date, y = median, color = type),
        size = 1.2
      ) +
      # Observed data points
      geom_point(
        data = tibble(date = data_list$dates, cases = data_list$cases),
        aes(x = date, y = cases),
        color = "red", size = 2, alpha = 0.8, shape = 16
      ) +
      scale_color_manual(values = c("Model Fit" = "steelblue", "Future Biological Incidence" = "green")) +
      scale_fill_manual(values = c("Model Fit" = "steelblue", "Future Biological Incidence" = "green")) +
      labs(
        title = "Observed Data, Model Fit, and Future Prediction",
        subtitle = "Red points: observed data; Blue: model fit; Green: future prediction",
        x = "Date",
        y = "Daily Cases"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
  ))
}

# ===============================================================================
# 6. MODEL VALIDATION AND DIAGNOSTICS
# ===============================================================================

# Perform model checking and validation
model_diagnostics <- function(fit, data_list) {
  cat("\n=== Model Diagnostics ===\n")

  # Posterior predictive checks for biological incidence
  cat("\n--- Biological Incidence Diagnostics ---\n")
  y_rep_inc <- rstan::extract(fit, pars = "incidence")$incidence
  y_obs_inc <- data_list$cases

  # Calculate Bayesian R-squared for incidence
  residuals_inc <- sweep(y_rep_inc, 2, y_obs_inc, "-")
  ss_res_inc <- rowSums(residuals_inc^2)
  ss_tot_inc <- sum((y_obs_inc - mean(y_obs_inc))^2)
  r_squared_inc <- 1 - ss_res_inc / ss_tot_inc

  cat(sprintf(
    "Biological Incidence - Bayesian R-squared: %.3f [%.3f, %.3f]\n",
    median(r_squared_inc), quantile(r_squared_inc, 0.025), quantile(r_squared_inc, 0.975)
  ))

  # Root Mean Square Error for incidence
  rmse_inc <- sqrt(rowMeans(residuals_inc^2))
  cat(sprintf(
    "Biological Incidence - RMSE: %.2f [%.2f, %.2f]\n",
    median(rmse_inc), quantile(rmse_inc, 0.025), quantile(rmse_inc, 0.975)
  ))

  # Mean Absolute Error for incidence
  mae_inc <- rowMeans(abs(residuals_inc))
  cat(sprintf(
    "Biological Incidence - MAE: %.2f [%.2f, %.2f]\n",
    median(mae_inc), quantile(mae_inc, 0.025), quantile(mae_inc, 0.975)
  ))

  # Posterior predictive checks for reported cases
  cat("\n--- Reported Cases Diagnostics ---\n")
  y_rep_rep <- rstan::extract(fit, pars = "reported_pred")$reported_pred
  y_obs_rep <- data_list$cases_reported

  # Calculate Bayesian R-squared for reported cases
  residuals_rep <- sweep(y_rep_rep, 2, y_obs_rep, "-")
  ss_res_rep <- rowSums(residuals_rep^2)
  ss_tot_rep <- sum((y_obs_rep - mean(y_obs_rep))^2)
  r_squared_rep <- 1 - ss_res_rep / ss_tot_rep

  cat(sprintf(
    "Reported Cases - Bayesian R-squared: %.3f [%.3f, %.3f]\n",
    median(r_squared_rep), quantile(r_squared_rep, 0.025), quantile(r_squared_rep, 0.975)
  ))

  # Root Mean Square Error for reported cases
  rmse_rep <- sqrt(rowMeans(residuals_rep^2))
  cat(sprintf(
    "Reported Cases - RMSE: %.2f [%.2f, %.2f]\n",
    median(rmse_rep), quantile(rmse_rep, 0.025), quantile(rmse_rep, 0.975)
  ))

  # Mean Absolute Error for reported cases
  mae_rep <- rowMeans(abs(residuals_rep))
  cat(sprintf(
    "Reported Cases - MAE: %.2f [%.2f, %.2f]\n",
    median(mae_rep), quantile(mae_rep, 0.025), quantile(mae_rep, 0.975)
  ))

  # Combined diagnostics
  cat("\n--- Combined Model Performance ---\n")
  total_rmse <- sqrt((ss_res_inc + ss_res_rep) / (length(y_obs_inc) + length(y_obs_rep)))
  cat(sprintf(
    "Combined RMSE: %.2f [%.2f, %.2f]\n",
    median(total_rmse), quantile(total_rmse, 0.025), quantile(total_rmse, 0.975)
  ))

  return(list(
    # Biological incidence diagnostics
    r_squared_incidence = r_squared_inc,
    rmse_incidence = rmse_inc,
    mae_incidence = mae_inc,
    # Reported cases diagnostics  
    r_squared_reported = r_squared_rep,
    rmse_reported = rmse_rep,
    mae_reported = mae_rep,
    # Combined
    combined_rmse = total_rmse
  ))
}

# ===============================================================================
# 6.5. RT (EFFECTIVE REPRODUCTION NUMBER) ESTIMATION
# ===============================================================================

# Calculate comprehensive Rt for observed + predicted periods
calculate_rt_comprehensive <- function(fit, data_list) {
  cat("\n=== Calculating Comprehensive Rt (Observed + Predicted) ===\n")

  # Extract posterior predictions for incidence (observed period)
  posterior_inc <- fit %>%
    as_draws_df() %>%
    select(starts_with("incidence")) %>%
    pivot_longer(everything(),
      names_to = "day",
      values_to = "value"
    ) %>%
    mutate(
      day = as.integer(str_extract(day, "\\d+")),
      date = data_list$dates[day]
    )

  # Calculate median incidence for observed period
  median_inc_observed <- posterior_inc %>%
    group_by(date, day) %>%
    summarise(
      median_incidence = median(value),
      .groups = "drop"
    ) %>%
    arrange(day) %>%
    mutate(period = "Observed")

  # Extract future predictions (if available)
  future_inc_data <- NULL
  if ("future_incidence" %in% names(rstan::extract(fit))) {
    future_inc <- rstan::extract(fit, pars = "future_incidence")$future_incidence
    
    # Calculate median for future predictions
    future_median <- apply(future_inc, 2, median)
    future_dates <- max(data_list$dates) + seq_len(length(future_median))
    future_inc_data <- data.frame(
      date = future_dates,
      day = max(median_inc_observed$day) + seq_len(length(future_median)),
      median_incidence = future_median,
      period = "Predicted"
    )
    cat(sprintf("Future incidence available for %d days\n", length(future_median)))
  }

  # Combine observed and predicted incidence
  if (!is.null(future_inc_data)) {
    combined_inc <- bind_rows(median_inc_observed, future_inc_data)
  } else {
    combined_inc <- median_inc_observed
    cat("No future predictions available, using observed data only\n")
  }

  # Prepare data for EpiEstim
  rt_data <- data.frame(
    dates = combined_inc$date,
    I = round(combined_inc$median_incidence),
    period = combined_inc$period
  )

  # Ensure positive incidence values (EpiEstim requirement)
  rt_data$I[rt_data$I <= 0] <- 1

  cat(sprintf("Total time series length: %d days\n", nrow(rt_data)))
  cat(sprintf(
    "Observed period: %d days, Predicted period: %d days\n",
    sum(rt_data$period == "Observed"), sum(rt_data$period == "Predicted")
  ))

  # Set time window for Rt estimation
  t_start <- seq(4, nrow(rt_data) - 2)
  t_end <- t_start + 2

  # Configure EpiEstim parameters
  conf <- make_config(
    list(
      mean_si = 3.5, # Mean serial interval for Chikungunya
      std_si = 2.6, # Standard deviation of serial interval
      t_start = t_start,
      t_end = t_end
    )
  )

  # Estimate Rt using parametric serial interval
  tryCatch(
    {
      rt_estimates <- estimate_R(
        rt_data[, c("dates", "I")], # EpiEstim needs only dates and I
        method = "parametric_si",
        config = conf
      )

      # Format results and add period information
      rt_result <- data.frame(
        date = rt_data$dates[rt_estimates$R$t_end],
        mean_rt = round(rt_estimates$R$`Mean(R)`, 4),
        lower_ci = round(rt_estimates$R$`Quantile.0.025(R)`, 4),
        upper_ci = round(rt_estimates$R$`Quantile.0.975(R)`, 4),
        # Use original period information from rt_data
        period = rt_data$period[rt_estimates$R$t_end],
        stringsAsFactors = FALSE
      ) %>%
        mutate(
          ci_string = paste0(mean_rt, " (", lower_ci, ", ", upper_ci, ")"),
          epidemic_phase = case_when(
            lower_ci > 1 ~ "Growth",
            upper_ci < 1 ~ "Decline",
            TRUE ~ "Plateau"
          )
        )

      cat(sprintf("Rt estimation successful for %d time points\n", nrow(rt_result)))
      cat(sprintf(
        "Observed Rt points: %d, Predicted Rt points: %d\n",
        sum(rt_result$period == "Observed"), sum(rt_result$period == "Predicted")
      ))

      # Create comprehensive Rt visualization
      rt_plot <- ggplot(rt_result, aes(x = date)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
        # Different ribbons for observed vs predicted
        geom_ribbon(
          data = filter(rt_result, period == "Observed"),
          aes(ymin = lower_ci, ymax = upper_ci),
          fill = "lightblue", alpha = 0.5
        ) +
        geom_ribbon(
          data = filter(rt_result, period == "Predicted"),
          aes(ymin = lower_ci, ymax = upper_ci),
          fill = "lightgreen", alpha = 0.5
        ) +
        # Different lines for observed vs predicted
        geom_line(
          data = filter(rt_result, period == "Observed"),
          aes(y = mean_rt), color = "blue", size = 1.2
        ) +
        geom_line(
          data = filter(rt_result, period == "Predicted"),
          aes(y = mean_rt), color = "green", size = 1.2
        ) +
        # Points colored by period
        geom_point(aes(y = mean_rt, shape = period, color = period), size = 2.5) +
        scale_color_manual(
          values = c("Observed" = "blue", "Predicted" = "green"),
          name = "Period"
        ) +
        scale_shape_manual(
          values = c("Observed" = 16, "Predicted" = 17),
          name = "Period"
        ) +
        scale_x_date(date_labels = "%m-%d", date_breaks = "2 days") +
        labs(
          title = "Comprehensive Rt Analysis: Observed + Predicted",
          subtitle = "Blue: Observed period Rt; Green: Predicted period Rt",
          x = "Date",
          y = "Effective Reproduction Number (Rt)",
          caption = "Dashed line shows epidemic threshold Rt=1; shaded areas are 95% CI"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "bottom"
        )

      return(list(
        rt_estimates = rt_estimates,
        rt_result = rt_result,
        rt_plot = rt_plot,
        rt_data = rt_data,
        combined_incidence = combined_inc
      ))
    },
    error = function(e) {
      cat("Error in comprehensive Rt estimation:", e$message, "\n")
      cat("Returning NULL for Rt analysis\n")
      return(NULL)
    }
  )
}

# Original function kept for compatibility
calculate_rt_from_model <- function(fit, data_list) {
  cat("\n=== Calculating Rt from Model Output (Observed Period Only) ===\n")

  # This now calls the comprehensive function but filters to observed period
  comprehensive_rt <- calculate_rt_comprehensive(fit, data_list)
  if (!is.null(comprehensive_rt)) {
    # Filter to observed period only for backward compatibility
    observed_rt_result <- comprehensive_rt$rt_result %>%
      filter(period == "Observed")
    
    # Create observed-only plot
    rt_plot <- ggplot(observed_rt_result, aes(x = date)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
        fill = "lightblue", alpha = 0.5
      ) +
      geom_line(aes(y = mean_rt), color = "blue", size = 1.2) +
      geom_point(aes(y = mean_rt, color = epidemic_phase), size = 2.5) +
      scale_color_manual(
        values = c("Growth" = "red", "Decline" = "green", "Plateau" = "orange"),
        name = "Epidemic Phase"
      ) +
      scale_x_date(date_labels = "%m-%d", date_breaks = "2 days") +
      labs(
        title = "Effective Reproduction Number (Rt) Time Series",
        subtitle = "Based on SEIR-SIR Model Predicted Incidence Data (Observed Period)",
        x = "Date",
        y = "Effective Reproduction Number (Rt)",
        caption = "Dashed line shows epidemic threshold Rt=1; shaded areas are 95% confidence intervals"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "bottom"
      )
    return(list(
      rt_estimates = comprehensive_rt$rt_estimates,
      rt_result = observed_rt_result,
      rt_plot = rt_plot,
      rt_data = comprehensive_rt$rt_data[comprehensive_rt$rt_data$period == "Observed", ]
    ))
  } else {
    return(NULL)
  }
}

# Alternative Rt calculation using observed data
calculate_rt_from_observed <- function(data_list) {
  cat("\n=== Calculating Rt from Observed Data ===\n")

  # Prepare observed data for EpiEstim
  rt_data <- data.frame(
    dates = data_list$dates,
    I = data_list$cases
  )
  # Set time window
  t_start <- seq(4, nrow(rt_data) - 2)
  t_end <- t_start + 2
  # Configure EpiEstim
  conf <- make_config(
    list(
      mean_si = 3.5,
      std_si = 2.6,
      t_start = t_start,
      t_end = t_end
    )
  )

  # Estimate Rt
  tryCatch(
    {
      rt_estimates <- estimate_R(
        rt_data,
        method = "parametric_si",
        config = conf
      )

      # Format results
      rt_result <- data.frame(
        date = rt_data$dates[rt_estimates$R$t_end],
        mean_rt = round(rt_estimates$R$`Mean(R)`, 4),
        lower_ci = round(rt_estimates$R$`Quantile.0.025(R)`, 4),
        upper_ci = round(rt_estimates$R$`Quantile.0.975(R)`, 4),
        stringsAsFactors = FALSE
      ) %>%
        mutate(
          ci_string = paste0(mean_rt, " (", lower_ci, ", ", upper_ci, ")"),
          epidemic_phase = case_when(
            lower_ci > 1 ~ "Growth",
            upper_ci < 1 ~ "Decline",
            TRUE ~ "Plateau"
          )
        )
      # Create visualization
      rt_plot <- ggplot(rt_result, aes(x = date)) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
        geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
          fill = "lightcoral", alpha = 0.5
        ) +
        geom_line(aes(y = mean_rt), color = "darkred", size = 1.2) +
        geom_point(aes(y = mean_rt, color = epidemic_phase), size = 2.5) +
        scale_color_manual(
          values = c("Growth" = "red", "Decline" = "green", "Plateau" = "orange"),
          name = "Epidemic Phase"
        ) +
        scale_x_date(date_labels = "%m-%d", date_breaks = "2 days") +
        labs(
          title = "Effective Reproduction Number (Rt) Time Series",
          subtitle = "Based on Observed Case Data",
          x = "Date",
          y = "Effective Reproduction Number (Rt)",
          caption = "Dashed line shows epidemic threshold Rt=1; shaded areas are 95% confidence intervals"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "bottom"
        )
      return(list(
        rt_estimates = rt_estimates,
        rt_result = rt_result,
        rt_plot = rt_plot,
        rt_data = rt_data
      ))
    },
    error = function(e) {
      cat("Error in observed Rt estimation:", e$message, "\n")
      return(NULL)
    }
  )
}

# ===============================================================================
# 8. PARAMETER EXTRACTION FOR REPORTING
# ===============================================================================

# Extract parameter summaries for reporting
extract_parameter_summaries <- function(fit) {
  cat("\n=== Extracting Parameter Summaries for Reporting ===\n")
  
  # Key epidemiological parameters
  key_params <- c("beta", "beta_hv", "sigma_h", "gamma_h", "vie", "beta_vh", "gamma_v", "gamma_rh", "awareness", "mu_v", "phi_incidence", "phi_reported")
  
  # Extract draws
  draws <- as_draws_df(fit)
  
  # Calculate parameter summaries
  param_summaries <- draws %>%
    select(all_of(key_params)) %>%
    summarise(
      across(everything(), list(
        mean = ~ round(mean(.x), 3),
        median = ~ round(median(.x), 3),
        sd = ~ round(sd(.x), 3),
        q025 = ~ round(quantile(.x, 0.025), 3),
        q975 = ~ round(quantile(.x, 0.975), 3)
      ))
    ) %>%
    pivot_longer(everything(), 
                names_to = c("parameter", "statistic"), 
                names_pattern = "(.+)_(.+)") %>%
    pivot_wider(names_from = statistic, values_from = value)
  
  # Calculate derived quantities
  derived_quantities <- param_summaries %>%
    filter(parameter %in% c("sigma_h", "gamma_h")) %>%
    mutate(
      parameter_name = case_when(
        parameter == "sigma_h" ~ "incubation_period",
        parameter == "gamma_h" ~ "infectious_period"
      ),
      mean_days = round(1 / mean, 1),
      q025_days = round(1 / q975, 1),  # Note: inverse relationship
      q975_days = round(1 / q025, 1)   # Note: inverse relationship
    )
  
  # Create formatted parameter descriptions
  param_descriptions <- param_summaries %>%
    mutate(
      parameter_name = case_when(
        parameter == "beta" ~ "Vector Biting Rate (beta)",
        parameter == "beta_hv" ~ "Vector-to-Human Transmission Probability (beta_hv)",
        parameter == "beta_vh" ~ "Human-to-Vector Transmission Probability (beta_vh)",
        parameter == "sigma_h" ~ "Latent Period Rate (sigma_h)",
        parameter == "gamma_h" ~ "Human Recovery Rate (gamma_h)",
        parameter == "vie" ~ "Human-Vector Contact Rate (vie)",
        parameter == "gamma_v" ~ "Vector Recovery Rate (gamma_v)",
        parameter == "gamma_rh" ~ "Reporting Rate (gamma_rh)",
        parameter == "awareness" ~ "Behavioral Response Intensity to Reports (awareness)",
        parameter == "mu_v" ~ "Vector Death Rate (mu_v)",
        parameter == "phi_incidence" ~ "Overdispersion Parameter for Incidence (phi_incidence)",
        parameter == "phi_reported" ~ "Overdispersion Parameter for Reported Cases (phi_reported)"
      ),
      formatted_value = sprintf("%.3f (%.3f - %.3f)", mean, q025, q975),
      unit = case_when(
        parameter %in% c("beta", "sigma_h", "gamma_h", "gamma_v", "gamma_rh", "mu_v") ~ "/day",
        parameter %in% c("beta_hv", "beta_vh", "vie", "awareness") ~ "",
        parameter %in% c("phi_incidence", "phi_reported") ~ ""
      ),
      full_description = paste0(parameter_name, " = ", formatted_value, unit)
    )
  
  # Calculate effective reproduction numbers
  # Extract parameters for Reff calculations
  beta_draws <- draws$beta
  beta_hv_draws <- draws$beta_hv
  beta_vh_draws <- draws$beta_vh
  gamma_h_draws <- draws$gamma_h
  gamma_v_draws <- draws$gamma_v
  vie_draws <- draws$vie
  
  # Calculate Reff (hm): human to mosquito transmission
  # Reff (hm) = (beta * beta_vh * vie) / gamma_v
  reff_hm_draws <- (beta_draws * beta_vh_draws * vie_draws) / gamma_v_draws
  
  # Calculate Reff (mh): mosquito to human transmission  
  # Reff (mh) = (beta * beta_hv * vie) / gamma_h
  reff_mh_draws <- (beta_draws * beta_hv_draws * vie_draws) / gamma_h_draws
  
  # Calculate overall Reff
  reff_total_draws <- reff_hm_draws * reff_mh_draws
  
  # Summarize Reff values
  reff_summaries <- list(
    reff_hm = list(
      mean = round(mean(reff_hm_draws), 2),
      median = round(median(reff_hm_draws), 2),
      q025 = round(quantile(reff_hm_draws, 0.025), 2),
      q975 = round(quantile(reff_hm_draws, 0.975), 2),
      formatted = sprintf("%.2f (%.2f - %.2f)", 
                         mean(reff_hm_draws), 
                         quantile(reff_hm_draws, 0.025), 
                         quantile(reff_hm_draws, 0.975))
    ),
    reff_mh = list(
      mean = round(mean(reff_mh_draws), 2),
      median = round(median(reff_mh_draws), 2),
      q025 = round(quantile(reff_mh_draws, 0.025), 2),
      q975 = round(quantile(reff_mh_draws, 0.975), 2),
      formatted = sprintf("%.2f (%.2f - %.2f)", 
                         mean(reff_mh_draws), 
                         quantile(reff_mh_draws, 0.025), 
                         quantile(reff_mh_draws, 0.975))
    ),
    reff_total = list(
      mean = round(mean(reff_total_draws), 2),
      median = round(median(reff_total_draws), 2),
      q025 = round(quantile(reff_total_draws, 0.025), 2),
      q975 = round(quantile(reff_total_draws, 0.975), 2),
      formatted = sprintf("%.2f (%.2f - %.2f)", 
                         mean(reff_total_draws), 
                         quantile(reff_total_draws, 0.025), 
                         quantile(reff_total_draws, 0.975))
    )
  )
  
  return(list(
    param_summaries = param_summaries,
    derived_quantities = derived_quantities,
    param_descriptions = param_descriptions,
    reff_summaries = reff_summaries
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

  # Step 1: Data preparation
  data_list <- setup_data()

  # Step 2: Model setup
  model_params <- setup_model_parameters()

  # Step 3: Prepare Stan data
  stan_data <- list(
    T = data_list$T,
    ts = data_list$ts,
    cases = data_list$cases,                    # Biological incidence (symptoms data)
    cases_reported = data_list$cases_reported,  # Reported cases 
    N_h = model_params$N_h,
    init_state = model_params$init_state
  )
  # Debug: Print what we're passing to Stan
  cat("\n=== DEBUG: Stan Data Structure ===\n")
  cat("T:", stan_data$T, "\n")
  cat("N_h:", stan_data$N_h, "\n")
  cat("ts length:", length(stan_data$ts), "\n")
  cat("cases length:", length(stan_data$cases), "\n")
  cat("cases_reported length:", length(stan_data$cases_reported), "\n")
  cat("init_state length:", length(stan_data$init_state), "\n")
  cat("init_state:", stan_data$init_state, "\n")
  cat("Sample cases (symptoms):", head(stan_data$cases, 5), "\n")
  cat("Sample cases_reported:", head(stan_data$cases_reported, 5), "\n")

  # Step 4: Run MCMC inference
  fit <- run_mcmc_inference(stan_data, n_chains = 8)

  # Step 4.5: Analyze chain convergence
  convergence_analysis <- analyze_chain_convergence(fit)

  # Step 5: Analyze results
  analysis_results <- analyze_parameters(fit)
  
  # Step 5.5: Extract parameter summaries for reporting
  param_summaries <- extract_parameter_summaries(fit)

  # Step 6: Create visualizations
  plots <- create_model_plots(fit, data_list)

  # Step 7: Model diagnostics
  diagnostics <- model_diagnostics(fit, data_list)

  # Step 8: Calculate Rt from multiple approaches
  rt_from_model <- calculate_rt_from_model(fit, data_list) # Observed period only
  rt_from_observed <- calculate_rt_from_observed(data_list) # Observed data only
  rt_comprehensive <- calculate_rt_comprehensive(fit, data_list) # Observed + Predicted

  # Display main plots
  print(plots$prediction_plot)
  if (!is.null(plots$future_prediction_plot)) {
    print(plots$future_prediction_plot)
  }
  if (!is.null(plots$future_reported_plot)) {
    print(plots$future_reported_plot)
  }
  # Display Rt plots if available
  if (!is.null(rt_from_model)) {
    print(rt_from_model$rt_plot)
  }
  if (!is.null(rt_from_observed)) {
    print(rt_from_observed$rt_plot)
  }
  if (!is.null(rt_comprehensive)) {
    print(rt_comprehensive$rt_plot)
  }

  cat("\n===============================================================================\n")
  cat("Analysis completed successfully!\n")
  cat("===============================================================================\n")

  # Save results to file
  results <- list(
    fit = fit,
    data = data_list,
    results = analysis_results,
    param_summaries = param_summaries,
    plots = plots,
    diagnostics = diagnostics,
    convergence = convergence_analysis,
    rt_from_model = rt_from_model,
    rt_from_observed = rt_from_observed,
    rt_comprehensive = rt_comprehensive
  )
  
  # Save results to file
  saveRDS(results, "fitting/model_results.rds")
  cat("Results saved to fitting/model_results.rds\n")
  
  return(results)
}
