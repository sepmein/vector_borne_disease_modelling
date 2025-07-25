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
  N_h <- 969.89e4 # Total human population (Foshan city)

  # More realistic initial conditions for dengue outbreak
  # Based on early case detection and vector dynamics
  init_state <- c(
    N_h - 10, # S_h0: Initially susceptible humans (assume some early exposure)
    5, # E_h0: Initially exposed humans (realistic for dengue)
    5, # I_h0: Initially infectious humans (matches early cases)
    0, # R_h0: Initially recovered humans
    5000, # S_v0: Much smaller vector population (more realistic)
    10, # I_v0: Some initially infectious vectors
    0 # R_v0: Initially recovered vectors
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
run_mcmc_inference <- function(stan_data, n_iter = 2000, n_chains = 8,
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
  key_params <- c("beta", "beta_hv", "sigma_h", "gamma_h", "vie", "beta_vh", "gamma_v")

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

  # Future prediction extraction and plot
  if ("future_incidence" %in% names(rstan::extract(fit))) {
    future_inc <- rstan::extract(fit, pars = "future_incidence")$future_incidence
    future_pred <- as.data.frame(t(apply(future_inc, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
    colnames(future_pred) <- c("lower_95", "lower_50", "median", "upper_50", "upper_95")
    future_pred$day <- 1:nrow(future_pred)
    future_pred$date <- max(data_list$dates) + future_pred$day

    p_future <- ggplot(future_pred, aes(x = date)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "green", alpha = 0.2) +
      geom_ribbon(aes(ymin = lower_50, ymax = upper_50), fill = "green", alpha = 0.4) +
      geom_line(aes(y = median), color = "green") +
      labs(title = "Future Prediction (Next 7 Days)", y = "Predicted Daily Cases", x = "Date") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)
      )
  } else {
    p_future <- NULL
    future_pred <- NULL
  }

  # Combine model fit and future prediction for plotting
  model_pred <- summ_inc %>%
    select(date, median, lower_50, upper_50, lower_95, upper_95) %>%
    mutate(type = "Model Fit")

  future_pred2 <- future_pred %>%
    select(date, median, lower_50, upper_50, lower_95, upper_95) %>%
    mutate(type = "Future Prediction")

  all_pred <- bind_rows(model_pred, future_pred2)

  return(list(
    prediction_plot = p1,
    parameter_plot = p2,
    summary_data = summ_inc,
    future_prediction_plot = p_future,
    future_prediction_summary = future_pred,
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
      scale_color_manual(values = c("Model Fit" = "steelblue", "Future Prediction" = "green")) +
      scale_fill_manual(values = c("Model Fit" = "steelblue", "Future Prediction" = "green")) +
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

  # Posterior predictive checks
  y_rep <- rstan::extract(fit, pars = "incidence")$incidence
  y_obs <- data_list$cases

  # Calculate Bayesian R-squared
  # (simplified version - could be enhanced)
  residuals <- sweep(y_rep, 2, y_obs, "-")
  ss_res <- rowSums(residuals^2)
  ss_tot <- sum((y_obs - mean(y_obs))^2)
  r_squared <- 1 - ss_res / ss_tot

  cat(sprintf(
    "Bayesian R-squared: %.3f [%.3f, %.3f]\n",
    median(r_squared), quantile(r_squared, 0.025), quantile(r_squared, 0.975)
  ))

  # Root Mean Square Error
  rmse <- sqrt(rowMeans(residuals^2))
  cat(sprintf(
    "RMSE: %.2f [%.2f, %.2f]\n",
    median(rmse), quantile(rmse, 0.025), quantile(rmse, 0.975)
  ))

  # Mean Absolute Error
  mae <- rowMeans(abs(residuals))
  cat(sprintf(
    "MAE: %.2f [%.2f, %.2f]\n",
    median(mae), quantile(mae, 0.025), quantile(mae, 0.975)
  ))

  return(list(
    r_squared = r_squared,
    rmse = rmse,
    mae = mae
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
        stringsAsFactors = FALSE
      ) %>%
        mutate(
          ci_string = paste0(mean_rt, " (", lower_ci, ", ", upper_ci, ")"),
          epidemic_phase = case_when(
            lower_ci > 1 ~ "Growth",
            upper_ci < 1 ~ "Decline",
            TRUE ~ "Plateau"
          ),
          # Determine period based on date
          period = ifelse(date <= max(data_list$dates), "Observed", "Predicted")
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
    cases = data_list$cases,
    N_h = model_params$N_h,
    init_state = model_params$init_state
  )
  # Debug: Print what we're passing to Stan
  cat("\n=== DEBUG: Stan Data Structure ===\n")
  cat("T:", stan_data$T, "\n")
  cat("N_h:", stan_data$N_h, "\n")
  cat("ts length:", length(stan_data$ts), "\n")
  cat("cases length:", length(stan_data$cases), "\n")
  cat("init_state length:", length(stan_data$init_state), "\n")
  cat("init_state:", stan_data$init_state, "\n")

  # Step 4: Run MCMC inference
  fit <- run_mcmc_inference(stan_data, n_chains = 8)

  # Step 4.5: Analyze chain convergence
  convergence_analysis <- analyze_chain_convergence(fit)

  # Step 5: Analyze results
  analysis_results <- analyze_parameters(fit)

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

  return(list(
    fit = fit,
    data = data_list,
    results = analysis_results,
    plots = plots,
    diagnostics = diagnostics,
    convergence = convergence_analysis,
    rt_from_model = rt_from_model,
    rt_from_observed = rt_from_observed,
    rt_comprehensive = rt_comprehensive
  ))
}