// seir_SEI_model.stan
functions {
  // 7-compartment SEIR-SEI model
  real[] seir_sir_ode(real t,
                      real[] y,
                      real[] theta,
                      real[] x_r,
                      int[]   x_i) {
    real S_h = y[1];
    real E_h = y[2];
    real I_h = y[3];
    real R_h = y[4];
    real Report_h = y[5];
    real S_v = y[6];
    real I_v = y[7];
    real R_v = y[8];

    real N_h = S_h + E_h + I_h + R_h;
    real N_v = S_v + I_v + R_v;

    real beta    = theta[1]; // mosquito biting rate
    real beta_hv = theta[2]; // transmission probability: mosquito to human
    real sigma_h = theta[3]; // exposed to infectious rate
    real gamma_h = theta[4]; // recovery rate for humans
    real vie     = theta[5]; // human to mosquito contact rate
    real beta_vh = theta[6]; // transmission probability: human to mosquito
    real gamma_v = theta[7]; // recovery rate for mosquitoes
    real gamma_rh = theta[8]; // reporting rate for infected humans
    real awareness = theta[9]; // awareness parameter: higher values = stronger behavioral response to reported cases
    real mu_v = theta[10]; // vector death rate (mortality)

    // Calculate awareness effect based on recent reporting intensity
    real recent_reports = fmax(1.0, Report_h / (t + 1.0)); // Average daily reports since start
    real awareness_effect = 1.0 / (1.0 + awareness * recent_reports);
    
    real dS_h = -beta * beta_hv * S_h * I_v / N_v * awareness_effect;
    real dE_h =  beta * beta_hv * S_h * I_v / N_v * awareness_effect - sigma_h * E_h;
    real dI_h =                sigma_h * E_h - gamma_h * I_h;
    real dR_h =                               gamma_h * I_h;
    real dReport_h = gamma_rh * I_h;

    // Vector birth rate maintains population size (birth_rate = death_rate for stability)
    real dS_v = mu_v * N_v - vie * beta_vh * S_v * I_h / N_h * awareness_effect - mu_v * S_v;
    real dI_v = vie * beta_vh * S_v * I_h / N_h * awareness_effect - gamma_v * I_v - mu_v * I_v;
    real dR_v = gamma_v * I_v - mu_v * R_v;

    return {dS_h, dE_h, dI_h, dR_h, dReport_h, dS_v, dI_v, dR_v};
  }
}

data {
  int<lower=1> T;             // number of days
  real        ts[T];          // observation time points (starting from 0)
  int<lower=0> cases[T];      // daily new cases (biological incidence)
  int<lower=0> cases_reported[T]; // daily new cases reported (observed data)
  real<lower=0> N_h;          // total human population
  real<lower=0> init_state[8]; // initial state: S_h, E_h, I_h, R_h, Report_h, S_v, I_v, R_v
}

parameters {
  real<lower=0> beta;
  real<lower=0> beta_hv;
  real<lower=0> sigma_h;
  real<lower=0> gamma_h;
  real<lower=0> vie;
  real<lower=0> beta_vh;
  real<lower=0> gamma_v;
  real<lower=0> gamma_rh;
  real<lower=0> awareness;
  real<lower=0> mu_v;          // vector death rate
  real<lower=0> phi_incidence;  // overdispersion for biological incidence
  real<lower=0> phi_reported;   // overdispersion for reported cases
}

transformed parameters {
  real y_hat[T,8];
  real incidence[T];           // biological incidence (E to I transition)
  real reported_pred[T];       // predicted reported cases

  // ODE solver
  y_hat = integrate_ode_rk45(
    seir_sir_ode,
    init_state, 
    0.0, 
    ts, 
    {beta, beta_hv, sigma_h, gamma_h, vie, beta_vh, gamma_v, gamma_rh, awareness, mu_v},
    rep_array(0.0,0),   /* x_r */
    rep_array(0,0)      /* x_i */
    );

  // daily incidence approximated as sigma_h * E_h (with positivity constraint)
  for (t in 1:T)
    incidence[t] = fmax(sigma_h * y_hat[t,2], 1e-10);

  // daily reported cases approximated as gamma_rh * I_h (with positivity constraint)
  for (t in 1:T)
    reported_pred[t] = fmax(gamma_rh * y_hat[t,3], 1e-10);
}

model {
  // Priors - CORRECTED from 5%-95% CI to proper mean/sd
  // Original intentions (keeping as comments):
  // beta: 5%-95% CI [?, ?] - mosquito biting rate
  // beta_hv: 5%-95% CI [?, ?] - mosquito to human transmission probability  
  // sigma_h: 5%-95% CI [?, ?] - human incubation rate (1/incubation_days)
  // gamma_h: 5%-95% CI [?, ?] - human recovery rate (1/infectious_days)
  // vie: 5%-95% CI [?, ?] - human-mosquito contact rate
  // beta_vh: 5%-95% CI [?, ?] - human to mosquito transmission probability
  // gamma_v: 5%-95% CI [?, ?] - mosquito recovery rate (1/infectious_days)
  
  // Priors 95% CI
  // beta    ~ normal(0.19, 0.79);
  // beta_hv ~ normal(0.001, 0.54);
  // sigma_h ~ normal(0.25, 0.8);
  // gamma_h ~ normal(1.0/8.0, 1.0/2.0);
  // vie     ~ normal(0.005, 0.75);
  // beta_vh ~ normal(0.19, 0.59);
  // gamma_v ~ normal(1.0/6.0, 1.0/2.0);
  // phi     ~ cauchy(0,5);
  // BETTER PRIORS using appropriate distributions for positive rates:
  beta    ~ normal(0.2, 0.5);   
  beta_hv ~ normal(0.24, 0.2);  
  sigma_h ~ normal(1/3, 0.08333);               
  gamma_h ~ normal(1/4, 0.125);               
  vie     ~ normal(0.23, 0.3);  
  beta_vh ~ normal(0.24, 0.24);                 
  gamma_v ~ normal(1/3.5, 0.119);               
  gamma_rh ~ normal(1/3.5, 0.5);               
  awareness ~ exponential(5);  // Positive values only, mean = 0.2, allows range 0.01-1.0
  mu_v ~ normal(1/10, 0.05);   // Vector death rate: mean 10-day lifespan, typical for mosquitoes
  phi_incidence ~ exponential(0.1);    // overdispersion for biological process        
  phi_reported ~ exponential(0.2);     // overdispersion for reporting process

  // Likelihood: negative binomial distribution for observed biological incidence
  for (t in 1:T)
    cases[t] ~ neg_binomial_2(incidence[t], phi_incidence);
  
  // Likelihood: negative binomial distribution for observed reported cases
  for (t in 1:T)
    cases_reported[t] ~ neg_binomial_2(reported_pred[t], phi_reported);
}

generated quantities {
  real future_incidence[7];
  real future_reported[7];     // ADDED: Missing declaration
  {
    real y_future[7, 8];       // FIXED: Changed from [7,7] to [7,8] for 8 state variables
    real future_ts[7];
    for (i in 1:7)
      future_ts[i] = ts[T] + i;
    y_future = integrate_ode_rk45(
      seir_sir_ode,
      y_hat[T], // state at last observed day
      ts[T],    // time at last observed day
      future_ts, // next 7 days
      {beta, beta_hv, sigma_h, gamma_h, vie, beta_vh, gamma_v, gamma_rh, awareness, mu_v}, // FIXED: Added gamma_rh, awareness, and mu_v
      rep_array(0.0,0),
      rep_array(0,0)
    );
    for (i in 1:7)
      future_incidence[i] = fmax(sigma_h * y_future[i,2], 1e-10);
    for (i in 1:7)
      future_reported[i] = fmax(gamma_rh * y_future[i,3], 1e-10);
  }
}