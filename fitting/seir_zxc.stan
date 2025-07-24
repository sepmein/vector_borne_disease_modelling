// seir_sir_model.stan
functions {
  // 7-compartment SEIR-SIR model
  real[] seir_sir_ode(real t,
                      real[] y,
                      real[] theta,
                      real[] x_r,
                      int[]   x_i) {
    real S_h = y[1];
    real E_h = y[2];
    real I_h = y[3];
    real R_h = y[4];
    real S_v = y[5];
    real I_v = y[6];
    real R_v = y[7];

    real N_h = S_h + E_h + I_h + R_h;
    real N_v = S_v + I_v + R_v;

    real beta    = theta[1]; // mosquito biting rate
    real beta_hv = theta[2]; // transmission probability: mosquito to human
    real sigma_h = theta[3]; // exposed to infectious rate
    real gamma_h = theta[4]; // recovery rate for humans
    real vie     = theta[5]; // human to mosquito contact rate
    real beta_vh = theta[6]; // transmission probability: human to mosquito
    real gamma_v = theta[7]; // recovery rate for mosquitoes

    real dS_h = -beta * beta_hv * S_h * I_v / N_v;
    real dE_h =  beta * beta_hv * S_h * I_v / N_v - sigma_h * E_h;
    real dI_h =                sigma_h * E_h - gamma_h * I_h;
    real dR_h =                               gamma_h * I_h;

    real dS_v = -vie  * beta_vh * S_v * I_h / N_h;
    real dI_v =  vie  * beta_vh * S_v * I_h / N_h - gamma_v * I_v;
    real dR_v =                                gamma_v * I_v;

    return {dS_h, dE_h, dI_h, dR_h, dS_v, dI_v, dR_v};
  }
}

data {
  int<lower=1> T;             // number of days
  real        ts[T];          // observation time points (starting from 0)
  int<lower=0> cases[T];      // daily new cases
  real<lower=0> N_h;          // total human population
  real<lower=0> init_state[7]; // initial state: S_h, E_h, I_h, R_h, S_v, I_v, R_v
}

parameters {
  real<lower=0> beta;
  real<lower=0> beta_hv;
  real<lower=0> sigma_h;
  real<lower=0> gamma_h;
  real<lower=0> vie;
  real<lower=0> beta_vh;
  real<lower=0> gamma_v;
  real<lower=0> phi;  // overdispersion (negative binomial distribution)
}

transformed parameters {
  real y_hat[T,7];
  real incidence[T];

  // ODE solver
  y_hat = integrate_ode_rk45(
    seir_sir_ode,
    init_state, 
    0.0, 
    ts, 
    {beta, beta_hv, sigma_h, gamma_h, vie, beta_vh, gamma_v},
    rep_array(0.0,0),   /* x_r */
    rep_array(0,0)      /* x_i */
  );

  // daily incidence approximated as sigma_h * E_h
  for (t in 1:T)
    incidence[t] = sigma_h * y_hat[t,2];
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
  // beta    ~ normal(0.19, 0.39);
  // beta_hv ~ normal(0.001, 0.54);
  // sigma_h ~ normal(0.25, 0.5);
  // gamma_h ~ normal(1.0/8.0, 1.0/2.0);
  // vie     ~ normal(0.005, 0.35);
  // beta_vh ~ normal(0.19, 0.39);
  // gamma_v ~ normal(1.0/6.0, 1.0/2.0);
  // phi     ~ cauchy(0,5);
  // BETTER PRIORS using appropriate distributions for positive rates:
  beta    ~ lognormal(log(0.2), 0.5);   // Mosquito biting rate (median=0.2, can't be negative)
  beta_hv ~ beta(0.24, 0.2);                 // Transmission prob (mean=0.2, bounded 0-1)
  sigma_h ~ gamma(1/3, (1/3 - 1/4));               // Incubation rate (mean=1/5 day, 5-day incubation)
  gamma_h ~ gamma(1/4, (1/4 - 1/8));               // Recovery rate (mean=1/7 day, 7-day infectious)
  vie     ~ lognormal(log(0.23), 0.3);  // Contact rate (median=0.05, can't be negative)
  beta_vh ~ beta(0.24, 0.24);                 // Transmission prob (mean=0.3, bounded 0-1)
  gamma_v ~ gamma(1/3.5, (1/3.5 - 1/6));               // Mosquito recovery (mean=1/10 day, 10-day infectious)
  phi     ~ exponential(0.1);           // Overdispersion (more realistic than Cauchy)

  // Likelihood: negative binomial distribution for observed daily cases
  for (t in 1:T)
    cases[t] ~ neg_binomial_2(incidence[t], phi);
}

generated quantities {
  // Posterior predictions can be added here
}
