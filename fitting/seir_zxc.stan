// seir_sir_model.stan
functions {
  // 7 室 SEIR-SIR 模型
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

    real beta    = theta[1]; // 蚊子叮人
    real beta_hv = theta[2]; // 蚊→人 传染概率
    real sigma_h = theta[3]; // 潜→传
    real gamma_h = theta[4]; // 传→愈
    real vie     = theta[5]; // 人→蚊 叮咬
    real beta_vh = theta[6]; // 人→蚊 传染概率
    real gamma_v = theta[7]; // 蚊 感→愈

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
  int<lower=1> T;             // 天数
  real        ts[T];          // 观测时间点（从 0 开始）
  int<lower=0> cases[T];      // 每日新增病例
  real<lower=0> N_h;          // 人群总量
  real<lower=0> init_state[7]; // 初始状态：S_h, E_h, I_h, R_h, S_v, I_v, R_v
}

parameters {
  real<lower=0> beta;
  real<lower=0> beta_hv;
  real<lower=0> sigma_h;
  real<lower=0> gamma_h;
  real<lower=0> vie;
  real<lower=0> beta_vh;
  real<lower=0> gamma_v;
  real<lower=0> phi;  // overdispersion（负二项分布）
}

transformed parameters {
  real y_hat[T,7];
  real incidence[T];

  // ODE 求解
  y_hat = integrate_ode_rk45(
    seir_sir_ode,
    init_state, 
    0.0, 
    ts, 
    {beta, beta_hv, sigma_h, gamma_h, vie, beta_vh, gamma_v},
    rep_array(0.0,0),   /* x_r */
    rep_array(0,0)      /* x_i */
  );

  // 每日新增 ≈ σ_h * E_h
  for (t in 1:T)
    incidence[t] = sigma_h * y_hat[t,2];
}

model {
  // ——先验（可根据实际调整）——
  beta    ~ normal(0.19, 0.39);
  beta_hv ~ normal(0.001, 0.54);
  sigma_h ~ normal(0.25, 0.5);
  gamma_h ~ normal(1/8, 1/2);
  vie     ~ normal(0.005, 0.35);
  beta_vh ~ normal(0.19, 0.39);
  gamma_v ~ normal(1/6, 1/2);
  phi     ~ cauchy(0,5);

  // ——似然：负二项分布拟合观测日增——
  for (t in 1:T)
    cases[t] ~ neg_binomial_2(incidence[t], phi);
}

generated quantities {
  // 这里可加后验预测等
}