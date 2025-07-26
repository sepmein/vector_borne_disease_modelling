# model_seir_sei.R
# SEIR-sei动力学模型定义与仿真
# Author: 
#
# This module defines the SEIR model for humans and sei model for vectors (mosquitoes),
# and provides a function to run the simulation.
#
# 人群：SEIR结构
# 蚊群：sei结构
#
# 用法 Usage:
#   result <- run_seir_sei_model(params, init, times)
#   # params: named vector of parameters
#   # init: named vector of initial states
#   # times: time vector
#
# 返回值：data.frame，包含各舱室人数随时间变化

library(deSolve)

# definition	value
# mosquito birth rate	83.75
# natural death rate for mosquito	0.03
# mosquito biting rate for transfer of infection from infected human to susceptible mosquito	0.30
# transmission probability per contact in susceptible mosquitoes	0.30
# progression rate from exposed to infectious mosquito population	0.29
# "mosquito biting rate for transfer of infection from
# infectious mosquito to susceptible human"	0.30
# transmission probability per contact in susceptible human	0.30
# progression rate of exposed to infected human population	0.25
# progression rate of infected to recovered human	0.16


# ODE function
dengue_seir_sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N_h <- S_h + E_h + I_h + R_h
    N_v <- S_v + I_v + R_v
    dS_h <- -beta * beta_hv * S_h * I_v / N_v
    dE_h <- beta * beta_hv * S_h * I_v / N_v - sigma_h * E_h
    dI_h <- sigma_h * E_h - gamma_h * I_h
    dR_h <- gamma_h * I_h
    dS_v <- -vie * beta_vh * S_v * I_h / N_h
    dI_v <- vie * beta_vh * S_v * I_h / N_h - gamma_v * I_v
    dR_v <- gamma_v * I_v
    return(list(c(dS_h, dE_h, dI_h, dR_h, dS_v, dI_v, dR_v)))
  })
}

# 初始化参数
# 人口总数，佛山市2024年人口数
N_h <- 969.89 * 10000

parameters <- c(
  # beta, mosquito biting rate
  beta = 0.25,
  # beta_hv, 易感染人群被蚊虫叮咬后感染的速率
  beta_hv = 0.24,
  # sigma, from exposed to infectious human
  sigma_h = 0.25,
  # gamma_h, from infectious to recovered human
  gamma_h = 0.16,
  # vie, Mosquito biting rate for transfer of infection from infected human population(I_h) to susceptible mosquito population(S_v)
  vie = 0.25,
  # beta_vh, from infectious human to susceptible mosquito
  beta_vh = 0.24,
  # gamma_v, from exposed mosquito to infectious mosquito
  gamma_v = 1 / 3.5

)

# 时间序列
times <- seq(0, 365, by = 1)

# 仿真主函数
run_seir_sir_model <- function(params, init, times) {
  out <- ode(y = init, times = times, func = dengue_seir_sir, parms = params)
  as.data.frame(out)
}

