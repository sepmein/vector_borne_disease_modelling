# model_seir_sir.R
# SEIR-SIR动力学模型定义与仿真
# Author: 
#
# This module defines the SEIR model for humans and SIR model for vectors (mosquitoes),
# and provides a function to run the simulation.
#
# 人群：SEIR结构
# 蚊群：SIR结构
#
# 用法 Usage:
#   result <- run_seir_sir_model(params, init, times)
#   # params: named vector of parameters
#   # init: named vector of initial states
#   # times: time vector
#
# 返回值：data.frame，包含各舱室人数随时间变化

library(deSolve)

# ODE function
dengue_seir_sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N_h <- S_h + E_h + I_h + R_h
    N_v <- S_v + I_v + R_v
    dS_h <- -beta_hv * S_h * I_v / N_v
    dE_h <- beta_hv * S_h * I_v / N_v - sigma_h * E_h
    dI_h <- sigma_h * E_h - gamma_h * I_h
    dR_h <- gamma_h * I_h
    dS_v <- -beta_vh * S_v * I_h / N_h
    dI_v <- beta_vh * S_v * I_h / N_h - gamma_v * I_v
    dR_v <- gamma_v * I_v
    return(list(c(dS_h, dE_h, dI_h, dR_h, dS_v, dI_v, dR_v)))
  })
}

# 仿真主函数
run_seir_sir_model <- function(params, init, times) {
  out <- ode(y = init, times = times, func = dengue_seir_sir, parms = params)
  as.data.frame(out)
} 