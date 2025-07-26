# seir_dengue_model.R
# 主脚本：矢量传播疾病动力学建模与Rt估算
# Main script: Vector-borne disease modelling and Rt estimation
#
# 1. 拟合潜伏期分布
# 2. 运行SEIR-SEI动力学模型
# 3. 估算实时再生数Rt

# 加载模块
source('fit_incubation_distribution.R')
source('model_seir_sir.R')
source('estimate_rt_epiestim.R')
library(ggplot2)
library(reshape2)

# 1. 拟合潜伏期分布（如有原始数据）
# incubation_data <- c(4, 5, 6, 7, 8, 5, 6, 7, 4, 8, 9, 6, 5, 7, 8, 6, 5, 7, 8, 6)
# fit_result <- fit_incubation_gamma(incubation_data)
# print(fit_result)
# mean_incubation <- fit_result$mean
# sd_incubation <- fit_result$sd

# 2. 运行SEIR-SEI模型
params <- c(
  beta_hv = 0.5,    # vector->human
  beta_vh = 0.5,    # human->vector
  sigma_h = 1/6,    # incubation rate (可用拟合均值)
  gamma_h = 1/7,    # recovery rate (human)
  gamma_v = 0       # recovery rate (vector)
)
init <- c(
  S_h = 990000,
  E_h = 5000,
  I_h = 5000,
  R_h = 0,
  S_v = 100000,
  I_v = 1000,
  R_v = 0
)
times <- seq(0, 180, by = 1)
out <- run_seir_sir_model(params, init, times)

# 画图：人群和蚊群
out_long_h <- melt(out[, c('time', 'S_h', 'E_h', 'I_h', 'R_h')], id = 'time')
out_long_v <- melt(out[, c('time', 'S_v', 'I_v', 'R_v')], id = 'time')
p1 <- ggplot(out_long_h, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  labs(title = 'SEIR Model for Humans (Dengue)',
       x = 'Time (days)',
       y = 'Number of Humans',
       color = 'Compartment') +
  theme_minimal()
p2 <- ggplot(out_long_v, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  labs(title = 'SIR Model for Vectors (Mosquitoes)',
       x = 'Time (days)',
       y = 'Number of Vectors',
       color = 'Compartment') +
  theme_minimal()
print(p1)
print(p2)

# 3. 估算Rt（用模型输出的每日新发病例数）
sigma_h <- params['sigma_h']
I_new <- sigma_h * out$E_h  # 每日新发病例数
mean_si <- 6  # 代间期均值
std_si <- 2   # 代间期标准差
estimate_rt_epiestim(I_new, mean_si, std_si)
