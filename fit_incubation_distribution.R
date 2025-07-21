# fit_incubation_distribution.R
# 潜伏期分布拟合模块
# Author: 
#
# 用于拟合潜伏期的统计分布（如伽马分布），并输出参数。
#
# Usage:
#   result <- fit_incubation_gamma(data)
#   # data: numeric vector of observed incubation periods
#   # result: list with shape, rate, mean, sd

library(MASS)

fit_incubation_gamma <- function(data) {
  fit <- fitdistr(data, densfun = "gamma")
  shape <- fit$estimate["shape"]
  rate <- fit$estimate["rate"]
  mean_incubation <- shape / rate
  sd_incubation <- sqrt(shape) / rate
  list(shape = shape, rate = rate, mean = mean_incubation, sd = sd_incubation)
}

# 示例 Example:
# data <- c(4, 5, 6, 7, 8, 5, 6, 7, 4, 8, 9, 6, 5, 7, 8, 6, 5, 7, 8, 6)
# result <- fit_incubation_gamma(data)
# print(result) 