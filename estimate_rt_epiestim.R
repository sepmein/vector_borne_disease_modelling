# estimate_rt_epiestim.R
# 用EpiEstim估算实时再生数Rt
# Author: 
#
# Usage:
#   estimate_rt_epiestim(I_new, mean_si, std_si)
#   # I_new: numeric vector of daily new cases
#   # mean_si: mean of serial interval (days)
#   # std_si: standard deviation of serial interval (days)
#
# 输出：EpiEstim结果对象，并自动绘图

library(EpiEstim)

estimate_rt_epiestim <- function(I_new, mean_si, std_si) {
  res <- estimate_R(
    incid = I_new,
    method = "parametric_si",
    config = make_config(list(
      mean_si = mean_si,
      std_si = std_si
    ))
  )
  plot(res)
  return(res)
}

# 示例 Example:
# I_new <- c(2, 3, 5, 7, 10, 15, 20, 18, 16, 14, 12, 10, 8, 6, 5, 4, 3, 2, 1, 0)
# estimate_rt_epiestim(I_new, mean_si = 6, std_si = 2) 