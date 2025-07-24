library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 1. 读取数据

dates <- as.Date(c("2025-07-08","2025-07-09","2025-07-10","2025-07-11",
                   "2025-07-12","2025-07-13","2025-07-14","2025-07-15",
                   "2025-07-16","2025-07-17","2025-07-18","2025-07-19",
                   "2025-07-20","2025-07-21","2025-07-22"))
cases <- c(1,1,3,8,20,48,116,280,104,296,393,473,540,373,536)
T <- length(cases)
ts <- as.numeric(dates - dates[1]) + 1

# 2. 初始状态假设
N_h <- 969.89e4
init_state <- c(
  N_h - 1,  # S_h0
  0,        # E_h0
  1,        # I_h0
  0,        # R_h0
  1e5,      # S_v0 （可调）
  0,        # I_v0
  0         # R_v0
)

stan_data <- list(
  T          = T,
  ts         = ts,
  cases      = cases,
  N_h        = N_h,
  init_state = init_state
)

# 3. 编译+采样
fit <- stan(
  file    = "fitting/seir_zxc.stan",
  data    = stan_data,
  iter    = 2000,
  chains  = 4,
  control = list(adapt_delta = 0.9)
)

# 4. 查看结果
print(fit, pars = c("beta","beta_hv","sigma_h","gamma_h","vie",
                    "beta_vh","gamma_v","phi"))

# 安装并加载必要包
# install.packages(c("tidybayes","ggplot2","dplyr"))
library(tidybayes)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(dplyr)
library(posterior)
library(stringr)

# 从 fit 提取后验 incidence
posterior_inc <- fit %>%
  as_draws_df() %>% 
  select(starts_with("incidence")) %>% 
  pivot_longer(everything(),
               names_to = "day",
               values_to = "value") %>%
  mutate(day = as.integer(str_remove(day, "incidence\\[")),
         date = dates[day])

# 计算中位数 & 95% 区间
summ_inc <- posterior_inc %>%
  group_by(date) %>%
  summarise(
    median = median(value),
    lower  = quantile(value, 0.025),
    upper  = quantile(value, 0.975)
  )

# 作图
ggplot() +
  geom_ribbon(data = summ_inc,
              aes(x = date, ymin = lower, ymax = upper),
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = summ_inc,
            aes(x = date, y = median),
            color = "steelblue", linewidth = 1) +
  geom_point(data = tibble(date = dates, cases = cases),
             aes(x = date, y = cases),
             color = "black", size = 1.5) +
  labs(
    title = "SEIR–SIR 模型后验预测与观测比较",
    x = "日期", y = "每日新增病例数",
    caption = "中位数（线）与 95% 后验区间（带）"
  ) +
  theme_minimal(base_family = "STKaiti")

# 提取所有参数
draws <- as_draws_df(fit)


