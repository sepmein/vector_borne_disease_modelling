if (!require(minpack.lm)) install.packages("minpack.lm")
if (!require(broom)) install.packages("broom")
library(minpack.lm)
library(broom)        # tidy, glance, augment
library(ggplot2)
library(dplyr)
library(showtext)
showtext_auto()

# import raw cases_reported_foshan.csv
cases_reported <- fread("data/cases_reported_foshan.csv")

# process cases_reported_foshan
cases_reported_foshan <- cases_reported[district == "佛山市"]

# Create a proper date column from year, month, and day
cases_reported_foshan[, date := as.Date(paste(year, month, day, sep = "-"))]

#    district  year month   day accumulated_cases accumulated_recover       date
#      <char> <int> <int> <int>             <int>               <int>     <Date>
# 1:   佛山市  2025     7    21              2658                  NA 2025-07-21
# 2:   佛山市  2025     7    20              2285                  NA 2025-07-20
# 3:   佛山市  2025     7    19              1873                 720 2025-07-19
# 4:   佛山市  2025     7    18              1199                  NA 2025-07-18
# 5:   佛山市  2025     7    15               478                  NA 2025-07-15

# plot cases_reported_foshan, bar chart, x axis is date, y axis is cases_reported
ggplot(cases_reported_foshan, aes(x = date, y = accumulated_cases)) +
  geom_bar(stat = "identity") +
  labs(title = "Cases Reported in Foshan", x = "Date", y = "Cumulative Cases Reported") +
  # 增加数字标签
  geom_text(aes(label = accumulated_cases), vjust = -0.5)

# 设置疫情起始日期和观测区间
start_date <- as.Date("2025-06-20")
first_obs_date <- as.Date("2025-07-15")
end_date <- max(cases_reported_foshan$date)

# 生成完整日期序列
dates_full <- data.frame(date = seq(start_date, end_date, by = "day"))

# 合并原始累计病例数据
cases_full <- left_join(dates_full, cases_reported_foshan[, c("date", "accumulated_cases")], by = "date")

# 只用6月20日和7月15日两点拟合指数模型
start_cases <- 1  # 假定6月20日累计病例数为1
end_cases <- cases_full$accumulated_cases[cases_full$date == first_obs_date]
fit_data <- data.frame(
    date = c(start_date, first_obs_date),
    t = c(0, as.numeric(first_obs_date - start_date)),
    accumulated_cases = c(start_cases, end_cases)
)
model <- nlsLM(
  accumulated_cases ~ a * exp(b * t),
  data = fit_data,
  start = list(a = fit_data$accumulated_cases[1], b = 0.1),
  control = nls.lm.control(maxiter = 500)
)

# 用模型外推6月20日到7月14日的累计病例数
cases_full <- cases_full %>%
  mutate(t = as.numeric(date - start_date))
missing_idx <- which(cases_full$date < first_obs_date)
cases_full$accumulated_cases_fitted <- NA
cases_full$accumulated_cases_fitted[missing_idx] <- predict(model, newdata = cases_full[missing_idx, ])
cases_full$accumulated_cases[missing_idx] <- cases_full$accumulated_cases_fitted[missing_idx]

# 计算每日新增病例数
cases_full <- cases_full %>%
  arrange(date) %>%
  mutate(daily_cases = accumulated_cases - lag(accumulated_cases, default = 0))

# 随机把数据取整
cases_full$daily_cases <- round(cases_full$daily_cases)

# 可视化
ggplot(cases_full, aes(x = date, y = daily_cases)) +
  geom_col(fill = "skyblue") +
  labs(title = "估算每日发病数（缺失区间指数外推）", x = "日期", y = "每日新增病例数") +
  # 增加数字标签
  geom_text(aes(label = daily_cases), vjust = -0.5) +
  theme(text = element_text(family = "PingFang SC"))
