# 加载必要的包
library(forecast)
library(tseries)
library(mgcv)
library(zoo)
# 设置随机种子以确保结果可重现
set.seed(123)

# 生成AR(1)模型数据：强正自相关 (系数phi=0.8)
real_ts <- arima.sim(
  model = list(order = c(2, 0, 1), ar = c(0.8, -0.3), ma = 0.4), 
  n = n
) + 0.05 * (1:n)  # 加上轻微趋势

# 转换为时间序列对象
real_ts <- ts(real_ts, start = 1, frequency = 12)
# 创建基础数据框
time_data <- data.frame(
  time = time(real_ts),
  reported_data = as.numeric(real_ts)
)
plot(time_data$time,time_data$reported_data)

#------场景1：季节性预测失败 - 模型在特定季节（如冬季）表现差--------
seasonal_failure <- real_ts

# 识别冬季月份（假设数据从1月开始，频率为12）
winter_months <- cycle(real_ts) %in% c(12, 1, 2)

# 在冬季月份添加系统性的预测偏差
bias_magnitude <- 2 * sd(real_ts)
seasonal_failure[winter_months] <- real_ts[winter_months] + bias_magnitude

time_data$data_bad_seasonal <- as.numeric(seasonal_failure)

#------场景2：残差具有自相关 - 模型未完全捕捉时间依赖结构---------
# 生成具有自相关的残差序列
residual_signal <- arima.sim(
  model = list(order = c(1, 0, 0), ar = 0.7), 
  n = n
)

pred_bad_residual <- real_ts + residual_signal * 0.5
time_data$data_bad_residual <- as.numeric(pred_bad_residual)

#------场景3：结构性断点预测失败 - 模型未适应数据生成过程的变化-------
pred_bad_breakpoint <- real_ts

# 创建一个结构性断点（例如，第60个时间点后数据生成过程改变）
break_point <- 70
# 在断点前拟合线性趋势
train_data <- data.frame(
  y = as.numeric(real_ts[1:break_point]),
  t = 1:break_point
)
linear_model <- lm(y ~ t, data = train_data)

# 断点后继续使用该线性模型预测
for(i in (break_point+1):n) {
  pred_bad_breakpoint[i] <- predict(linear_model, 
                                    newdata = data.frame(t = i))
}
time_data$data_bad_breakpoint <- pred_bad_breakpoint
plot(time_data$reported_data,time_data$data_bad_breakpoint)

#---------场景4：朴素预测模型（直接使用上一个点的真实值做预测值）------------
# 朴素预测（使用前一期真实值）
pred_naive <- rep(NA, n)  # 初始化预测向量

# 第一个点没有前一期，可以用第一个真实值本身或均值
pred_naive[1] <- real_ts[1]  # 第一种选择：用第一个真实值

# 对于第2到第n个点，使用前一期真实值作为预测
for(i in 2:n) {
  pred_naive[i] <- real_ts[i-1]
}

# 添加到数据框
time_data$data_bad_naive <- as.numeric(pred_naive)
#---------场景5：完全失败模拟------------
mean_data <- mean(time_data$reported_data)
pred_bad <- rep(mean_data,100)
pred_bad <- pred_bad + rnorm(100, mean = 0, sd = 0.1)
time_data$data_bad <- pred_bad
