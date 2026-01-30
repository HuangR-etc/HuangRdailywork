# 加载时间序列分析必要的包
library(forecast)  # 用于时间序列分析和预测
library(tseries)   # 时间序列检验
library(zoo)       # 时间序列处理
library(Metrics)   # 评估指标

packages <- c("forecast","tseries","zoo","Metrics")
versions <- sapply(packages, function(x) as.character(packageVersion(x)))
versions

# 修改后的calculate_all_metrics函数 - 移除空间指标，增加时间序列指标
calculate_all_metrics <- function(time_data, time_points, 
                                  response_var = "reported_data", 
                                  predictor_var = "data",
                                  delta = 0.01) {
  
  # 提取观测值和预测值
  y_true <- time_data[[response_var]]
  y_pred <- time_data[[predictor_var]]
  n <- length(y_true)
  
  # 检查预测值是否为常数
  y_pred_sd <- sd(y_pred, na.rm = TRUE)
  is_constant_pred <- y_pred_sd == 0 || is.na(y_pred_sd)
  
  # 1. 基础误差指标
  residuals <- y_true - y_pred
  
  # OLS回归
  ols_r2 <- NA
  adjusted_r2 <- NA
  if (!is_constant_pred) {
    tryCatch({
      raw_model <- lm(as.formula(paste(response_var, "~", predictor_var)), data = time_data)
      ols_r2 <- summary(raw_model)$r.squared
      p <- 1
      adjusted_r2 <- 1 - (1 - ols_r2) * (n - 1) / (n - p - 1)
    }, error = function(e) {
      ols_r2 <<- NA
      adjusted_r2 <<- NA
    })
  }
  
  # 计算基于残差概率密度的加权RMSE
  weighted_rmse <- NA
  density_weights <- rep(1, n)  # 默认权重为1
  
  # 步骤1：计算残差的概率密度
  if (length(unique(residuals)) > 1) {  # 确保残差不全是相同的值
    tryCatch({
      # 计算残差的高斯核密度估计
      kde_density <- density(residuals, kernel = "gaussian", na.rm = TRUE)
      
      # 计算每个残差点处的概率密度
      get_density <- approxfun(kde_density$x, kde_density$y)
      prob_density <- get_density(residuals)
      
      # 处理可能的NA值（例如，残差在核密度估计范围外）
      prob_density[is.na(prob_density)] <- min(prob_density, na.rm = TRUE)
      
      # 计算原始权重（逆概率密度 + 平滑参数）
      raw_weights <- 1 / (prob_density + delta)
      
      # Min-Max归一化权重到[0,1]范围
      density_weights <- (raw_weights - min(raw_weights)) / 
        (max(raw_weights) - min(raw_weights))
      
      # 计算加权RMSE
      weighted_mse <- sum(density_weights * residuals^2, na.rm = TRUE) / 
        sum(density_weights, na.rm = TRUE)
      weighted_rmse <- sqrt(weighted_mse)
      
    }, error = function(e) {
      cat("加权RMSE计算错误:", e$message, "\n")
      weighted_rmse <<- sqrt(mean(residuals^2, na.rm = TRUE))  # 回退到标准RMSE
      density_weights <<- rep(1, n)
    })
  } else {
    # 如果残差几乎相同，退回到标准RMSE
    weighted_rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
    density_weights <- rep(1, n)
  }
  
  # 误差指标
  mse <- mean(residuals^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(residuals))
  median_ae <- median(abs(residuals))
  
  # 相关性指标
  pearson_corr <- NA
  spearman_corr <- NA
  if (!is_constant_pred) {
    tryCatch({
      pearson_corr <- cor(y_true, y_pred, method = "pearson")
      spearman_corr <- cor(y_true, y_pred, method = "spearman")
    }, error = function(e) {
      pearson_corr <<- NA
      spearman_corr <<- NA
    })
  }
  
  # 百分比误差指标
  safe_division <- function(a, b) {
    ifelse(b == 0, ifelse(a == 0, 0, Inf), a / b)
  }
  mape <- mean(abs(safe_division(residuals, y_true))) * 100
  
  # SMAPE指标
  numerator <- abs(residuals)
  denominator <- (abs(y_true) + abs(y_pred)) / 2
  denominator[denominator == 0] <- 1e-10
  smape_val <- mean(numerator / denominator) * 100
  cn_smape <- 1 - (smape_val / 200)
  
  # 新增指标: 百分比标准误差 (%SEE) 
  # 步骤1: 计算标准估计误差 (SEE)
  # p代表模型参数个数
  p_for_see <- 2  # 对应于简单线性回归 (截距 + 斜率)
  see <- sqrt(sum(residuals^2) / (n - p_for_see))
  # 步骤2: 计算观测值平均值
  y_mean <- mean(y_true)
  # 步骤3: 计算百分比标准误差
  percent_see <- (see / y_mean) * 100
  
  # 新增指标: 观测值:预测值比值及其汇总统计
  # 核心公式: R_i = y_i / yhat_i
  R_i <- safe_division(y_true, y_pred)  # 使用safe_division避免除以0
  
  # 汇总统计量
  R_whole <- sum(y_true)/sum(y_pred)
  R_arithmetic_mean <- mean(R_i, na.rm = TRUE)  # 算术平均值
  R_geometric_mean <- exp(mean(log(R_i), na.rm = TRUE))  # 几何平均值
  R_median <- median(R_i, na.rm = TRUE)  # 中位数比值
  
  
  # 2. 时间序列自相关指标
  
  # 真实值的自相关系数
  acf_real_lag1 <- NA
  tryCatch({
    acf_real <- acf(y_true, plot = FALSE, lag.max = 1, na.action = na.pass)
    acf_real_lag1 <- acf_real$acf[2]
  }, error = function(e) {
    acf_real_lag1 <<- NA
  })
  
  # 预测值的自相关系数
  acf_pred_lag1 <- NA
  if (!is_constant_pred) {
    tryCatch({
      acf_pred <- acf(y_pred, plot = FALSE, lag.max = 1, na.action = na.pass)
      acf_pred_lag1 <- acf_pred$acf[2]
    }, error = function(e) {
      acf_pred_lag1 <<- NA
    })
  }
  
  # 残差的自相关系数
  acf_residual_lag1 <- NA
  tryCatch({
    acf_resid <- acf(residuals, plot = FALSE, lag.max = 1, na.action = na.pass)
    acf_residual_lag1 <- acf_resid$acf[2]
  }, error = function(e) {
    acf_residual_lag1 <<- NA
  })
  
  # Ljung-Box检验（残差自相关检验）
  lb_test_pvalue <- NA
  tryCatch({
    lb_test <- Box.test(residuals, lag = 1, type = "Ljung-Box")
    lb_test_pvalue <- lb_test$p.value
  }, error = function(e) {
    lb_test_pvalue <<- NA
  })
  
  # 残差的ADF检验（平稳性检验）
  adf_test_pvalue <- NA
  tryCatch({
    adf_test <- adf.test(residuals, alternative = "stationary")
    adf_test_pvalue <- adf_test$p.value
  }, error = function(e) {
    adf_test_pvalue <<- NA
  })
  
  # 时间序列模型的R²（ARIMA模型）
  arima_r2 <- NA
  tryCatch({
    # 拟合ARIMA模型
    arima_model <- auto.arima(y_true, seasonal = FALSE)
    arima_residuals <- residuals(arima_model)
    arima_r2 <- max(0, 1 - var(arima_residuals)/var(y_true))
  }, error = function(e) {
    arima_r2 <<- NA
  })
  
  # 3. 时间序列预测指标
  # MASE (Mean Absolute Scaled Error) - 时间序列预测常用指标
  mase_value <- NA
  tryCatch({
    # 朴素预测作为基准
    naive_forecast <- c(NA, y_true[1:(n-1)])
    mae_naive <- mean(abs(y_true[2:n] - naive_forecast[2:n]), na.rm = TRUE)
    mase_value <- mae / mae_naive
  }, error = function(e) {
    mase_value <<- NA
  })
  
  # 返回所有指标
  return(list(
    # 基础指标
    ols_r2 = ols_r2,
    adjusted_r2 = adjusted_r2,
    mse = mse,
    rmse = rmse,
    # 新增指标：weighted_rmse
    weighted_rmse = weighted_rmse,  # 基于残差概率密度的加权RMSE
    mae = mae,
    median_ae = median_ae,
    
    # 相关性指标
    pearson_corr = pearson_corr,
    spearman_corr = spearman_corr,
    
    # 百分比误差指标
    mape = mape,
    smape = smape_val,
    cn_smape = cn_smape,
    
    # 新增: 百分比标准误差
    see = see,
    y_mean = y_mean,
    percent_see = percent_see,
    
    # 新增: 观测值与预测值比值
    R_whole = R_whole,
    R_ari = R_arithmetic_mean, #算数均值
    R_geo = R_geometric_mean, #几何均值
    R_median = R_median, #中位数
    
    # 时间序列自相关指标
    acf_real_lag1 = acf_real_lag1,
    acf_pred_lag1 = acf_pred_lag1,
    acf_residual_lag1 = acf_residual_lag1,
    lb_test_pvalue = lb_test_pvalue,
    adf_test_pvalue = adf_test_pvalue,
    arima_r2 = arima_r2,
    mase = mase_value,
    
    # 添加一个标志，指示预测值是否为常数
    is_constant_pred = is_constant_pred
  ))
}

# 模仿scikit-learn的TimeSeriesSplit逻辑的时间序列分割分析函数
perform_time_split_analysis <- function(time_data, n_splits = 5, 
                                        response_var = "reported_data", predictor_var = "data",
                                        scenario_name = "Unknown", sim_id = 1) {
  
  results <- list()
  n_total <- nrow(time_data)
  
  # 确保分割次数合理
  n_splits <- min(n_splits, n_total - 1)
  if (n_splits < 2) {
    stop("时间序列太短，无法进行有效分割")
  }
  
  # 生成时间点
  time_points <- 1:n_total
  
  # 模仿scikit-learn的TimeSeriesSplit逻辑
  # 每次分割，训练集逐步扩展，测试集紧随其后
  for (split_idx in 1:n_splits) {
    # 使用tryCatch进行错误处理
    result <- tryCatch({
      # 计算训练集和测试集的大小
      # 类似scikit-learn，训练集逐步扩展
      train_size <- round(n_total * (split_idx / (n_splits + 1)))
      test_size <- round(n_total * (1 / (n_splits + 1)))
      
      # 确保最小训练集大小
      min_train_size <- 20
      train_size <- max(train_size, min_train_size)
      
      # 确保有足够的测试数据
      if (train_size + test_size > n_total) {
        test_size <- n_total - train_size
      }
      
      if (test_size < 5) {
        next  # 测试集太小，跳过此次分割
      }
      
      # 确定训练集和测试集索引
      train_indices <- 1:train_size
      test_indices <- (train_size + 1):(train_size + test_size)
      
      # 创建训练集和测试集
      train_data <- time_data[train_indices, ]
      test_data <- time_data[test_indices, ]
      
      # 对应的时序点
      train_time <- time_points[train_indices]
      test_time <- time_points[test_indices]
      
      # 计算训练集指标
      metrics_train <- calculate_all_metrics(train_data, train_time,
                                             response_var = response_var, predictor_var = predictor_var)
      
      # 计算测试集指标
      metrics_test <- calculate_all_metrics(test_data, test_time,
                                            response_var = response_var, predictor_var = predictor_var)
      
      # 构建结果行
      result_row <- data.frame(
        scenario = scenario_name,
        simulation_id = sim_id,
        split_iteration = split_idx,
        train_size = length(train_indices),
        test_size = length(test_indices),
        train_start = min(train_indices),
        train_end = max(train_indices),
        test_start = min(test_indices),
        test_end = max(test_indices),
        
        # 训练集指标
        ols_r2_train = metrics_train$ols_r2,
        adjusted_r2_train = metrics_train$adjusted_r2,
        mse_train = metrics_train$mse,
        rmse_train = metrics_train$rmse,
        weighted_rmse_train = metrics_train$weighted_rmse,
        mae_train = metrics_train$mae,
        median_ae_train = metrics_train$median_ae,
        pearson_corr_train = metrics_train$pearson_corr,
        spearman_corr_train = metrics_train$spearman_corr,
        mape_train = metrics_train$mape,
        smape_train = metrics_train$smape,
        cn_smape_train = metrics_train$cn_smape,
        acf_real_lag1_train = metrics_train$acf_real_lag1,
        acf_pred_lag1_train = metrics_train$acf_pred_lag1,
        acf_residual_lag1_train = metrics_train$acf_residual_lag1,
        lb_test_pvalue_train = metrics_train$lb_test_pvalue,
        adf_test_pvalue_train = metrics_train$adf_test_pvalue,
        arima_r2_train = metrics_train$arima_r2,
        mase_train = metrics_train$mase,
        # 新增: 百分比标准误差
        see_train = metrics_train$see,
        percent_see_train = metrics_train$percent_see,
        
        # 新增: 观测值与预测值比值
        R_whole_train = metrics_train$R_whole,
        R_ari_train = metrics_train$R_ari, #算数均值
        R_geo_train = metrics_train$R_geo, #几何均值
        R_median_train = metrics_train$R_median, #中位数
        
        # 测试集指标
        ols_r2_test = metrics_test$ols_r2,
        adjusted_r2_test = metrics_test$adjusted_r2,
        mse_test = metrics_test$mse,
        rmse_test = metrics_test$rmse,
        weighted_rmse_test = metrics_test$weighted_rmse,
        mae_test = metrics_test$mae,
        median_ae_test = metrics_test$median_ae,
        pearson_corr_test = metrics_test$pearson_corr,
        spearman_corr_test = metrics_test$spearman_corr,
        mape_test = metrics_test$mape,
        smape_test = metrics_test$smape,
        cn_smape_test = metrics_test$cn_smape,
        acf_real_lag1_test = metrics_test$acf_real_lag1,
        acf_pred_lag1_test = metrics_test$acf_pred_lag1,
        acf_residual_lag1_test = metrics_test$acf_residual_lag1,
        lb_test_pvalue_test = metrics_test$lb_test_pvalue,
        adf_test_pvalue_test = metrics_test$adf_test_pvalue,
        arima_r2_test = metrics_test$arima_r2,
        mase_test = metrics_test$mase,
        # 新增: 百分比标准误差
        see_test = metrics_test$see,
        percent_see_test = metrics_test$percent_see,
        
        # 新增: 观测值与预测值比值
        R_whole_test = metrics_test$R_whole,
        R_ari_test = metrics_test$R_ari, #算数均值
        R_geo_test = metrics_test$R_geo, #几何均值
        R_median_test = metrics_test$R_median
      )
      
      # 返回成功的结果
      list(success = TRUE, result = result_row)
      
    }, error = function(e) {
      list(success = FALSE, result = NULL)
    })
    
    # 如果成功，保存结果
    if (result$success) {
      results[[split_idx]] <- result$result
    }
  }
  
  if (length(results) < n_splits) {
    cat("警告: 场景", scenario_name, "模拟", sim_id, "只完成了", length(results), "次成功分割\n")
  }
  
  return(do.call(rbind, results))
}

# 修改后的主模拟函数（时间序列版本）
run_single_simulation <- function(sim_id) {
  cat("正在进行第", sim_id, "次时间序列模拟...\n")
  
  # 设置随机种子
  set.seed(123 + sim_id)
  
  # 1. 定义时间序列长度
  n <- 200  # 时间序列长度
  
  # 2. 生成ARMA(2,1)模型数据：具有自相关的时间序列
  real_ts <- arima.sim(
    model = list(order = c(2, 0, 1), ar = c(0.8, -0.3), ma = 0.4), 
    n = n
  ) + 0.05 * (1:n)  # 加上轻微趋势
  
  # 转换为时间序列对象
  real_ts <- ts(real_ts, start = 1, frequency = 12)
  
  # 3. 创建基础数据框
  time_data <- data.frame(
    time = 1:n,
    reported_data = as.numeric(real_ts)
  )
  
  # 4. 场景0：好预测值（微小随机误差）
  time_data$data_good <- time_data$reported_data + rnorm(n, mean = 0, sd = 0.1)
  
  # 5. 场景1：季节性预测失败 - 模型在特定季节（如冬季）表现差
  # 识别冬季月份（假设数据从1月开始，频率为12）
  winter_months <- (time_data$time %% 12) %in% c(12, 1, 2, 0)  # 12月、1月、2月
  
  # 在冬季月份添加系统性的预测偏差
  bias_magnitude <- 2 * sd(time_data$reported_data)
  time_data$data_bad_seasonal <- time_data$reported_data
  time_data$data_bad_seasonal[winter_months] <- time_data$data_bad_seasonal[winter_months] + bias_magnitude
  
  # 6. 场景2：残差具有自相关 - 模型未完全捕捉时间依赖结构
  # 生成具有自相关的残差序列
  
  residual_scale <- 3  # 控制残差的总体变异程度
  
  # 场景2：残差时间自相关很强
  residual_strong_ar <- arima.sim(
    model = list(order = c(1, 0, 0), ar = 0.9),  # AR(1)系数=0.9，强自相关
    n = n
  )
  # 标准化残差序列，使其方差为1，然后缩放
  residual_strong_ar <- as.numeric(scale(residual_strong_ar)) * residual_scale
  time_data$data_bad_residual_strong <- time_data$reported_data + residual_strong_ar
  
  # 场景4：残差时间自相关极弱（无自相关）
  # 使用ar=0生成纯白噪声，确保完全无自相关
  residual_weak_ar <- arima.sim(
    model = list(order = c(1, 0, 0), ar = 0),  # AR(1)系数=0，无自相关
    n = n
  )
  # 标准化残差序列，使其方差为1，然后缩放
  residual_weak_ar <- as.numeric(scale(residual_weak_ar)) * residual_scale
  time_data$data_bad_residual_weak <- time_data$reported_data + residual_weak_ar
  
  # 8. 场景5：朴素预测模型（直接使用上一个点的真实值做预测值）
  time_data$data_bad_naive <- c(NA, time_data$reported_data[1:(n-1)])
  time_data$data_bad_naive[1] <- time_data$reported_data[1]  # 第一个点用自身值
  
  # 9. 场景6：完全失败模型
  mean_data <- mean(time_data$reported_data)
  time_data$data_bad <- rep(mean_data, n) + rnorm(n, mean = 0, sd = 0.1)
  
  # 10. 计算各场景指标（完整数据集）
  scenarios <- c("Good", "UC1", "UC2_st","UC2_wk", "UC3", "Bad")
  predictor_vars <- c("data_good", "data_bad_seasonal", "data_bad_residual_strong",
                      "data_bad_residual_weak",
                     "data_bad_naive", "data_bad")
  
  # 移除包含NA的行（特别是朴素预测的第一个点）
  time_data_clean <- time_data[complete.cases(time_data), ]
  
  metrics_list <- list()
  for (i in 1:length(scenarios)) {
    metrics <- calculate_all_metrics(time_data_clean, time_data_clean$time,
                                     predictor_var = predictor_vars[i])
    metrics_list[[scenarios[i]]] <- metrics
  }
  
  # 11. 对每个场景进行TimeSeriesSplit分析
  split_results_list <- list()
  
  for (i in 1:length(scenarios)) {
    cat("  场景", scenarios[i], "TimeSeriesSplit分析...\n")
    
    split_results <- perform_time_split_analysis(
      time_data = time_data_clean,
      n_splits = 5,  # 使用5折时间序列交叉验证
      response_var = "reported_data",
      predictor_var = predictor_vars[i],
      scenario_name = scenarios[i],
      sim_id = sim_id
    )
    
    split_results_list[[scenarios[i]]] <- split_results
  }
  
  # 合并所有分割结果
  all_split_results_for_sim <- do.call(rbind, split_results_list)
  
  return(list(
    full_data_metrics = metrics_list,
    split_analysis_results = all_split_results_for_sim,
    time_data = time_data_clean
  ))
}

# 执行批量模拟
n_simulations <- 50
n_iterations <- 50
cat("开始时间序列自相关模拟，总共", n_simulations, "次模拟...\n")
start_time <- Sys.time()

# 使用lapply进行串行计算
all_sim_results <- lapply(1:n_simulations, function(i) {
  result <- run_single_simulation(i)
  cat("完成第", i, "次时间序列模拟\n")
  return(result)
})

end_time <- Sys.time()
cat("时间序列模拟完成，耗时:", round(end_time - start_time, 2), "秒\n")

# 整理结果
# 完整数据集结果（进行敏感性分析）
full_metrics_list <- lapply(all_sim_results, function(x) x$full_data_metrics)

# 分割数据集结果（进行泛化能力分析）
split_results_list <- lapply(all_sim_results, function(x) x$split_analysis_results)
all_split_results_df <- do.call(rbind, split_results_list)
#----------8.结果分析----------
# 8. 计算平均指标
calculate_average_metrics <- function(all_results, n_sim) {
  scenarios <- c("Good", "UC1",  "UC2_st","UC2_wk", "UC3", "Bad")
  
  # 获取所有指标名称
  metric_names <- names(all_results[[1]][[1]])
  
  # 创建存储平均结果的列表
  avg_results <- list()
  
  for (scenario in scenarios) {
    # 提取该场景在所有模拟中的指标
    scenario_metrics <- lapply(all_results, function(x) x[[scenario]])
    
    # 计算每个指标的平均值和标准误差
    avg_metrics <- list()
    for (metric in metric_names) {
      values <- sapply(scenario_metrics, function(x) {
        val <- x[[metric]]
        if (is.null(val) || is.na(val)) NA else val
      })
      
      # 移除NA值
      valid_values <- values[!is.na(values)]
      n_valid <- length(valid_values)
      
      if (n_valid > 0) {
        avg_metrics[[paste0(metric, "_mean")]] <- mean(valid_values)
        avg_metrics[[paste0(metric, "_se")]] <- sd(valid_values) / sqrt(n_valid)
        avg_metrics[[paste0(metric, "_n")]] <- n_valid
      } else {
        avg_metrics[[paste0(metric, "_mean")]] <- NA
        avg_metrics[[paste0(metric, "_se")]] <- NA
        avg_metrics[[paste0(metric, "_n")]] <- 0
      }
    }
    avg_results[[scenario]] <- avg_metrics
  }
  
  return(avg_results)
}

# 计算完整数据集的平均指标
avg_metrics <- calculate_average_metrics(full_metrics_list, n_simulations)

# # 9. 创建完整数据集的比较表格
create_final_comparison_table <- function(avg_metrics) {
  scenarios <- c("Good", "UC1",  "UC2_st","UC2_wk", "UC3", "Bad")
  
  # 获取所有指标名称（去掉后缀）
  all_names <- names(avg_metrics[[1]])
  metric_base_names <- unique(gsub("_(mean|se|n)$", "", all_names))
  
  # 创建结果数据框
  result_df <- data.frame(Scenario = scenarios)
  
  # 为每个指标添加均值和标准误差列
  for (metric in metric_base_names) {
    mean_col <- paste0(metric, "_mean")
    se_col <- paste0(metric, "_se")
    n_col <- paste0(metric, "_n")
    
    # 提取各场景的均值和标准误差
    means <- sapply(scenarios, function(s) avg_metrics[[s]][[mean_col]])
    ses <- sapply(scenarios, function(s) avg_metrics[[s]][[se_col]])
    ns <- sapply(scenarios, function(s) avg_metrics[[s]][[n_col]])
    
    result_df[[paste0(metric, "_mean")]] <- means
    result_df[[paste0(metric, "_se")]] <- ses
    result_df[[paste0(metric, "_n")]] <- ns
  }
  
  return(result_df)
}

final_comparison <- create_final_comparison_table(avg_metrics)

#-------9. 改进的可视化方案（基于多次模拟结果）----------

# 9.1 准备数据 - 从多次模拟结果中提取均值

# 创建简化的比较数据框（只包含均值）
comparison_means <- data.frame(Scenario = final_comparison$Scenario)

# 提取所有指标的均值列
metric_names <- gsub("_mean", "", grep("_mean$", names(final_comparison), value = TRUE))

for (metric in metric_names) {
  mean_col <- paste0(metric, "_mean")
  if (mean_col %in% names(final_comparison)) {
    comparison_means[[metric]] <- final_comparison[[mean_col]]
  }
}

# 提取标准误差用于误差线
comparison_se <- data.frame(Scenario = final_comparison$Scenario)
for (metric in metric_names) {
  se_col <- paste0(metric, "_se")
  if (se_col %in% names(final_comparison)) {
    comparison_se[[metric]] <- final_comparison[[se_col]]
  }
}

# 定义指标分类（重新组织，包含系统发育指标）
# 越大越好指标
higher_better_metrics <- c(
  "ols_r2", "adjusted_r2", "pearson_corr", "spearman_corr", "cn_smape",
  "pic_r2", "pic_pearson", "pic_spearman", "whitened_r2", "resid_r2", "lik_r2"
)

# 越小越好指标
lower_better_metrics <- c(
  "mse", "rmse", "weighted_rmse", "mae", "median_ae", "mape", "smape",
  "pic_rmse", "pic_mae", "whitened_rmse", "whitened_mae",
  "see", "percent_see"
)

# 越接近1越好指标
closer_to_one_metrics <- c("R_whole", "R_ari", "R_geo", "R_median")

# 筛选实际存在的指标
available_higher <- higher_better_metrics[higher_better_metrics %in% names(comparison_means)]
available_lower <- lower_better_metrics[lower_better_metrics %in% names(comparison_means)]
available_closer_to_one <- closer_to_one_metrics[closer_to_one_metrics %in% names(comparison_means)]

cat("基于", n_simulations, "次模拟的平均结果\n")
cat("可用的'越大越好'指标:", length(available_higher), "\n")
cat("可用的'越小越好'指标:", length(available_lower), "\n")
cat("可用的'越接近1越好'指标:", length(available_closer_to_one), "\n")

# 设置颜色方案
good_color <- "#c9cb05"  # 好模型
bad_colors <- c("#bab1d8","#9e7cba","#8076b5","#624c7c","#ac5aa1","#27447c")  # 紫色系 - 坏模型
scenario_colors <- c(good_color, bad_colors)

# 创建输出目录
output_dir <- "time_model_comparison_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("创建输出目录:", output_dir, "\n")
}

# 9.2 创建绘图函数
cat("\n=== 绘制带误差线的柱状图 ===\n")

# 创建带误差线的柱状图函数
plot_metric_with_errorbars <- function(metric, means_df, se_df, metric_type = "higher_better") {
  scenarios <- means_df$Scenario
  values <- means_df[[metric]]
  errors <- se_df[[metric]]
  
  # 设置图形参数
  ylab <- ""
  main_title <- paste(metric)
  
  # 根据指标类型设置不同的y轴范围
  if (metric_type == "closer_to_one") {
    # 对于接近1的指标，确保y轴包含1
    y_max <- max(values + errors, 1, na.rm = TRUE) * 1.15
    y_min <- min(values - errors, 0, na.rm = TRUE) * 0.85
    y_lim <- c(y_min, y_max)
  } else if (metric_type == "higher_better") {
    y_lim <- c(0, max(values + errors, na.rm = TRUE) * 1.15)
  } else {  # lower_better
    y_lim <- c(0, max(values + errors, na.rm = TRUE) * 1.15)
  }
  
  # 创建条形图
  bar_pos <- barplot(values, names.arg = scenarios, 
                     col = scenario_colors,
                     main = main_title,
                     ylab = ylab,
                     ylim = y_lim,
                     border = "white",
                     las = 2)
  
  # 添加误差线
  arrows(bar_pos, values - errors, bar_pos, values + errors, 
         angle = 90, code = 3, length = 0.05, lwd = 1.5)
  
  # 添加数值标签（均值）
  text(x = bar_pos, y = values + errors, 
       label = sprintf("%.3f", values), 
       pos = 3, cex = 0.6, col = "darkblue", font = 2)
}

# 9.3 创建导出图形的函数
create_figure_pdf <- function(metrics, metric_type, filename, n_cols = 3) {
  # 计算布局
  n_metrics <- length(metrics)
  n_rows <- ceiling(n_metrics / n_cols)
  
  # 设置图形尺寸（适合论文）
  # 每个子图宽度2.5英寸，高度2英寸
  fig_width <- n_cols * 2.5
  fig_height <- n_rows * 2.0
  
  # 打开PDF设备
  pdf(file.path(output_dir, filename), 
      width = fig_width, 
      height = fig_height)
  
  # 设置图形参数
  par(mfrow = c(n_rows, n_cols), 
      mar = c(3.5, 3, 2, 1),  # 减小边距
      oma = c(2,0,0,0), # 整个图形的外边界：下=2行，左=0，上=0，右=0
      mgp = c(1.5, 0.5, 0),   # 调整坐标轴标签位置
      cex.main = 0.9,         # 减小标题字体
      cex.axis = 0.7)         # 减小坐标轴字体
  
  # 绘制所有指标
  for (metric in metrics) {
    if (all(is.na(comparison_means[[metric]]))) next
    
    tryCatch({
      plot_metric_with_errorbars(metric, comparison_means, comparison_se, metric_type)
    }, error = function(e) {
      plot(1, type = "n", main = paste(metric, "\n绘图错误"), xlab = "", ylab = "")
    })
  }
  
  # 关闭设备
  dev.off()
  
  cat("已保存图形到:", filename, "尺寸:", round(fig_width, 1), "x", round(fig_height, 1), "英寸\n")
}

# 9.4 创建单个图形（适合在R中查看）
create_figure_screen <- function(metrics, metric_type, n_cols = 3) {
  # 计算布局
  n_metrics <- length(metrics)
  n_rows <- ceiling(n_metrics / n_cols)
  
  # 设置图形参数
  par(mfrow = c(n_rows, n_cols), 
      mar = c(4, 4, 3, 1), 
      mgp = c(2, 0.8, 0))
  
  # 绘制所有指标
  for (metric in metrics) {
    if (all(is.na(comparison_means[[metric]]))) next
    
    tryCatch({
      plot_metric_with_errorbars(metric, comparison_means, comparison_se, metric_type)
    }, error = function(e) {
      plot(1, type = "n", main = paste(metric, "\n绘图错误"), xlab = "", ylab = "")
    })
  }
}

# 9.5 智能调整每行子图数量函数
smart_adjust_ncols <- function(metrics, max_per_page = 12, min_ncols = 2) {
  n_metrics <- length(metrics)
  
  if (n_metrics <= 4) {
    # 指标较少时，每行显示2个
    return(2)
  } else if (n_metrics <= 9) {
    # 指标适中时，每行显示3个
    return(3)
  } else {
    # 指标很多时，每行显示4个
    return(4)
  }
}

# 9.6 导出图形
cat("\n=== 导出图形到文件 ===\n")

# 9.6.1 导出"越大越好"指标
if (length(available_higher) > 0) {
  cat("--- 导出'越大越好'指标 ---\n")
  
  # 智能调整每行子图数量
  n_cols_higher <- smart_adjust_ncols(available_higher)
  cat("每行显示子图数量:", n_cols_higher, "\n")
  
  # 在屏幕上显示
  create_figure_screen(available_higher, "越大越好", n_cols = n_cols_higher)
  
  # 导出PDF（适合论文）
  create_figure_pdf(available_higher, "越大越好", "higher_better_metrics.pdf", n_cols = n_cols_higher)
} else {
  cat("没有可用的'越大越好'指标\n")
}

# 9.6.2 导出"越小越好"指标
if (length(available_lower) > 0) {
  cat("--- 导出'越小越好'指标 ---\n")
  
  # 智能调整每行子图数量
  n_cols_lower <- smart_adjust_ncols(available_lower)
  cat("每行显示子图数量:", n_cols_lower, "\n")
  
  # 在屏幕上显示
  create_figure_screen(available_lower, "越小越好", n_cols = n_cols_lower)
  
  # 导出PDF
  create_figure_pdf(available_lower, "越小越好", "lower_better_metrics.pdf", n_cols = n_cols_lower)
} else {
  cat("没有可用的'越小越好'指标\n")
}

# 9.6.3 导出"越接近1越好"指标
if (length(available_closer_to_one) > 0) {
  cat("--- 导出'越接近1越好'指标 ---\n")
  
  # 对于接近1的指标，通常数量较少，每行显示2个比较合适
  n_cols_closer <- 2
  
  # 在屏幕上显示
  n_rows <- ceiling(length(available_closer_to_one) / n_cols_closer)
  par(mfrow = c(n_rows, n_cols_closer), 
      mar = c(4, 4, 3, 1), 
      mgp = c(2, 0.8, 0))
  
  for (metric in available_closer_to_one) {
    if (all(is.na(comparison_means[[metric]]))) next
    
    tryCatch({
      plot_metric_with_errorbars(metric, comparison_means, comparison_se, "closer_to_one")
    }, error = function(e) {
      plot(1, type = "n", main = paste(metric, "\n绘图错误"), xlab = "", ylab = "")
    })
  }
  
  # 导出PDF
  n_rows <- ceiling(length(available_closer_to_one) / n_cols_closer)
  fig_width <- n_cols_closer * 2.5
  fig_height <- n_rows * 2.0
  
  pdf(file.path(output_dir, "closer_to_one_metrics.pdf"), 
      width = fig_width, 
      height = fig_height)
  
  par(mfrow = c(n_rows, n_cols_closer), 
      mar = c(3.5, 3, 2, 1),
      mgp = c(1.5, 0.5, 0),
      cex.main = 0.9,
      cex.axis = 0.7)
  
  for (metric in available_closer_to_one) {
    if (all(is.na(comparison_means[[metric]]))) next
    
    tryCatch({
      plot_metric_with_errorbars(metric, comparison_means, comparison_se, "closer_to_one")
    }, error = function(e) {
      plot(1, type = "n", main = paste(metric, "\n绘图错误"), xlab = "", ylab = "")
    })
  }
  
  dev.off()
  cat("已保存图形到: closer_to_one_metrics.pdf, 尺寸:", round(fig_width, 1), "x", round(fig_height, 1), "英寸\n")
} else {
  cat("没有可用的'越接近1越好'指标\n")
}

# 重置图形参数
par(mfrow = c(1, 1))

# 9.7 创建组合图（将所有指标放在一个文件中）
cat("\n=== 创建所有指标的组合图 ===\n")

# 创建包含所有指标的PDF文件
create_all_metrics_pdf <- function() {
  all_metrics <- c(available_higher, available_lower, available_closer_to_one)
  metric_types <- c(
    rep("higher_better", length(available_higher)),
    rep("lower_better", length(available_lower)),
    rep("closer_to_one", length(available_closer_to_one))
  )
  
  n_metrics <- length(all_metrics)
  
  # 智能调整列数
  if (n_metrics <= 9) {
    n_cols <- 3
  } else if (n_metrics <= 16) {
    n_cols <- 4
  } else {
    n_cols <- 5
  }
  
  n_rows <- ceiling(n_metrics / n_cols)
  
  # 计算图形尺寸
  fig_width <- n_cols * 2.2  # 稍微减小宽度，使图形更紧凑
  fig_height <- n_rows * 1.8  # 稍微减小高度
  
  pdf(file.path(output_dir, "all_metrics_comparison.pdf"), 
      width = fig_width, 
      height = fig_height)
  
  par(mfrow = c(n_rows, n_cols), 
      mar = c(3, 2.5, 2, 1),
      mgp = c(1.2, 0.4, 0),
      cex.main = 0.8,
      cex.axis = 0.7)
  
  for (i in 1:n_metrics) {
    metric <- all_metrics[i]
    metric_type <- metric_types[i]
    
    if (all(is.na(comparison_means[[metric]]))) next
    
    tryCatch({
      plot_metric_with_errorbars(metric, comparison_means, comparison_se, metric_type)
    }, error = function(e) {
      plot(1, type = "n", main = paste(metric, "\n绘图错误"), xlab = "", ylab = "")
    })
  }
  
  dev.off()
  cat("已保存所有指标组合图到: all_metrics_comparison.pdf, 尺寸:", round(fig_width, 1), "x", round(fig_height, 1), "英寸\n")
}

# 如果至少有一个指标，创建组合图
if (length(available_higher) + length(available_lower) + length(available_closer_to_one) > 0) {
  create_all_metrics_pdf()
}

# 9.8 添加汇总说明
cat("\n=== 指标分类说明 ===\n")
cat("1. 越大越好指标:", paste(available_higher, collapse = ", "), "\n")
cat("2. 越小越好指标:", paste(available_lower, collapse = ", "), "\n")
cat("3. 越接近1越好指标:", paste(available_closer_to_one, collapse = ", "), "\n")
cat("\n图形已保存到目录:", output_dir, "\n")
cat("文件列表:\n")
cat("  - higher_better_metrics.pdf: 越大越好指标\n")
cat("  - lower_better_metrics.pdf: 越小越好指标\n")
cat("  - closer_to_one_metrics.pdf: 越接近1越好指标\n")
cat("  - all_metrics_comparison.pdf: 所有指标组合图\n")
#-------泛化能力部分的可视化--------
library(ggplot2)
library(dplyr)
calculate_delta_metrics <- function(split_results_df) {
  # 定义各类指标
  r2_metrics <- c("ols", "adjusted", "pic", "whitened", "resid", "lik")
  error_metrics <- c("mse", "rmse","weighted_rmse", "mae", "median_ae", 
                     "pic_rmse", "pic_mae", "whitened_rmse", "whitened_mae",
                     "see","percent_see","R_whole","R_ari","R_geo","R_median")
  corr_metrics <- c("pearson", "spearman", "pic_pearson", "pic_spearman")
  percentage_metrics <- c("mape", "smape", "cn_smape")
  
  # 创建结果数据框副本
  result_with_delta <- split_results_df
  
  # 计算delta = (test - train) / [(train + test)/2]
  delta_columns <- c()
  
  # 辅助函数：安全计算delta
  safe_calculate_delta <- function(train_vals, test_vals) {
    # 检查是否有NA值
    if(any(is.na(train_vals)) || any(is.na(test_vals))) {
      return(rep(NA, length(train_vals)))
    }
    
    # 计算分母
    denominator <- (train_vals + test_vals) / 2
    
    # 安全处理小分母：只对非NA值进行操作
    non_na_idx <- !is.na(denominator)
    denominator_safe <- denominator
    
    if(any(non_na_idx)) {
      small_denom_idx <- non_na_idx & (abs(denominator[non_na_idx]) < 1e-10)
      if(any(small_denom_idx)) {
        # 只对有问题的值进行修正
        denominator_safe[non_na_idx][small_denom_idx] <- 
          sign(denominator[non_na_idx][small_denom_idx]) * 1e-10
      }
    }
    
    # 计算delta
    delta <- (test_vals - train_vals) / denominator_safe
    return(delta)
  }
  
  # 1. R²类指标
  for(metric in r2_metrics) {
    train_col <- paste0(metric, "_r2_train")
    test_col <- paste0(metric, "_r2_test")
    delta_col <- paste0("delta_", metric, "_r2")
    
    if(train_col %in% colnames(result_with_delta) && test_col %in% colnames(result_with_delta)) {
      result_with_delta[[delta_col]] <- safe_calculate_delta(
        result_with_delta[[train_col]], 
        result_with_delta[[test_col]]
      )
      delta_columns <- c(delta_columns, delta_col)
    }
  }
  
  # 2. 误差类指标
  for(metric in error_metrics) {
    train_col <- paste0(metric, "_train")
    test_col <- paste0(metric, "_test")
    delta_col <- paste0("delta_", metric)
    
    if(train_col %in% colnames(result_with_delta) && test_col %in% colnames(result_with_delta)) {
      result_with_delta[[delta_col]] <- safe_calculate_delta(
        result_with_delta[[train_col]], 
        result_with_delta[[test_col]]
      )
      delta_columns <- c(delta_columns, delta_col)
    }
  }
  
  # 3. 相关性指标
  for(metric in corr_metrics) {
    # 检查两种可能的列名格式
    train_col1 <- paste0(metric, "_corr_train")
    test_col1 <- paste0(metric, "_corr_test")
    train_col2 <- paste0(metric, "_train")
    test_col2 <- paste0(metric, "_test")
    
    if(train_col1 %in% colnames(result_with_delta) && test_col1 %in% colnames(result_with_delta)) {
      result_with_delta[[paste0("delta_", metric, "_corr")]] <- safe_calculate_delta(
        result_with_delta[[train_col1]], 
        result_with_delta[[test_col1]]
      )
      delta_columns <- c(delta_columns, paste0("delta_", metric, "_corr"))
    } else if(train_col2 %in% colnames(result_with_delta) && test_col2 %in% colnames(result_with_delta)) {
      result_with_delta[[paste0("delta_", metric)]] <- safe_calculate_delta(
        result_with_delta[[train_col2]], 
        result_with_delta[[test_col2]]
      )
      delta_columns <- c(delta_columns, paste0("delta_", metric))
    }
  }
  
  # 4. 百分比误差指标
  for(metric in percentage_metrics) {
    train_col <- paste0(metric, "_train")
    test_col <- paste0(metric, "_test")
    delta_col <- paste0("delta_", metric)
    
    if(train_col %in% colnames(result_with_delta) && test_col %in% colnames(result_with_delta)) {
      result_with_delta[[delta_col]] <- safe_calculate_delta(
        result_with_delta[[train_col]], 
        result_with_delta[[test_col]]
      )
      delta_columns <- c(delta_columns, delta_col)
    }
  }
  
  return(list(
    data = result_with_delta,
    delta_columns = delta_columns
  ))
}

# 先检查数据中是否有NA值
cat("检查分割结果数据中的NA值情况:\n")
cat("总行数:", nrow(all_split_results_df), "\n")

# 检查关键列是否有NA
check_na_columns <- function(df) {
  na_counts <- sapply(df, function(x) sum(is.na(x)))
  na_columns <- na_counts[na_counts > 0]
  if(length(na_columns) > 0) {
    cat("包含NA值的列:\n")
    print(na_columns)
  } else {
    cat("没有发现NA值\n")
  }
}

check_na_columns(all_split_results_df)

# 应用修复后的函数计算delta指标
cat("\n应用修复后的函数计算delta指标...\n")
delta_results <- calculate_delta_metrics(all_split_results_df)
delta_data <- delta_results$data
delta_columns <- delta_results$delta_columns

cat("成功计算了", length(delta_columns), "个delta指标\n")

# 检查delta指标中的NA情况
if(length(delta_columns) > 0) {
  delta_na_counts <- sapply(delta_columns, function(col) {
    if(col %in% names(delta_data)) {
      sum(is.na(delta_data[[col]]))
    } else {
      NA
    }
  })
  cat("各delta指标的NA值数量:\n")
  print(delta_na_counts)
}
# 计算每个场景和每个指标的delta方差
calculate_delta_variance <- function(delta_data, delta_columns) {
  scenarios <- unique(delta_data$scenario)
  
  variance_results <- data.frame()
  
  for(scenario in scenarios) {
    scenario_data <- delta_data[delta_data$scenario == scenario, ]
    
    for(metric in delta_columns) {
      values <- scenario_data[[metric]]
      values <- values[!is.na(values)]  # 移除NA值
      
      if(length(values) > 1) {
        variance_val <- var(values)
        mean_val <- mean(values)
        abs_mean_val <- mean(abs(values))
        
        # 确定指标类型
        if(grepl("_r2$", metric)) {
          metric_type <- "R-squared"
        } else if(grepl("mse|rmse|weighted_rmse|mae|median_ae|see|percent_see|R_whole|R_ari|R_geo|R_median", metric)) {
          metric_type <- "Error"
        } else if(grepl("pearson|spearman", metric)) {
          metric_type <- "Correlation"
        } else if(grepl("mape|smape|cn_smape", metric)) {
          metric_type <- "Percentage Error"
        } else if(grepl("lambda", metric)) {
          metric_type <- "Phylogenetic Signal"
        } else {
          metric_type <- "Other"
        }
        
        # 简化指标名称用于显示
        clean_metric_name <- gsub("^delta_", "", metric)
        
        variance_results <- rbind(variance_results, data.frame(
          Scenario = scenario,
          Metric = metric,
          CleanMetric = clean_metric_name,
          Type = metric_type,
          Variance = variance_val,
          MeanDelta = mean_val,
          AbsMeanDelta = abs_mean_val,
          N = length(values)
        ))
      }
    }
  }
  
  return(variance_results)
}

# 计算方差
variance_results <- calculate_delta_variance(delta_data, delta_columns)

# 创建保存结果的子目录
output_dir <- "time_model_comparison_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("创建输出目录:", output_dir, "\n")
}

#-------修改后的方差排序图绘制逻辑（按场景选择刻度类型）----------

# 按场景分别绘制方差排序柱状图并保存为PDF
plot_and_save_variance_barchart <- function(variance_results, output_dir) {
  scenarios <- unique(variance_results$Scenario)
  saved_files <- list()
  
  for(scenario in scenarios) {
    # 筛选当前场景的数据
    scenario_data <- variance_results[variance_results$Scenario == scenario, ]
    
    # 按方差排序
    scenario_ranked <- scenario_data[order(-scenario_data$Variance), ]
    
    # 设置颜色（每个场景使用单一颜色）
    scenario_colors <- c(
      "Good" = "#c9cb05", 
      "UC1" = "#bab1d8", 
      "UC2_st" = "#9e7cba", 
      "UC2_wk" = "#624c7c",
      "UC3" = "#9e7cba", 
      "Bad" = "#ac5aa1"
    )
    
    current_color <- scenario_colors[scenario]
    
    # 确保方差值为正值（对数刻度需要正值）
    # 将零值替换为一个很小的正数
    min_positive_variance <- min(scenario_ranked$Variance[scenario_ranked$Variance > 0], na.rm = TRUE)
    if(is.infinite(min_positive_variance) || is.na(min_positive_variance)) {
      min_positive_variance <- 1e-10
    }
    scenario_ranked$Variance_adj <- ifelse(scenario_ranked$Variance <= 0, 
                                           min_positive_variance * 0.1, 
                                           scenario_ranked$Variance)
    
    # 判断是否使用对数刻度
    use_log_scale <- !(scenario %in% c("Good", "UC1"))
    
    if(use_log_scale) {
      cat("场景", scenario, "使用对数刻度\n")
      
      # 计算对数刻度下的合适范围
      log_values <- log10(scenario_ranked$Variance_adj)
      log_max <- max(log_values, na.rm = TRUE)
      log_min <- min(log_values, na.rm = TRUE)
      
      # 扩展范围用于标签显示
      log_range <- log_max - log_min
      log_upper_limit <- log_max + (log_range * 0.05)  # 增加5%的空间用于标签
      
      # 转换为实际值
      upper_limit <- 10^(log_upper_limit)
      
      p <- ggplot(scenario_ranked, aes(x = reorder(CleanMetric, Variance_adj), y = Variance_adj)) +
        geom_bar(stat = "identity", fill = current_color, alpha = 0.8, width = 0.7) +
        coord_flip() +
        labs(
          title = paste("Scenario", scenario, "- Delta Metric Variance Ranking"),
          x = "Metric",
          y = "Variance (log10 scale)"
        ) +
        # 使用对数刻度
        scale_y_log10(
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = function(x) {
            ifelse(x >= 1000, 
                   formatC(x, format = "e", digits = 1),
                   formatC(x, format = "f", digits = 3))
          },
          expand = expansion(mult = c(0, 0.1)),
          limits = c(NA, upper_limit)  # 设置上限
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.y = element_text(size = 11, color = "black", face = "bold"),
          axis.text.x = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
          plot.margin = margin(1, 2, 1, 1.2, "cm"),  # 增加右边距
          panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
          panel.grid.minor = element_blank()
        )
      
    } else {
      cat("场景", scenario, "使用线性刻度\n")
      
      # 对于线性刻度，直接使用原始方差值（调整后的值）
      # 计算线性刻度下的合适范围
      y_max <- max(scenario_ranked$Variance_adj, na.rm = TRUE)
      y_min <- min(scenario_ranked$Variance_adj, na.rm = TRUE)
      
      # 扩展范围用于标签显示
      y_range <- y_max - y_min
      upper_limit <- y_max + (y_range * 0.1)  # 增加10%的空间用于标签
      
      p <- ggplot(scenario_ranked, aes(x = reorder(CleanMetric, Variance_adj), y = Variance_adj)) +
        geom_bar(stat = "identity", fill = current_color, alpha = 0.8, width = 0.7) +
        coord_flip() +
        labs(
          title = paste("Scenario", scenario, "- Delta Metric Variance Ranking"),
          x = "Metric",
          y = "Variance"
        ) +
        # 使用线性刻度
        scale_y_continuous(
          labels = function(x) {
            # 对于线性刻度，使用更友好的格式
            if(any(x >= 1000, na.rm = TRUE)) {
              # 如果有大数值，部分使用科学计数法
              return(ifelse(x >= 1000, 
                            formatC(x, format = "e", digits = 1),
                            sprintf("%.3f", x)))
            } else {
              # 否则使用常规小数格式
              return(sprintf("%.3f", x))
            }
          },
          expand = expansion(mult = c(0, 0.1)),
          limits = c(NA, upper_limit)  # 设置上限
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.y = element_text(size = 11, color = "black", face = "bold"),
          axis.text.x = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
          plot.margin = margin(1, 2, 1, 1.2, "cm"),  # 增加右边距
          panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
          panel.grid.minor = element_blank()
        )
    }
    
    # 添加数值标签（根据刻度类型调整格式）
    if(use_log_scale) {
      p <- p + geom_text(
        aes(label = ifelse(Variance_adj >= 1000, 
                           sprintf("%.2e", Variance_adj), 
                           sprintf("%.3f", Variance_adj))), 
        hjust = -0.1, 
        size = 3.5, 
        color = "darkblue"
      )
    } else {
      # 对于线性刻度，使用更简单的格式
      p <- p + geom_text(
        aes(label = ifelse(Variance_adj >= 0.001 & Variance_adj < 1000, 
                           sprintf("%.4f", Variance_adj),
                           sprintf("%.2e", Variance_adj))), 
        hjust = -0.1, 
        size = 3.5, 
        color = "darkblue"
      )
    }
    
    pdf_file <- file.path(output_dir, paste0("scenario_", scenario, "_delta_variance_", 
                                             ifelse(use_log_scale, "log", "linear"), ".pdf"))
    
    # 保存为PDF
    ggsave(pdf_file, p, width = 10, height = 7, dpi = 300)
    
    cat("已保存:", pdf_file, "\n")
    saved_files[[scenario]] <- pdf_file
    
    # 也显示在R中查看
    print(p)
  }
  
  return(saved_files)
}

# 保存方差数据为CSV文件
save_variance_data <- function(variance_results, output_dir) {
  # 保存详细数据
  csv_file <- file.path(output_dir, "delta_variance_results.csv")
  write.csv(variance_results, csv_file, row.names = FALSE)
  cat("已保存方差数据:", csv_file, "\n")
  
  # 创建汇总表格（按场景和指标类型）
  summary_data <- variance_results %>%
    group_by(Scenario, Type) %>%
    summarise(
      Avg_Variance = mean(Variance, na.rm = TRUE),
      Max_Variance = max(Variance, na.rm = TRUE),
      Min_Variance = min(Variance, na.rm = TRUE),
      Avg_AbsDelta = mean(AbsMeanDelta, na.rm = TRUE),
      n_Metrics = n(),
      .groups = 'drop'
    )
  
  summary_file <- file.path(output_dir, "delta_variance_summary.csv")
  write.csv(summary_data, summary_file, row.names = FALSE)
  cat("已保存汇总数据:", summary_file, "\n")
  
  return(list(details = csv_file, summary = summary_file))
}

# 执行绘图和保存
cat("=== 保存Delta指标可视化图表和数据 ===\n")

# 方差排序柱状图（按场景分开保存为PDF）
cat("1. 按场景绘制方差排序柱状图（Good和UC1使用线性刻度，其他使用对数刻度）...\n")
pdf_files <- plot_and_save_variance_barchart(variance_results, output_dir)

# 保存数据为CSV
cat("\n2. 保存方差数据为CSV文件...\n")
csv_files <- save_variance_data(variance_results, output_dir)

# 输出总结信息
cat("\n=== 处理完成 ===")
cat("\n输出目录:", output_dir)
cat("\n保存的PDF文件:")
for(scenario in names(pdf_files)) {
  cat("\n  - 场景", scenario, ":", basename(pdf_files[[scenario]]))
}
cat("\n保存的CSV文件:")
cat("\n  - 详细数据:", basename(csv_files$details))
cat("\n  - 汇总数据:", basename(csv_files$summary))
cat("\n")
#-------新增：收集方差排序结果，与另一分析流程格式保持一致----------
# 定义数据集ID（根据实际情况修改）
dataset_id <- "simulation_data"  # 可以根据需要修改，如"simulation_1", "simulation_2"等

# 为每个场景生成方差排名结果并保存
cat("\n=== 收集方差排序结果，与下游分析流程保持一致 ===\n")

# 初始化一个列表来存储所有场景的排名结果
all_ranking_results <- list()

# 对每个场景进行处理
scenarios <- unique(variance_results$Scenario)

for(scenario in scenarios) {
  cat("处理场景:", scenario, "\n")
  
  # 筛选当前场景的数据
  scenario_data <- variance_results[variance_results$Scenario == scenario, ]
  
  # 1. 计算当前场景的方差排名（方差越小越好，所以用升序排名）
  # 注意：使用负号将方差转换为降序，这样方差小的排名靠前
  scenario_data$Variance_Rank <- rank(scenario_data$Variance, ties.method = "min")
  
  # 2. 按方差从小到大排序（方差越小越好）
  scenario_data_sorted <- scenario_data[order(scenario_data$Variance), ]
  
  # 3. 创建排名结果数据框（与另一分析流程格式完全相同）
  ranking_results <- data.frame(
    Metric = scenario_data_sorted$CleanMetric,  # 使用简化后的指标名称
    Type = scenario_data_sorted$Type,
    Variance = scenario_data_sorted$Variance,
    Rank = 1:nrow(scenario_data_sorted),  # 重新编号从1开始
    Dataset = paste0(dataset_id, "_", scenario)  # 格式：数据集ID_场景名
  )
  
  # 4. 打印当前场景的排名结果摘要
  cat("场景", scenario, "方差排名前5的指标:\n")
  print(head(ranking_results, 5))
  cat("\n")
  
  # 5. 保存当前场景的排名结果到文件
  ranking_filename <- paste0("time_var_rank_", dataset_id, "_", scenario, ".csv")
  write.csv(ranking_results, ranking_filename, row.names = FALSE)
  cat("场景", scenario, "排名结果已保存到:", ranking_filename, "\n\n")
  
  # 6. 将当前场景的结果添加到总列表中
  all_ranking_results[[scenario]] <- ranking_results
}