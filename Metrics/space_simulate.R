setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/UC")
# 加载空间分析必要的包
library(gstat)    # 用于空间模拟和变异函数
library(sp)       # 处理空间数据
library(spdep)    # 空间自相关分析
library(ape)      # 莫兰指数计算
library(nlme)     # 空间线性模型
library(fields)   # 用于距离计算和可视化
library(cluster)  # 聚类

# calculate_all_metrics函数
calculate_all_metrics <- function(spatial_data, coords_matrix, 
                                  response_var = "reported_data", 
                                  predictor_var = "data",
                                  delta = 0.01) {     # 新增参数：平滑参数
  
  # 提取观测值和预测值
  y_true <- spatial_data[[response_var]]
  y_pred <- spatial_data[[predictor_var]]
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
      raw_model <- lm(as.formula(paste(response_var, "~", predictor_var)), data = spatial_data)
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
  
  # 普通误差指标
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
  
  # 2. 空间自相关指标
  
  # 计算距离矩阵的逆矩阵（用于莫兰指数）
  spatial_dist <- as.matrix(dist(coords_matrix))
  spatial_dist_inv <- 1/spatial_dist
  diag(spatial_dist_inv) <- 0
  
  # 真实值的莫兰指数
  moran_real <- NA
  tryCatch({
    moran_result_real <- Moran.I(y_true, spatial_dist_inv)
    moran_real <- moran_result_real$observed
  }, error = function(e) {
    moran_real <<- NA
  })
  
  # 预测值的莫兰指数
  moran_pred <- NA
  if (!is_constant_pred) {
    tryCatch({
      moran_result_pred <- Moran.I(y_pred, spatial_dist_inv)
      moran_pred <- moran_result_pred$observed
    }, error = function(e) {
      moran_pred <<- NA
    })
  }
  
  # 残差的莫兰指数
  moran_residual <- NA
  tryCatch({
    moran_result_resid <- Moran.I(residuals, spatial_dist_inv)
    moran_residual <- moran_result_resid$observed
  }, error = function(e) {
    moran_residual <<- NA
  })
  
  # 空间误差模型的R²
  spatial_gls_r2 <- NA
  tryCatch({
    # 创建空间数据框
    spatial_df <- spatial_data
    spatial_df$x_coord <- coords_matrix[,1]
    spatial_df$y_coord <- coords_matrix[,2]
    
    # 拟合空间误差模型
    gls_model <- gls(as.formula(paste(response_var, "~", predictor_var)), 
                     data = spatial_df,
                     correlation = corSpher(form = ~ x_coord + y_coord, nugget = TRUE),
                     method = "REML")
    
    # 计算空间模型的R²
    null_model <- gls(as.formula(paste(response_var, "~ 1")), 
                      data = spatial_df,
                      method = "REML")
    
    spatial_gls_r2 <- max(0, 1 - var(residuals(gls_model))/var(residuals(null_model)))
  }, error = function(e) {
    spatial_gls_r2 <<- NA
  })
  
  # 变异函数分析指标
  variogram_range <- NA
  variogram_sill <- NA
  variogram_nugget <- NA
  tryCatch({
    # 计算残差的实验变异函数
    spatial_df$residuals <- residuals
    coordinates(spatial_df) <- ~ x_coord + y_coord
    emp_variogram <- variogram(residuals ~ 1, data = spatial_df)
    
    # 拟合理论变异函数模型
    fit_variogram <- fit.variogram(emp_variogram, model = vgm(1, "Sph", 0.3, 0.1))
    
    if (nrow(fit_variogram) > 1) {
      variogram_range <- fit_variogram$range[2]
      variogram_sill <- fit_variogram$psill[2]
      variogram_nugget <- fit_variogram$psill[1]
    }
  }, error = function(e) {
    # 如果变异函数拟合失败，保持NA值
  })
  
  # 返回所有指标
  return(list(
    # 基础指标
    ols_r2 = ols_r2,
    adjusted_r2 = adjusted_r2,
    mse = mse,
    rmse = rmse,
    mae = mae,
    median_ae = median_ae,
    # 新增指标：weighted_rmse
    weighted_rmse = weighted_rmse,  # 基于残差概率密度的加权RMSE
    
    # 相关性指标
    pearson_corr = pearson_corr,
    spearman_corr = spearman_corr,
    
    # 百分比误差指标
    mape = mape,
    smape = smape_val,
    cn_smape = cn_smape,
    
    # 空间自相关指标
    moran_real = moran_real,
    moran_pred = moran_pred,
    moran_residual = moran_residual,
    spatial_gls_r2 = spatial_gls_r2,
    variogram_range = variogram_range,
    variogram_sill = variogram_sill,
    variogram_nugget = variogram_nugget,
    
    # 添加一个标志，指示预测值是否为常数
    is_constant_pred = is_constant_pred
  ))
}

# 随机分割分析函数（空间版本）
perform_random_split_analysis <- function(spatial_data, coords_matrix, n_iterations = 50, train_ratio = 0.8, 
                                          response_var = "reported_data", predictor_var = "data",
                                          scenario_name = "Unknown",
                                          sim_id = 1) {
  
  results <- list()
  successful_iterations <- 0
  iteration_count <- 0
  
  while (successful_iterations < n_iterations && iteration_count < n_iterations * 2) {
    iteration_count <- iteration_count + 1
    
    # 设置随机种子
    current_seed <- 123 + sim_id * 1000 + iteration_count
    set.seed(current_seed)
    
    # 使用tryCatch进行错误处理
    result <- tryCatch({
      # 随机划分数据
      n_total <- nrow(spatial_data)
      n_train <- round(n_total * train_ratio)
      n_test <- n_total - n_train
      
      train_indices <- sample(1:n_total, n_train, replace = FALSE)
      test_indices <- setdiff(1:n_total, train_indices)
      
      # 创建训练集和测试集
      train_data <- spatial_data[train_indices, ]
      test_data <- spatial_data[test_indices, ]
      
      # 对应的坐标
      train_coords <- coords_matrix[train_indices, ]
      test_coords <- coords_matrix[test_indices, ]
      
      # 计算训练集指标
      metrics_train <- calculate_all_metrics(train_data, train_coords,
                                             response_var = response_var, predictor_var = predictor_var)
      
      # 计算测试集指标
      metrics_test <- calculate_all_metrics(test_data, test_coords,
                                            response_var = response_var, predictor_var = predictor_var)
      
      # 构建结果行
      result_row <- data.frame(
        scenario = scenario_name,
        simulation_id = sim_id,
        split_iteration = successful_iterations + 1,
        random_seed = current_seed,
        train_size = n_train,
        test_size = n_test,
        
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
        moran_real_train = metrics_train$moran_real,
        moran_pred_train = metrics_train$moran_pred,
        moran_residual_train = metrics_train$moran_residual,
        spatial_gls_r2_train = metrics_train$spatial_gls_r2,
        variogram_range_train = metrics_train$variogram_range,
        variogram_sill_train = metrics_train$variogram_sill,
        variogram_nugget_train = metrics_train$variogram_nugget,
        
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
        moran_real_test = metrics_test$moran_real,
        moran_pred_test = metrics_test$moran_pred,
        moran_residual_test = metrics_test$moran_residual,
        spatial_gls_r2_test = metrics_test$spatial_gls_r2,
        variogram_range_test = metrics_test$variogram_range,
        variogram_sill_test = metrics_test$variogram_sill,
        variogram_nugget_test = metrics_test$variogram_nugget
      )
      
      # 返回成功的结果
      list(success = TRUE, result = result_row)
      
    }, error = function(e) {
      list(success = FALSE, result = NULL)
    })
    
    # 如果成功，保存结果
    if (result$success) {
      successful_iterations <- successful_iterations + 1
      results[[successful_iterations]] <- result$result
    }
  }
  
  if (successful_iterations < n_iterations) {
    cat("警告: 场景", scenario_name, "模拟", sim_id, "只完成了", successful_iterations, "次成功分割\n")
  }
  
  return(do.call(rbind, results))
}

# 主模拟函数（空间版本）
run_single_simulation <- function(sim_id) {
  cat("正在进行第", sim_id, "次空间模拟...\n")
  
  # 设置随机种子
  set.seed(123 + sim_id)
  
  # 1. 定义空间场景
  n_locations <- 100
  
  # 在单位平面内随机生成空间坐标
  spatial_coords <- data.frame(
    x = runif(n_locations),
    y = runif(n_locations)
  )
  coords_matrix <- as.matrix(spatial_coords)
  
  # 为每个点分配ID
  spatial_coords$location_id <- paste0("loc", 1:n_locations)
  
  # 2. 计算空间距离矩阵
  spatial_dist <- as.matrix(dist(coords_matrix))
  
  # 3. 定义主导空间过程
  dominant_variogram <- vgm(
    psill = 1.0,
    model = "Sph",
    range = 0.3,
    nugget = 0.1
  )
  
  # 4. 生成主导空间过程
  sp_points <- SpatialPoints(spatial_coords[, c("x", "y")])
  spatial_field <- predict(gstat(
    formula = z ~ 1,
    locations = ~ x + y,
    dummy = TRUE,
    beta = 0,
    model = dominant_variogram
  ), newdata = sp_points, nsim = 1)
  
  # 创建基础数据框
  base_data <- data.frame(
    location_id = spatial_coords$location_id,
    x_coord = spatial_coords$x,
    y_coord = spatial_coords$y,
    reported_data = spatial_field@data$sim1  # 具有空间自相关的"真实值"
  )
  
  # 计算主导过程的方差，用于控制残差场的强度
  dominant_variance <- var(base_data$reported_data)
  cat("主导过程方差:", round(dominant_variance, 3), "\n")
  
  # 5. 生成好预测值（微小随机误差）
  base_data$data_good <- base_data$reported_data + rnorm(n_locations, mean = 0, sd = 0.1)
  
  # 6. 场景2：空间聚集误差
  spatial_clusters <- pam(spatial_dist, k = 4)
  biased_region <- sample(1:4, 1)
  bias_magnitude <- 2 * sd(base_data$reported_data)
  
  base_data$data_bad_region <- base_data$reported_data
  base_data$data_bad_region[spatial_clusters$clustering == biased_region] <- 
    base_data$data_bad_region[spatial_clusters$clustering == biased_region] + bias_magnitude
  
  # 7. 场景3-5：残差具有不同强度的空间自相关
  # 目标：残差场方差占主导过程方差的80%（确保不掩盖主导过程）
  residual_variance_ratio <- 0.8
  target_residual_variance <- dominant_variance * residual_variance_ratio
  
  # 7.1 强空间自相关残差场 (90%结构方差，10%块金)
  strong_resid_model <- vgm(
    psill = target_residual_variance * 0.9,  # 90%结构方差
    model = "Sph",
    range = 0.3,  # 稍小于主导过程的变程
    nugget = target_residual_variance * 0.1  # 10%块金
  )
  strong_resid_field <- predict(gstat(
    formula = z ~ 1,
    locations = ~ x + y,
    dummy = TRUE,
    beta = 0,
    model = strong_resid_model
  ), newdata = sp_points, nsim = 1)
  base_data$data_strong_spatialerr <- base_data$reported_data + strong_resid_field@data$sim1
  
  # 7.3 极弱空间自相关残差场 (纯块金模型，无空间自相关)
  weak_resid_model <- vgm(
    psill = 0,  # 无结构方差
    model = "Nug",  # 块金模型
    range = 0,  # 无变程
    nugget = target_residual_variance  # 全部为块金方差
  )
  weak_resid_field <- predict(gstat(
    formula = z ~ 1,
    locations = ~ x + y,
    dummy = TRUE,
    beta = 0,
    model = weak_resid_model
  ), newdata = sp_points, nsim = 1)
  base_data$data_weak_spatialerr <- base_data$reported_data + weak_resid_field@data$sim1
  
  # 8. 场景6：边界效应误差
  boundary_dist <- pmin(base_data$x_coord, 1-base_data$x_coord,
                        base_data$y_coord, 1-base_data$y_coord)
  boundary_effect <- boundary_dist/max(boundary_dist)
  base_data$data_boundary_bias <- base_data$reported_data + 
    rnorm(nrow(base_data), 0, sd(base_data$reported_data)*0.5) * (1-boundary_effect)
  
  # 9. 场景7：完全失败模型
  mean_data <- mean(base_data$reported_data)
  base_data$data_bad <- rep(mean_data, n_locations) + rnorm(n_locations, mean = 0, sd = 0.1)
  
  # 10. 计算各场景指标（完整数据集）
  scenarios <- c("Good", "UC1", "UC2_st", "UC2_wk", "UC3", "Bad")
  predictor_vars <- c("data_good", "data_bad_region", "data_strong_spatialerr", 
                      "data_weak_spatialerr", 
                      "data_boundary_bias", "data_bad")
  
  metrics_list <- list()
  for (i in 1:length(scenarios)) {
    metrics <- calculate_all_metrics(base_data, coords_matrix,
                                     predictor_var = predictor_vars[i])
    metrics_list[[scenarios[i]]] <- metrics
  }
  
  # 11. 对每个场景进行分割分析
  split_results_list <- list()
  
  for (i in 1:length(scenarios)) {
    cat("  场景", scenarios[i], "分割分析...\n")
    
    split_results <- perform_random_split_analysis(
      spatial_data = base_data,
      coords_matrix = coords_matrix,
      n_iterations = 20,  # 减少迭代次数以加快测试
      train_ratio = 0.8,
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
    base_data = base_data,
    coords_matrix = coords_matrix
  ))
}

# 执行批量模拟（减少模拟次数以加快测试）
n_simulations <- 30
cat("开始空间自相关模拟，总共", n_simulations, "次模拟...\n")
start_time <- Sys.time()

# 使用lapply进行串行计算
all_sim_results <- lapply(1:n_simulations, function(i) {
  result <- run_single_simulation(i)
  cat("完成第", i, "次空间模拟\n")
  return(result)
})

end_time <- Sys.time()
cat("空间模拟完成，耗时:", round(end_time - start_time, 2), "秒\n")

# 整理结果
#完整数据集结果（进行敏感性分析）
full_metrics_list <- lapply(all_sim_results, function(x) x$full_data_metrics)
#分割数据集结果（进行泛化能力分析）
split_results_list <- lapply(all_sim_results, function(x) x$split_analysis_results)
all_split_results_df <- do.call(rbind, split_results_list)

#----------结果分析----------
# 8. 计算平均指标
calculate_average_metrics <- function(all_results, n_sim) {
  scenarios <- c("Good", "UC1", "UC2_st", "UC2_wk","UC3", "Bad")
  
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
  scenarios <- c("Good", "UC1",  "UC2_st", "UC2_wk", "UC3", "Bad")
  
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

# 定义指标分类（与之前相同）
higher_better_metrics <- c(
  "ols_r2", "adjusted_r2", "pearson_corr", "spearman_corr", "cn_smape"
)

lower_better_metrics <- c(
  "mse", "rmse", "weighted_rmse", "mae", "median_ae", "mape", "smape"
)

# 筛选实际存在的指标
available_higher <- higher_better_metrics[higher_better_metrics %in% names(comparison_means)]
available_lower <- lower_better_metrics[lower_better_metrics %in% names(comparison_means)]

cat("基于", n_simulations, "次模拟的平均结果\n")
cat("可用的'越大越好'指标:", length(available_higher), "\n")
cat("可用的'越小越好'指标:", length(available_lower), "\n")

# 设置颜色方案
good_color <- "#c9cb05"  # 好模型
bad_colors <- c("#bab1d8","#9e7cba","#8076b5","#624c7c","#ac5aa1","#27447c")  # 紫色系 - 坏模型
scenario_colors <- c(good_color, bad_colors)

# 9.2 改进的柱状图（带误差线）
cat("\n=== 绘制带误差线的柱状图 ===\n")

# 创建带误差线的柱状图函数
plot_metric_with_errorbars <- function(metric, means_df, se_df, metric_type = "higher_better") {
  scenarios <- means_df$Scenario
  values <- means_df[[metric]]
  errors <- se_df[[metric]]
  
  # 设置y轴标签和标题
  ylab <- ""
  if (metric_type == "higher_better") {
    main_title <- paste(metric)
  } else {
    main_title <- paste(metric)
  }
  
  # 创建条形图
  bar_pos <- barplot(values, names.arg = scenarios, 
                     col = scenario_colors,
                     main = main_title,
                     ylab = ylab,
                     ylim = c(0, max(values + errors, na.rm = TRUE) * 1.15),
                     border = "white",
                     las = 2)
  
  # 添加误差线
  arrows(bar_pos, values - errors, bar_pos, values + errors, 
         angle = 90, code = 3, length = 0.05, lwd = 1.5)
  
  # 添加数值标签（均值）
  text(x = bar_pos, y = values + errors, 
       label = sprintf("%.3f", values), 
       pos = 3, cex = 0.6, col = "darkblue", font = 2)
  
  # 添加标准误差标签
  text(x = bar_pos, y = values - errors - max(values, na.rm = TRUE) * 0.02, 
       label = sprintf("±%.3f", errors), 
       pos = 1, cex = 0.5, col = "darkred")
}

# 9.2.1 绘制"越大越好"指标
cat("--- 绘制'越大越好'指标（带误差线）---\n")

n_higher <- length(available_higher)
n_cols <- 3  # 减少每行列数以便更好地显示
n_rows <- ceiling(n_higher / n_cols)

# 设置图形参数
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), mgp = c(2, 0.8, 0))

for (metric in available_higher) {
  if (all(is.na(comparison_means[[metric]]))) next
  
  tryCatch({
    plot_metric_with_errorbars(metric, comparison_means, comparison_se, "higher_better")
  }, error = function(e) {
    plot(1, type = "n", main = paste(metric, "\n绘图错误"), xlab = "", ylab = "")
  })
}

# mtext(paste("'越大越好'型指标比较（基于", n_simulations, "次模拟）"), 
#       side = 3, outer = TRUE, line = -1.5, cex = 1.2)

# 9.2.2 绘制"越小越好"指标
cat("--- 绘制'越小越好'指标（带误差线）---\n")

n_lower <- length(available_lower)
n_cols <- 3
n_rows <- ceiling(n_lower / n_cols)

par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), mgp = c(2, 0.8, 0))

for (metric in available_lower) {
  if (all(is.na(comparison_means[[metric]]))) next
  
  tryCatch({
    plot_metric_with_errorbars(metric, comparison_means, comparison_se, "lower_better")
  }, error = function(e) {
    plot(1, type = "n", main = paste(metric, "\n绘图错误"), xlab = "", ylab = "")
  })
}

# mtext(paste("'越小越好'型指标比较（基于", n_simulations, "次模拟）"), 
#       side = 3, outer = TRUE, line = -1.5, cex = 1.2)

# 重置图形参数
par(mfrow = c(1, 1))


#-------泛化能力部分的可视化--------
library(ggplot2)
library(dplyr)

calculate_delta_metrics <- function(split_results_df) {
  # 定义各类指标
  r2_metrics <- c("ols", "adjusted")
  error_metrics <- c("mse", "rmse","weighted_rmse", "mae", "median_ae")
  corr_metrics <- c("pearson", "spearman")
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
        } else if(grepl("mse|rmse|weighted_rmse|mae|median_ae", metric)) {
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

#-------修改后的方差排序图绘制逻辑----------

# 按场景分别绘制方差排序柱状图
plot_variance_barchart_by_scenario <- function(variance_results) {
  scenarios <- unique(variance_results$Scenario)
  plots <- list()
  
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
      "UC3" = "#ac5aa1", 
      "Bad" = "#27447c"
    )
    
    current_color <- scenario_colors[scenario]
    
    # 创建图形
    p <- ggplot(scenario_ranked, aes(x = reorder(CleanMetric, Variance), y = Variance)) +
      geom_bar(stat = "identity", fill = current_color, alpha = 0.8) +
      coord_flip() +
      labs(
        title = paste("场景", scenario, "- Delta指标方差排序"),
        subtitle = "Delta = (测试集 - 训练集) / 均值",
        x = "指标",
        y = "方差"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "none"
      ) +
      # 添加数值标签
      geom_text(aes(label = sprintf("%.4f", Variance)), 
                hjust = -0.1, size = 3, color = "darkblue")
    
    plots[[scenario]] <- p
  }
  
  return(plots)
}

# 执行分开绘图
cat("=== 按场景分开绘制Delta指标可视化图表 ===\n")

# 方差排序柱状图（按场景分开）
cat("1. 按场景绘制方差排序柱状图...\n")
variance_plots_by_scenario <- plot_variance_barchart_by_scenario(variance_results)

# 显示每个场景的方差排序图
for(scenario in names(variance_plots_by_scenario)) {
  cat("显示场景", scenario, "的方差排序图...\n")
  print(variance_plots_by_scenario[[scenario]])
}
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
  ranking_filename <- paste0("space_var_rank_", dataset_id, "_", scenario, ".csv")
  write.csv(ranking_results, ranking_filename, row.names = FALSE)
  cat("场景", scenario, "排名结果已保存到:", ranking_filename, "\n\n")
  
  # 6. 将当前场景的结果添加到总列表中
  all_ranking_results[[scenario]] <- ranking_results
}