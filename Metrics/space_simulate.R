setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/UC")
# 加载空间分析必要的包
library(gstat)    # 用于空间模拟和变异函数
library(sp)       # 处理空间数据
library(ape)      # 莫兰指数计算
library(nlme)     # 空间线性模型
library(cluster)  # 聚类
packages <- c("gstat", "sp", "spdep", "ape", "nlme", "fields", "cluster")
versions <- sapply(packages, function(x) as.character(packageVersion(x)))
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
  
  # 新增指标: 百分比标准误差 (%SEE)
  # 步骤1: 计算标准估计误差 (SEE)
  # 注意: 这里p代表模型参数个数
  p_for_see <- 2  # 对应于简单线性回归 (截距 + 斜率)
  see <- sqrt(sum(residuals^2) / (n - p_for_see))
  # 步骤2: 计算观测值平均值
  y_mean <- mean(y_true)
  # 步骤3: 计算百分比标准误差
  percent_see <- (see / y_mean) * 100
  
  # 新增指标: 观测值:预测值比值及其汇总统计
  R_whole <- sum(y_true)/sum(y_pred)
  R_i <- safe_division(y_true, y_pred)  # 使用safe_division避免除以0
  
  # 汇总统计量
  R_arithmetic_mean <- mean(R_i, na.rm = TRUE)  # 算术平均值
  R_geometric_mean <- exp(mean(log(R_i), na.rm = TRUE))  # 几何平均值
  R_median <- median(R_i, na.rm = TRUE)  # 中位数比值
  
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
    
    # 新增: 百分比标准误差
    see = see,
    y_mean = y_mean,
    percent_see = percent_see,
    
    # 新增: 观测值与预测值比值
    R_whole = R_whole,
    R_ari = R_arithmetic_mean, #算数均值
    R_geo = R_geometric_mean, #几何均值
    R_median = R_median, #中位数
    
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
        moran_real_test = metrics_test$moran_real,
        moran_pred_test = metrics_test$moran_pred,
        moran_residual_test = metrics_test$moran_residual,
        spatial_gls_r2_test = metrics_test$spatial_gls_r2,
        variogram_range_test = metrics_test$variogram_range,
        variogram_sill_test = metrics_test$variogram_sill,
        variogram_nugget_test = metrics_test$variogram_nugget,
        # 新增: 百分比标准误差
        see_test = metrics_test$see,
        percent_see_test = metrics_test$percent_see,
        
        # 新增: 观测值与预测值比值
        R_whole_test = metrics_test$R_whole,
        R_ari_test = metrics_test$R_ari, #算数均值
        R_geo_test = metrics_test$R_geo, #几何均值
        R_median = metrics_test$R_median
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

# 新增：空间分块分割分析函数
perform_spatial_block_split_kfold <- function(spatial_data, coords_matrix, n_blocks = 4,
                                              response_var = "reported_data", predictor_var = "data",
                                              scenario_name = "Unknown", sim_id = 1) {
  
  # 计算空间距离矩阵
  spatial_dist <- as.matrix(dist(coords_matrix))
  
  # 使用k-medoids聚类（PAM）基于空间距离分成n_blocks块
  pam_result <- pam(spatial_dist, k = n_blocks)
  cluster_assignments <- pam_result$clustering
  
  # 记录块的大小
  block_sizes <- table(cluster_assignments)
  cat("  场景", scenario_name, "模拟", sim_id, ": 空间分块大小 = ", 
      paste(block_sizes, collapse = ", "), "\n")
  
  # 存储每次分割的结果
  all_results <- list()
  
  # 为每个点添加ID
  spatial_data$point_id <- paste0("point", 1:nrow(spatial_data))
  
  # 对每一块进行交叉验证
  for (test_block in 1:n_blocks) {
    # 划分训练集和测试集
    test_points <- which(cluster_assignments == test_block)
    train_points <- which(cluster_assignments != test_block)
    
    # 确保训练集和测试集不重叠
    if(length(intersect(train_points, test_points)) > 0) {
      stop("训练集和测试集有重叠！")
    }
    
    # 创建训练集和测试集
    train_data <- spatial_data[train_points, ]
    test_data <- spatial_data[test_points, ]
    
    # 确保训练集和测试集都不为空
    if(nrow(train_data) == 0 || nrow(test_data) == 0) {
      warning("训练集或测试集为空！跳过这次分割。")
      next
    }
    
    # 获取训练集和测试集的坐标
    train_coords <- coords_matrix[train_points, ]
    test_coords <- coords_matrix[test_points, ]
    
    # 计算训练集和测试集的指标
    metrics_train <- calculate_all_metrics(train_data, train_coords,
                                           response_var = response_var, predictor_var = predictor_var)
    
    metrics_test <- calculate_all_metrics(test_data, test_coords,
                                          response_var = response_var, predictor_var = predictor_var)
    
    # 构建结果行
    result_row <- data.frame(
      scenario = scenario_name,
      simulation_id = sim_id,
      split_iteration = test_block,  # 使用块编号作为分割迭代
      split_method = "spatial_block",  # 标记为空间分块分割
      train_size = nrow(train_data),
      test_size = nrow(test_data),
      train_proportion = round(nrow(train_data)/nrow(spatial_data), 3),
      test_proportion = round(nrow(test_data)/nrow(spatial_data), 3),
      block_id = test_block,
      n_blocks = n_blocks,
      
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
      # 新增: 百分比标准误差
      see_train = metrics_train$see,
      percent_see_train = metrics_train$percent_see,
      
      # 新增: 观测值与预测值比值
      R_whole_train = metrics_train$R_whole,
      R_ari_train = metrics_train$R_ari,
      R_geo_train = metrics_train$R_geo,
      R_median_train = metrics_train$R_median,
      
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
      variogram_nugget_test = metrics_test$variogram_nugget,
      # 新增: 百分比标准误差
      see_test = metrics_test$see,
      percent_see_test = metrics_test$percent_see,
      
      # 新增: 观测值与预测值比值
      R_whole_test = metrics_test$R_whole,
      R_ari_test = metrics_test$R_ari,
      R_geo_test = metrics_test$R_geo,
      R_median_test = metrics_test$R_median
    )
    
    all_results[[test_block]] <- result_row
  }
  
  # 合并所有结果
  if (length(all_results) > 0) {
    return(do.call(rbind, all_results))
  } else {
    return(NULL)
  }
}

# 更新主模拟函数，支持两种分割方法
run_single_simulation_with_both_methods <- function(sim_id, n_random_iterations = 50, 
                                                    n_blocks = 4, train_ratio = 0.8) {
  cat("正在进行第", sim_id, "次空间模拟（包含随机分割和空间分块分割）...\n")
  
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
  
  # 存储所有结果
  random_split_results_list <- list()
  block_split_results_list <- list()
  
  # 对每个场景进行分析
  for (i in 1:length(scenarios)) {
    scenario <- scenarios[i]
    predictor_var <- predictor_vars[i]
    
    cat("  场景", scenario, "分析...\n")
    
    # 1. 随机分割分析
    tryCatch({
      random_results <- perform_random_split_analysis(
        spatial_data = base_data,
        coords_matrix = coords_matrix,
        n_iterations = n_random_iterations,
        train_ratio = train_ratio,
        response_var = "reported_data",
        predictor_var = predictor_var,
        scenario_name = scenario,
        sim_id = sim_id
      )
      
      random_split_results_list[[scenario]] <- random_results
    }, error = function(e) {
      cat("    场景", scenario, "随机分割失败:", e$message, "\n")
    })
    
    # 2. 空间分块分割分析
    tryCatch({
      block_results <- perform_spatial_block_split_kfold(
        spatial_data = base_data,
        coords_matrix = coords_matrix,
        n_blocks = n_blocks,
        response_var = "reported_data",
        predictor_var = predictor_var,
        scenario_name = scenario,
        sim_id = sim_id
      )
      
      if (!is.null(block_results)) {
        block_split_results_list[[scenario]] <- block_results
      } else {
        cat("    场景", scenario, "空间分块分割返回NULL结果\n")
      }
    }, error = function(e) {
      cat("    场景", scenario, "空间分块分割失败:", e$message, "\n")
    })
  }
  
  # 计算完整数据集的指标
  metrics_list <- list()
  for (i in 1:length(scenarios)) {
    metrics <- calculate_all_metrics(base_data, coords_matrix,
                                     predictor_var = predictor_vars[i])
    metrics_list[[scenarios[i]]] <- metrics
  }
  
  # 分别合并结果
  random_split_df <- NULL
  block_split_df <- NULL
  
  if (length(random_split_results_list) > 0) {
    random_split_df <- do.call(rbind, random_split_results_list)
    rownames(random_split_df) <- NULL
  }
  
  if (length(block_split_results_list) > 0) {
    block_split_df <- do.call(rbind, block_split_results_list)
    rownames(block_split_df) <- NULL
  }
  
  return(list(
    full_data_metrics = metrics_list,
    random_split_results = random_split_df,
    block_split_results = block_split_df,
    base_data = base_data,
    coords_matrix = coords_matrix
  ))
}

# 执行批量模拟
n_simulations <- 50
n_random_iterations <- 50
n_blocks <- 4

cat("开始空间模拟，总共", n_simulations, "次模拟（包含随机分割和空间分块分割）...\n")
start_time <- Sys.time()

# 使用lapply进行串行计算
all_sim_results <- lapply(1:n_simulations, function(i) {
  result <- run_single_simulation_with_both_methods(
    i, 
    n_random_iterations = n_random_iterations,
    n_blocks = n_blocks
  )
  cat("完成第", i, "次模拟\n")
  return(result)
})

end_time <- Sys.time()
cat("空间模拟完成，耗时:", round(end_time - start_time, 2), "秒\n")

# 整理结果
# 提取完整数据集结果
full_metrics_list <- lapply(all_sim_results, function(x) x$full_data_metrics)

# 提取随机分割结果
random_split_results_list <- lapply(all_sim_results, function(x) x$random_split_results)
all_random_split_results <- do.call(rbind, random_split_results_list)

# 提取空间分块分割结果
block_split_results_list <- lapply(all_sim_results, function(x) x$block_split_results)
all_block_split_results <- do.call(rbind, block_split_results_list)

# 查看结果结构
cat("\n=== 完整数据集结果（敏感性分析） ===\n")
cat("总共", length(full_metrics_list), "次模拟\n")

cat("\n=== 随机分割结果（泛化能力分析） ===\n")
if (!is.null(all_random_split_results)) {
  cat("总行数:", nrow(all_random_split_results), "\n")
  cat("场景分布:\n")
  print(table(all_random_split_results$scenario))
  cat("分割方法分布:\n")
  print(table(all_random_split_results$split_method))
} else {
  cat("随机分割结果为空\n")
}

cat("\n=== 空间分块分割结果（泛化能力分析） ===\n")
if (!is.null(all_block_split_results)) {
  cat("总行数:", nrow(all_block_split_results), "\n")
  cat("场景分布:\n")
  print(table(all_block_split_results$scenario))
  cat("分割方法分布:\n")
  print(table(all_block_split_results$split_method))
} else {
  cat("空间分块分割结果为空\n")
}

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
output_dir <- "space_model_comparison_plots"
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
#-------切换数据集--------
# all_split_results_df <- all_random_split_results
all_split_results_df <- all_block_split_results
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
output_dir <- "space_model_comparison_plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("创建输出目录:", output_dir, "\n")
}

#-------修改后的方差排序图绘制逻辑（全部使用对数刻度）----------

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
    scenario_ranked$Variance_adj <- ifelse(scenario_ranked$Variance <= 0, 
                                           min_positive_variance * 0.1, 
                                           scenario_ranked$Variance)
    
    # 使用对数刻度版本
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
      ) +
      # 添加数值标签，使用合适的格式
      geom_text(
        aes(label = ifelse(Variance_adj >= 1000, 
                           sprintf("%.2e", Variance_adj), 
                           sprintf("%.3f", Variance_adj))), 
        hjust = -0.1, 
        size = 3.5, 
        color = "darkblue"
      )
    
    pdf_file <- file.path(output_dir, paste0("scenario_", scenario, "_delta_variance_log.pdf"))
    
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
cat("1. 按场景绘制方差排序柱状图（全部使用对数刻度）...\n")
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
  ranking_filename <- paste0("space_var_rank_", dataset_id, "_", scenario, ".csv")
  write.csv(ranking_results, ranking_filename, row.names = FALSE)
  cat("场景", scenario, "排名结果已保存到:", ranking_filename, "\n\n")
  
  # 6. 将当前场景的结果添加到总列表中
  all_ranking_results[[scenario]] <- ranking_results
}