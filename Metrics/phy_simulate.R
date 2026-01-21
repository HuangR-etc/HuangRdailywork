# 加载必要的包
library(ape)
library(phytools)
library(nlme)
library(geiger)
library(rr2)
library(caper)
library(cluster)
packages <- c("ape", "phytools", "nlme","geiger","rr2","caper" , "cluster")
versions <- sapply(packages, function(x) as.character(packageVersion(x)))
versions
# 计算各个指标的calculate_all_metrics函数，增加错误处理
calculate_all_metrics <- function(final_data, tree_subset, comp_data, 
                                  response_var = "reported_data", 
                                  predictor_var = "data",
                                  delta = 0.01) {
  
  # 提取观测值和预测值
  y_true <- final_data[[response_var]]
  y_pred <- final_data[[predictor_var]]
  n <- length(y_true)
  
  # 检查预测值是否为常数
  y_pred_sd <- sd(y_pred, na.rm = TRUE)
  is_constant_pred <- y_pred_sd == 0 || is.na(y_pred_sd)
  
  # 1. 基础误差指标
  residuals <- y_true - y_pred
  
  # OLS回归 - 添加错误处理
  ols_r2 <- NA
  adjusted_r2 <- NA
  if (!is_constant_pred) {
    tryCatch({
      raw_model <- lm(as.formula(paste(response_var, "~", predictor_var)), data = final_data)
      ols_r2 <- summary(raw_model)$r.squared
      # 调整R² (p=1个预测变量)
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
  # 误差指标（这些应该可以正常计算）
  mse <- mean(residuals^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(residuals))
  median_ae <- median(abs(residuals))
  
  # 相关性指标 - 添加错误处理
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
  R_whole <- sum(y_true)/sum(y_pred)
  R_i <- safe_division(y_true, y_pred)  # 使用safe_division避免除以0
  
  # 汇总统计量
  R_arithmetic_mean <- mean(R_i, na.rm = TRUE)  # 算术平均值
  R_geometric_mean <- exp(mean(log(R_i), na.rm = TRUE))  # 几何平均值
  R_median <- median(R_i, na.rm = TRUE)  # 中位数比值
  
  
  # 2. 系统发育相关指标
  
  # 系统发育信号lambda
  lambda <- NA
  tryCatch({
    response_values <- y_true
    names(response_values) <- final_data$GCF
    lambda_result <- phytools::phylosig(tree_subset, response_values, method = "lambda", test = FALSE)
    lambda <- as.numeric(lambda_result$lambda)
  }, error = function(e) {
    lambda <<- NA
  })
  
  # PIC相关指标
  pic_obs <- NA
  pic_pred <- NA
  pic_pearson <- NA
  pic_spearman <- NA
  pic_rmse <- NA
  pic_mae <- NA
  pic_r2 <- NA
  
  if (!is_constant_pred) {
    tryCatch({
      pic_obs <- pic(y_true, tree_subset)
      pic_pred <- pic(y_pred, tree_subset)
      
      if (length(pic_obs) > 1 && length(pic_pred) > 1) {
        pic_pearson <- cor(pic_obs, pic_pred, method = "pearson")
        pic_spearman <- cor(pic_obs, pic_pred, method = "spearman")
        
        pic_residuals <- pic_obs - pic_pred
        pic_rmse <- sqrt(mean(pic_residuals^2))
        pic_mae <- mean(abs(pic_residuals))
        
        # PIC R² (无截距模型)
        pic_model <- lm(pic_obs ~ pic_pred - 1)
        pic_r2 <- summary(pic_model)$r.squared
      }
    }, error = function(e) {
      # 如果PIC计算失败，保持NA值
    })
  }
  
  # 系统发育白化指标
  calculate_whitened_metrics <- function(y, yhat, tip_labels, tree, lambda) {
    # 如果lambda是NA，返回NA
    if (is.na(lambda)) {
      return(list(
        whitened_r2 = NA,
        whitened_rmse = NA,
        whitened_mae = NA
      ))
    }
    
    tryCatch({
      # 辅助函数：根据lambda调整协方差矩阵
      adjust_vcv_by_lambda <- function(tree, lambda) {
        V_bm <- ape::vcv(tree)
        diag_vals <- diag(V_bm)
        V_lambda <- V_bm
        n <- nrow(V_bm)
        
        for (i in 1:(n-1)) {
          for (j in (i+1):n) {
            V_lambda[i, j] <- V_lambda[j, i] <- V_bm[i, j] * lambda
          }
        }
        eps <- 1e-8
        V_lambda <- V_lambda + diag(eps, n)
        return(V_lambda)
      }
      
      names(y) <- names(yhat) <- tip_labels
      V <- adjust_vcv_by_lambda(tree, lambda)
      V <- V[tip_labels, tip_labels, drop = FALSE]
      
      # Cholesky 白化
      eps <- 1e-8
      W <- solve(V + diag(eps, nrow(V)))
      L <- chol(W)
      
      # 白化变换
      y_star <- as.numeric(L %*% y)
      yhat_star <- as.numeric(L %*% yhat)
      
      # 在白化空间计算指标
      whitened_residuals <- y_star - yhat_star
      whitened_rmse <- sqrt(mean(whitened_residuals^2))
      whitened_mae <- mean(abs(whitened_residuals))
      
      # 白化R²
      ybar_star <- mean(y_star)
      RSS <- sum(whitened_residuals^2)
      TSS <- sum((y_star - ybar_star)^2)
      whitened_r2 <- 1 - RSS/TSS
      
      return(list(
        whitened_r2 = whitened_r2,
        whitened_rmse = whitened_rmse,
        whitened_mae = whitened_mae
      ))
    }, error = function(e) {
      return(list(
        whitened_r2 = NA,
        whitened_rmse = NA,
        whitened_mae = NA
      ))
    })
  }
  
  whitened_results <- calculate_whitened_metrics(
    y = y_true,
    yhat = y_pred,
    tip_labels = final_data$GCF,
    tree = tree_subset,
    lambda = lambda
  )
  
  # 3. Ives相关指标
  resid_r2 <- NA
  lik_r2 <- NA
  
  if (!is_constant_pred) {
    tryCatch({
      final_data$weights <- diag(ape::vcv(tree_subset))
      
      gls_model <- gls(as.formula(paste(response_var, "~", predictor_var)), 
                       data = final_data,
                       correlation = corBrownian(1, form = ~GCF, phy = tree_subset),
                       weights = varFixed(~weights),
                       method = "ML")
      
      # Ives提出的R2
      resid_r2 <- R2_resid(gls_model)
      lik_r2 <- R2_lik(gls_model)
    }, error = function(e) {
      # 如果GLS失败，保持NA值
    })
  }
  
  
  # 返回所有指标
  return(list(
    # 基础指标
    ols_r2 = ols_r2,
    adjusted_r2 = adjusted_r2,
    mse = mse,
    rmse = rmse,
    weighted_rmse = weighted_rmse,
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
    R_whole = R_whole,#整体比值
    R_ari = R_arithmetic_mean, #算数均值
    R_geo = R_geometric_mean, #几何均值
    R_median = R_median, #中位数
    
    # 系统发育指标
    lambda = lambda,
    pic_r2 = pic_r2,
    pic_pearson = pic_pearson,
    pic_spearman = pic_spearman,
    pic_rmse = pic_rmse,
    pic_mae = pic_mae,
    
    # 白化指标
    whitened_r2 = whitened_results$whitened_r2,
    whitened_rmse = whitened_results$whitened_rmse,
    whitened_mae = whitened_results$whitened_mae,
    
    # Ives指标
    resid_r2 = resid_r2,
    lik_r2 = lik_r2,
    
    # 添加一个标志，指示预测值是否为常数
    is_constant_pred = is_constant_pred
  ))
}
# 修剪树和数据的辅助函数
prune_tree_and_data <- function(tree, data) {
  # 确保数据中的GCF在树中存在
  common_species <- intersect(data$GCF, tree$tip.label)
  
  if (length(common_species) < 2) {
    stop("树和数据中共同的物种数量不足")
  }
  
  # 修剪树和数据
  tree_pruned <- keep.tip(tree, common_species)
  data_pruned <- data[data$GCF %in% common_species, ]
  
  # 确保数据顺序与树一致
  data_ordered <- data_pruned[match(tree_pruned$tip.label, data_pruned$GCF), ]
  
  # 创建比较数据对象
  comp_data <- comparative.data(tree_pruned, data_ordered, names.col = "GCF")
  
  return(list(
    tree = tree_pruned,
    data = data_ordered,
    comp_data = comp_data
  ))
}

# 随机分割分析函数（修改版）
perform_random_split_analysis <- function(final_data, tree, n_iterations = 50, train_ratio = 0.8, 
                                          response_var = "reported_data", predictor_var = "data",
                                          scenario_name = "Unknown", sim_id = 1) {
  
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
      n_total <- nrow(final_data)
      n_train <- round(n_total * train_ratio)
      n_test <- n_total - n_train
      
      train_indices <- sample(1:n_total, n_train, replace = FALSE)
      test_indices <- setdiff(1:n_total, train_indices)
      
      # 创建训练集和测试集
      train_data <- final_data[train_indices, ]
      test_data <- final_data[test_indices, ]
      
      # 对训练集进行计算
      prune_train <- prune_tree_and_data(tree, train_data)
      tree_train <- prune_train$tree
      train_data_ordered <- prune_train$data
      comp_data_train <- prune_train$comp_data
      
      metrics_train <- calculate_all_metrics(train_data_ordered, tree_train, comp_data_train,
                                             response_var = response_var, predictor_var = predictor_var)
      
      # 对测试集进行计算
      prune_test <- prune_tree_and_data(tree, test_data)
      tree_test <- prune_test$tree
      test_data_ordered <- prune_test$data
      comp_data_test <- prune_test$comp_data
      
      metrics_test <- calculate_all_metrics(test_data_ordered, tree_test, comp_data_test,
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
        lambda_train = metrics_train$lambda,
        pic_r2_train = metrics_train$pic_r2,
        pic_pearson_train = metrics_train$pic_pearson,
        pic_spearman_train = metrics_train$pic_spearman,
        pic_rmse_train = metrics_train$pic_rmse,
        pic_mae_train = metrics_train$pic_mae,
        whitened_r2_train = metrics_train$whitened_r2,
        whitened_rmse_train = metrics_train$whitened_rmse,
        whitened_mae_train = metrics_train$whitened_mae,
        resid_r2_train = metrics_train$resid_r2,
        lik_r2_train = metrics_train$lik_r2,
        
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
        lambda_test = metrics_test$lambda,
        pic_r2_test = metrics_test$pic_r2,
        pic_pearson_test = metrics_test$pic_pearson,
        pic_spearman_test = metrics_test$pic_spearman,
        pic_rmse_test = metrics_test$pic_rmse,
        pic_mae_test = metrics_test$pic_mae,
        whitened_r2_test = metrics_test$whitened_r2,
        whitened_rmse_test = metrics_test$whitened_rmse,
        whitened_mae_test = metrics_test$whitened_mae,
        resid_r2_test = metrics_test$resid_r2,
        lik_r2_test = metrics_test$lik_r2,
        
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
      successful_iterations <- successful_iterations + 1
      results[[successful_iterations]] <- result$result
    }
  }
  
  if (successful_iterations < n_iterations) {
    cat("警告: 场景", scenario_name, "模拟", sim_id, "只完成了", successful_iterations, "次成功分割\n")
  }
  
  return(do.call(rbind, results))
}

# 设置重复模拟次数
n_simulations <- 50
n_iterations <- 50
set.seed(123)

# 存储所有结果
all_sim_results <- list()
all_split_results <- list()

# 批量模拟函数（修改版，包含分割分析）
run_single_simulation <- function(sim_id) {
  cat("正在进行第", sim_id, "次模拟...\n")
  
  # 设置随机种子
  set.seed(123 + sim_id)
  
  # 1. 生成基础数据
  n_species <- 128
  
  # 生成系统发育树
  tree <- rcoal(n_species)
  tree$tip.label <- paste0("species", 1:n_species)
  
  # 使用BM模型模拟真实性状值
  real_trait <- rTraitCont(tree, model = "BM", sigma = 1)
  names(real_trait) <- tree$tip.label
  
  # 创建基础数据框
  base_data <- data.frame(
    GCF = tree$tip.label,
    reported_data = real_trait
  )
  
  # 生成好预测值（微小随机误差）
  pred_good <- real_trait + rnorm(n_species, mean = 0, sd = 0.1)
  base_data$data_good <- pred_good
  
  # 2. 场景1：系统发育聚集误差 (PH-UC1)
  distance_matrix <- cophenetic(tree)
  pam_result <- pam(distance_matrix, k = 4)
  folds <- pam_result$clustering
  names(folds) <- tree$tip.label
  
  biased_cluster <- sample(1:4, 1)
  bias_magnitude <- 2 * sd(real_trait)
  
  pred_bad_uc1 <- real_trait
  for (i in 1:n_species) {
    species_name <- tree$tip.label[i]
    if (folds[species_name] == biased_cluster) {
      pred_bad_uc1[i] <- pred_bad_uc1[i] + bias_magnitude
    }
  }
  base_data$data_bad_uc1 <- pred_bad_uc1
  
  # 3. 场景2：残差具有系统发育信号 (PH-UC2)
  # 3a. 强系统发育信号残差 (sigma=1)
  # 生成具有强系统发育信号的残差
  residual_strong_raw <- rTraitCont(tree, model = "BM", sigma = 1)
  # 标准化到目标方差
  target_variance <- 0.5^2  # 目标方差
  current_variance <- var(residual_strong_raw)
  scaling_factor <- sqrt(target_variance / current_variance)
  residual_strong <- residual_strong_raw * scaling_factor
  
  # 3b. 弱系统发育信号残差 (独立同分布噪声)
  # 生成独立同分布的残差
  residual_weak <- rnorm(n_species, mean = 0, sd = sqrt(target_variance))
  # 确保残差独立（无系统发育信号）
  pred_bad_uc2_strong <- real_trait + residual_strong
  pred_bad_uc2_weak <- real_trait + residual_weak
  
  base_data$data_bad_uc2_st <- pred_bad_uc2_strong
  base_data$data_bad_uc2_wk <- pred_bad_uc2_weak
  
  # 4. 场景3：远缘类群预测失败
  distance_matrix <- cophenetic(tree)
  hc <- hclust(as.dist(distance_matrix), method = "complete")
  outgroup_cut <- cutree(hc, k = 2)
  group_sizes <- table(outgroup_cut)
  outgroup_id <- which.min(group_sizes)
  outgroup_species <- which(outgroup_cut == outgroup_id)
  
  # 使用GLS估计系统发育均值
  phylo_intercept_model <- gls(reported_data ~ 1, 
                               data = base_data,
                               correlation = corBrownian(1, form = ~GCF, phy = tree),
                               method = "ML")
  phylogenetic_mean <- as.numeric(coef(phylo_intercept_model))
  
  pred_bad_uc3 <- real_trait
  pred_bad_uc3[outgroup_species] <- phylogenetic_mean
  base_data$data_bad_uc3 <- pred_bad_uc3
  
  # 5. 场景4：完全失败模型
  pred_bad_uc4 <- rep(phylogenetic_mean, n_species)
  names(pred_bad_uc4) <- tree$tip.label
  pred_bad_uc4 <- pred_bad_uc4 + rnorm(n_species, mean = 0, sd = 0.1)
  base_data$data_bad_uc4 <- pred_bad_uc4
  
  # 6. 计算各场景指标（完整数据集）
  scenarios <- c("Good", "UC1", "UC2_st","UC2_wk", "UC3", "UC4")
  predictor_vars <- c("data_good", "data_bad_uc1", "data_bad_uc2_st", "data_bad_uc2_wk",
                      "data_bad_uc3", "data_bad_uc4")
  
  metrics_list <- list()
  for (i in 1:length(scenarios)) {
    metrics <- calculate_all_metrics(base_data, tree, base_data,
                                     predictor_var = predictor_vars[i])
    metrics_list[[scenarios[i]]] <- metrics
  }
  
  # 7. 对每个场景进行分割分析
  split_results_list <- list()
  
  for (i in 1:length(scenarios)) {
    cat("  场景", scenarios[i], "分割分析...\n")
    
    split_results <- perform_random_split_analysis(
      final_data = base_data,
      tree = tree,
      n_iterations = n_iterations,
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
    base_data = base_data,  # 保存原始数据用于后续分析
    tree = tree
  ))
}

# 执行批量模拟
cat("开始批量模拟，总共", n_simulations, "次模拟，每次", n_iterations, "次分割...\n")
start_time <- Sys.time()

# 使用lapply进行串行计算
all_sim_results <- lapply(1:n_simulations, function(i) {
  result <- run_single_simulation(i)
  cat("完成第", i, "次模拟\n")
  return(result)
})

end_time <- Sys.time()
cat("批量模拟完成，耗时:", round(end_time - start_time, 2), "秒\n")

# 整理结果
# 提取完整数据集的指标
full_metrics_list <- lapply(all_sim_results, function(x) x$full_data_metrics)

# 提取分割分析结果
split_results_list <- lapply(all_sim_results, function(x) x$split_analysis_results)
all_split_results_df <- do.call(rbind, split_results_list)

#----------结果分析----------
# 8. 计算平均指标
calculate_average_metrics <- function(all_results, n_sim) {
  scenarios <- c("Good", "UC1", "UC2_st","UC2_wk", "UC3", "UC4")
  
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
  scenarios <- c("Good", "UC1", "UC2_st","UC2_wk", "UC3", "UC4")
  
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
output_dir <- "phy_model_comparison_plots"
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
      "UC3" = "#9e7cba", 
      "UC4" = "#ac5aa1"
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
  ranking_filename <- paste0("phy_var_rank_", dataset_id, "_", scenario, ".csv")
  write.csv(ranking_results, ranking_filename, row.names = FALSE)
  cat("场景", scenario, "排名结果已保存到:", ranking_filename, "\n\n")
  
  # 6. 将当前场景的结果添加到总列表中
  all_ranking_results[[scenario]] <- ranking_results
}