# 加载必要的包
library(ape)
library(phytools)
library(nlme)
library(geiger)
library(rr2)
library(caper)
library(cluster)

# 计算各个指标的calculate_all_metrics函数，增加错误处理
calculate_all_metrics <- function(final_data, tree_subset, comp_data, 
                                  response_var = "reported_data", 
                                  predictor_var = "data") {
  
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
    mae = mae,
    median_ae = median_ae,
    
    # 相关性指标
    pearson_corr = pearson_corr,
    spearman_corr = spearman_corr,
    
    # 百分比误差指标
    mape = mape,
    smape = smape_val,
    cn_smape = cn_smape,
    
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
# 设置重复模拟次数
n_simulations <- 100  # 可以根据需要调整
set.seed(123)

# 1. 生成基础数据
n_species <- 128

# 生成系统发育树
tree <- rcoal(n_species)
tree$tip.label <- paste0("species", 1:n_species)

# 使用BM模型模拟真实性状值
real_trait <- rTraitCont(tree, model = "BM", sigma = 1)
# 创建存储所有模拟结果的列表
all_sim_results <- vector("list", n_simulations)

# 批量模拟函数
run_single_simulation <- function(sim_id) {
  cat("正在进行第", sim_id, "次模拟...\n")
  
  # 设置随机种子（保证每次模拟可重现）
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
  residual_signal <- rTraitCont(tree, model = "BM", sigma = 0.5)
  pred_bad_uc2 <- real_trait + residual_signal
  base_data$data_bad_uc2 <- pred_bad_uc2
  
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
  
  # 6. 计算各场景指标
  scenarios <- c("Good", "UC1", "UC2", "UC3", "UC4")
  predictor_vars <- c("data_good", "data_bad_uc1", "data_bad_uc2", 
                      "data_bad_uc3", "data_bad_uc4")
  
  metrics_list <- list()
  for (i in 1:length(scenarios)) {
    metrics <- calculate_all_metrics(base_data, tree, base_data,
                                     predictor_var = predictor_vars[i])
    metrics_list[[scenarios[i]]] <- metrics
  }
  
  return(metrics_list)
}

# 执行批量模拟
cat("开始批量模拟，总共", n_simulations, "次...\n")
start_time <- Sys.time()

# 使用lapply进行并行计算（可选）
# 如果需要并行计算，可以使用parallel包
all_sim_results <- lapply(1:n_simulations, run_single_simulation)

end_time <- Sys.time()
cat("批量模拟完成，耗时:", round(end_time - start_time, 2), "秒\n")

# 8. 计算平均指标
calculate_average_metrics <- function(all_results, n_sim) {
  scenarios <- c("Good", "UC1", "UC2", "UC3", "UC4")
  
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

# 计算平均指标
avg_metrics <- calculate_average_metrics(all_sim_results, n_simulations)

# 9. 创建最终比较表格
create_final_comparison_table <- function(avg_metrics) {
  scenarios <- c("Good", "UC1", "UC2", "UC3", "UC4")
  
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
  "ols_r2", "adjusted_r2", "pearson_corr", "spearman_corr", "cn_smape",
  "pic_r2", "pic_pearson", "pic_spearman", "whitened_r2", "resid_r2", "lik_r2"
)

lower_better_metrics <- c(
  "mse", "rmse", "mae", "median_ae", "mape", "smape",
  "pic_rmse", "pic_mae", "whitened_rmse", "whitened_mae"
)

# 筛选实际存在的指标
available_higher <- higher_better_metrics[higher_better_metrics %in% names(comparison_means)]
available_lower <- lower_better_metrics[lower_better_metrics %in% names(comparison_means)]

cat("基于", n_simulations, "次模拟的平均结果\n")
cat("可用的'越大越好'指标:", length(available_higher), "\n")
cat("可用的'越小越好'指标:", length(available_lower), "\n")

# 设置颜色方案
good_color <- "#c9cb05"  # 好模型
bad_colors <- c("#bab1d8","#8076b5","#9e7cba", "#ac5aa1")  # 紫色系 - 坏模型
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

