# 加载必要的包
library(ape)
library(phytools)
library(nlme)
library(geiger)
library(rr2)
library(caper)
# 设置随机种子保证结果可重现
set.seed(123)

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

# 查看基础数据
head(base_data)

# 2. 场景1：系统发育聚集误差 (PH-UC1) - 使用PAM聚类
cat("=== 场景1：系统发育聚集误差 ===\n")

# 使用PAM聚类方法
distance_matrix <- cophenetic(tree)
pam_result <- pam(distance_matrix, k = 4)
folds <- pam_result$clustering
names(folds) <- tree$tip.label

# 随机选择一个类群添加偏差
biased_cluster <- sample(1:4, 1)
#一个常数，等于二倍真实值标准差
bias_magnitude <- 2 * sd(real_trait)

pred_bad_uc1 <- real_trait
for (i in 1:n_species) {
  species_name <- tree$tip.label[i]
  if (folds[species_name] == biased_cluster) {
    pred_bad_uc1[i] <- pred_bad_uc1[i] + bias_magnitude
  }
}
base_data$data_bad_uc1 <- pred_bad_uc1

# # 绘制着色发育树
# cluster_colors <- c("red", "blue", "green", "purple")
# tip_colors <- cluster_colors[folds[tree$tip.label]]
# 
# par(mfrow = c(1, 2))
# plot(tree, show.tip.label = FALSE, main = "PAM聚类结果")
# tiplabels(pch = 16, col = tip_colors)
# legend("topleft", legend = paste("类群", 1:4), 
#        col = cluster_colors, pch = 16, cex = 0.8)
# 
# plot(tree, show.tip.label = FALSE, main = paste("PH-UC1: 类群", biased_cluster, "被添加偏差"))
# biased_tip_colors <- ifelse(folds[tree$tip.label] == biased_cluster, "red", "black")
# tiplabels(pch = 16, col = biased_tip_colors)
# legend("topleft", legend = c("被添加偏差的类群", "其他类群"), 
#        col = c("red", "black"), pch = 16, cex = 0.8)
# par(mfrow = c(1, 1))
head(base_data)
# 3. 场景2：残差具有系统发育信号 (PH-UC2)
cat("=== 场景2：残差具有系统发育信号 ===\n")

residual_signal <- rTraitCont(tree, model = "BM", sigma = 0.5)
pred_bad_uc2 <- real_trait + residual_signal
base_data$data_bad_uc2 <- pred_bad_uc2

# 4. 场景3：近期分化类群预测失败 (PH-UC3)
cat("=== 场景3：近期分化类群预测失败 ===\n")

# 计算节点深度并分类物种
node_depths <- node.depth.edgelength(tree)
tip_depths <- node_depths[1:n_species]
depth_threshold <- quantile(tip_depths, 0.7)

basal_species <- which(tip_depths <= depth_threshold)
derived_species <- which(tip_depths > depth_threshold)

cat("基部物种数量:", length(basal_species), "\n")
cat("衍生物种数量:", length(derived_species), "\n")

# 计算系统发育均值（使用GLS模型的截距）
# 方法1：使用ape::root函数获取根节点状态
# 方法2：使用GLS模型估计系统发育均值（更稳健）

# 使用GLS估计系统发育均值
# 创建一个只有截距的GLS模型
phylo_intercept_model <- gls(reported_data ~ 1, 
                             data = base_data,
                             correlation = corBrownian(1, form = ~GCF,phy = tree),
                             method = "ML")

# 获取系统发育均值（截距）
phylogenetic_mean <- as.numeric(coef(phylo_intercept_model))
# 生成预测值
pred_bad_uc3 <- real_trait
pred_bad_uc3[derived_species] <- mean(real_trait)
base_data$data_bad_uc3 <- pred_bad_uc3

# # 可视化
# par(mfrow = c(1, 2))
# plot(tree, show.tip.label = FALSE, main = "基部vs衍生物种分类")
# tiplabels(pch = 16, col = ifelse(tip_depths > depth_threshold, "red", "blue"))
# legend("topleft", legend = c("基部物种", "衍生物种"), col = c("blue", "red"), pch = 16)
# 
# plot(tree, show.tip.label = FALSE, main = "场景3: 衍生物种预测失败")
# tiplabels(pch = 16, col = ifelse(tip_depths > depth_threshold, "red", "black"))
# legend("topleft", legend = c("预测失败", "预测准确"), col = c("red", "black"), pch = 16)
# par(mfrow = c(1, 1))

# 5. 场景4：远缘类群预测失败 (PH-UC4)
cat("=== 场景4：远缘类群预测失败 ===\n")
distance_matrix <- cophenetic(tree)
hc <- hclust(as.dist(distance_matrix), method = "complete")

# 切割成2个主要簇，选择较小的作为外类群
outgroup_cut <- cutree(hc, k = 2)
group_sizes <- table(outgroup_cut)
outgroup_id <- which.min(group_sizes)  # 选择较小的簇作为外类群

outgroup_species <- which(outgroup_cut == outgroup_id)
ingroup_species <- which(outgroup_cut != outgroup_id)

cat("外类群物种数量:", length(outgroup_species), "\n")
cat("内类群物种数量:", length(ingroup_species), "\n")

# 生成预测值：外类群预测失败（使用系统发育均值）
pred_bad_uc4 <- real_trait
pred_bad_uc4[outgroup_species] <- phylogenetic_mean  # 使用系统发育均值替换普通均值
base_data$data_bad_uc4 <- pred_bad_uc4

# # 场景4可视化
# par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
# 
# # 左图：显示内外类群划分
# plot(tree, show.tip.label = FALSE, main = "内类群 vs 外类群")
# tip_colors_uc4 <- ifelse(outgroup_cut == outgroup_id, "orange", "purple")
# tiplabels(pch = 16, col = tip_colors_uc4)
# legend("topleft", legend = c("内类群", "外类群"), 
#        col = c("purple", "orange"), pch = 16, cex = 0.8)
# 
# # 右图：显示预测失败情况
# plot(tree, show.tip.label = FALSE, main = "场景4: 外类群预测失败")
# tip_colors_fail <- ifelse(outgroup_cut == outgroup_id, "red", "black")
# tiplabels(pch = 16, col = tip_colors_fail)
# legend("topleft", legend = c("预测失败", "预测准确"), 
#        col = c("red", "black"), pch = 16, cex = 0.8)
# 
# # 重置图形参数
# par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# 6. 场景5：完全失败模型 (PH-UC4)
cat("=== 场景5：完全失败模型 ===\n")

# 生成预测值：所有物种都预测为系统发育均值
pred_bad_uc5 <- rep(phylogenetic_mean, n_species)  # 使用系统发育均值替换普通均值
names(pred_bad_uc5) <- tree$tip.label
base_data$data_bad_uc5 <- pred_bad_uc5
# 查看完整数据
head(base_data)

# 7. 计算各场景的指标
# 修改后的calculate_all_metrics函数，增加错误处理
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
cat("开始计算各场景的指标...\n")

# 好预测模型的指标
cat("\n--- 好预测模型指标 ---\n")
metrics_good <- calculate_all_metrics(base_data, tree, base_data, 
                                      predictor_var = "data_good")

# 场景1指标
cat("\n--- 场景1：系统发育聚集误差指标 ---\n")
metrics_uc1 <- calculate_all_metrics(base_data, tree, base_data,
                                     predictor_var = "data_bad_uc1")

# 场景2指标  
cat("\n--- 场景2：残差具有系统发育信号指标 ---\n")
metrics_uc2 <- calculate_all_metrics(base_data, tree, base_data,
                                     predictor_var = "data_bad_uc2")

# 场景3指标
cat("\n--- 场景3：近期分化类群预测失败指标 ---\n")
metrics_uc3 <- calculate_all_metrics(base_data, tree, base_data,
                                     predictor_var = "data_bad_uc3")

# 场景4指标
cat("\n--- 场景4：远缘类群预测失败指标 ---\n")
metrics_uc4 <- calculate_all_metrics(base_data, tree, base_data,
                                     predictor_var = "data_bad_uc4")

# 场景5指标
cat("\n--- 场景5：完全失败模型指标 ---\n")
metrics_uc5 <- calculate_all_metrics(base_data, tree, base_data,
                                     predictor_var = "data_bad_uc5")

#-------8. 创建完整的指标比较表格-------------
cat("\n=== 各场景指标比较 ===\n")
scenarios <- c("Good", "UC1", "UC2", "UC3", "UC4", "UC5")
all_metrics <- list(metrics_good, metrics_uc1, metrics_uc2, metrics_uc3, metrics_uc4, metrics_uc5)

# 提取所有指标名称（排除NULL值）
get_all_metric_names <- function(metrics_list) {
  all_names <- unique(unlist(lapply(metrics_list, names)))
  # 过滤掉可能的NULL或NA值
  valid_names <- character()
  for (name in all_names) {
    if (!is.null(metrics_list[[1]][[name]]) && !is.na(metrics_list[[1]][[name]])) {
      valid_names <- c(valid_names, name)
    }
  }
  return(valid_names)
}

metric_names <- get_all_metric_names(all_metrics)
cat("可用的指标数量:", length(metric_names), "\n")
cat("指标名称:", paste(metric_names, collapse = ", "), "\n")

# 创建完整的比较表格
# 创建完整的比较表格
create_comparison_table <- function(metrics_list, metric_names, scenarios) {
  result <- data.frame(Scenario = scenarios)
  
  for (metric in metric_names) {
    metric_values <- sapply(metrics_list, function(x) {
      value <- x[[metric]]
      if (is.null(value) || is.na(value)) NA else value
    })
    result[[metric]] <- metric_values
  }
  
  return(result)
}

comparison_df <- create_comparison_table(all_metrics, metric_names, scenarios)

# 修正：只对数值列进行四舍五入
print_rounded_comparison <- function(df, digits = 4) {
  # 复制数据框
  result <- df
  
  # 对数值列进行四舍五入
  numeric_cols <- sapply(result, is.numeric)
  result[numeric_cols] <- round(result[numeric_cols], digits)
  
  return(result)
}

# 打印四舍五入后的比较表格
cat("=== 各场景指标比较（四舍五入到4位小数） ===\n")
print(print_rounded_comparison(comparison_df, 4))

# 为了后续分析，我们还需要一个只包含数值指标的版本
# 移除非数值列和常数预测标志
numeric_metric_names <- setdiff(metric_names, c("Scenario", "is_constant_pred"))
comparison_numeric <- comparison_df[, c("Scenario", numeric_metric_names)]

#-------9. 改进的可视化方案----------
# 9. 改进的可视化方案 - 按指标类型分开绘图

# 定义指标分类
higher_better_metrics <- c(
  "ols_r2", "adjusted_r2", "pearson_corr", "spearman_corr", "cn_smape",
  "pic_r2", "pic_pearson", "pic_spearman", "whitened_r2", "resid_r2", "lik_r2"
)

lower_better_metrics <- c(
  "mse", "rmse", "mae", "median_ae", "mape", "smape",
  "pic_rmse", "pic_mae", "whitened_rmse", "whitened_mae"
)

# 检查哪些指标在数据中实际存在
available_higher <- higher_better_metrics[higher_better_metrics %in% names(comparison_numeric)]
available_lower <- lower_better_metrics[lower_better_metrics %in% names(comparison_numeric)]

cat("可用的'越大越好'指标:", length(available_higher), "\n")
cat("可用的'越小越好'指标:", length(available_lower), "\n")

# 设置颜色方案
good_color <- "#2E8B57"  # 海绿色 - 好模型
bad_colors <- c("#E74C3C", "#C0392B", "#A93226", "#922B21", "#7B241C")  # 红色系 - 坏模型
scenario_colors <- c(good_color, bad_colors)

# 9.1 绘制"越大越好"指标对比图
cat("\n--- 绘制'越大越好'指标对比图 ---\n")

# 计算面板布局
n_higher <- length(available_higher)
n_cols <- 4  # 每行4个图
n_rows <- ceiling(n_higher / n_cols)

# 设置图形参数
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), mgp = c(2, 0.8, 0), cex.main = 1.1)

for (metric in available_higher) {
  # 提取当前指标的数据
  values <- comparison_numeric[[metric]]
  scenarios <- comparison_numeric$Scenario
  
  # 跳过全为NA的指标
  if (all(is.na(values))) {
    plot(1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
         main = paste(metric, "\n(无数据)"), bty = "n")
    next
  }
  
  # 创建条形图
  bar_plot <- barplot(values, names.arg = scenarios, 
                      col = scenario_colors,
                      main = metric,
                      ylab = "指标值",
                      ylim = c(0, max(values, na.rm = TRUE) * 1.1),
                      border = "white",
                      las = 2)
  
  # 添加数值标签
  text(x = bar_plot, y = values, 
       label = round(values, 3), 
       pos = 3, cex = 0.7, col = "darkblue", font = 2)
  
  # 添加参考线（好模型的值）
  good_value <- values[scenarios == "Good"]
  if (!is.na(good_value)) {
    abline(h = good_value, col = "blue", lty = 2, lwd = 1.5)
    text(x = max(bar_plot) * 0.8, y = good_value * 1.05, 
         labels = "好模型参考线", cex = 0.6, col = "blue")
  }
  
  # 添加指标性能说明
  mtext("↑ 值越大越好", side = 3, line = 0, cex = 0.6, col = "darkgreen")
}

# 添加总标题
mtext("模型评估指标比较 - '越大越好'型指标", 
      side = 3, outer = TRUE, line = -1.5, cex = 1.3, font = 2)

# 9.2 绘制"越小越好"指标对比图
cat("\n--- 绘制'越小越好'指标对比图 ---\n")

# 计算面板布局
n_lower <- length(available_lower)
n_cols <- 4
n_rows <- ceiling(n_lower / n_cols)

# 设置图形参数
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), mgp = c(2, 0.8, 0), cex.main = 1.1)

for (metric in available_lower) {
  # 提取当前指标的数据
  values <- comparison_numeric[[metric]]
  scenarios <- comparison_numeric$Scenario
  
  # 跳过全为NA的指标
  if (all(is.na(values))) {
    plot(1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
         main = paste(metric, "\n(无数据)"), bty = "n")
    next
  }
  
  # 创建条形图
  bar_plot <- barplot(values, names.arg = scenarios, 
                      col = scenario_colors,
                      main = metric,
                      ylab = "指标值",
                      ylim = c(0, max(values, na.rm = TRUE) * 1.1),
                      border = "white",
                      las = 2)
  
  # 添加数值标签
  text(x = bar_plot, y = values, 
       label = round(values, 3), 
       pos = 3, cex = 0.7, col = "darkblue", font = 2)
  
  # 添加参考线（好模型的值）
  good_value <- values[scenarios == "Good"]
  if (!is.na(good_value)) {
    abline(h = good_value, col = "blue", lty = 2, lwd = 1.5)
    text(x = max(bar_plot) * 0.8, y = good_value * 1.05, 
         labels = "好模型参考线", cex = 0.6, col = "blue")
  }
  
  # 添加指标性能说明
  mtext("↓ 值越小越好", side = 3, line = 0, cex = 0.6, col = "darkred")
}

# 添加总标题
mtext("模型评估指标比较 - '越小越好'型指标", 
      side = 3, outer = TRUE, line = -1.5, cex = 1.3, font = 2)

# 重置图形参数
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# 9.3 创建总结性表格，显示各指标区分好坏模型的能力
cat("\n=== 各指标区分能力分析 ===\n")

analyze_discrimination_power <- function(comparison_df, higher_metrics, lower_metrics) {
  results <- data.frame(
    Metric = character(),
    Type = character(),
    Good_Value = numeric(),
    Bad_Avg_Value = numeric(),
    Discrimination_Ratio = numeric(),
    Discrimination_Power = character(),
    stringsAsFactors = FALSE
  )
  
  # 分析"越大越好"型指标
  for (metric in higher_metrics) {
    if (metric %in% names(comparison_df)) {
      values <- comparison_df[[metric]]
      good_value <- values[comparison_df$Scenario == "Good"]
      bad_values <- values[comparison_df$Scenario != "Good"]
      bad_avg <- mean(bad_values, na.rm = TRUE)
      
      # 计算区分度比率
      if (!is.na(good_value) && !is.na(bad_avg) && bad_avg != 0) {
        ratio <- good_value / bad_avg
      } else {
        ratio <- NA
      }
      
      # 评估区分能力
      if (!is.na(ratio)) {
        if (ratio > 1.5) power <- "强"
        else if (ratio > 1.1) power <- "中等"
        else if (ratio > 1.0) power <- "弱"
        else power <- "无"
      } else {
        power <- "无法计算"
      }
      
      results <- rbind(results, data.frame(
        Metric = metric,
        Type = "越大越好",
        Good_Value = round(good_value, 4),
        Bad_Avg_Value = round(bad_avg, 4),
        Discrimination_Ratio = round(ratio, 4),
        Discrimination_Power = power
      ))
    }
  }
  
  # 分析"越小越好"型指标
  for (metric in lower_metrics) {
    if (metric %in% names(comparison_df)) {
      values <- comparison_df[[metric]]
      good_value <- values[comparison_df$Scenario == "Good"]
      bad_values <- values[comparison_df$Scenario != "Good"]
      bad_avg <- mean(bad_values, na.rm = TRUE)
      
      # 计算区分度比率
      if (!is.na(good_value) && !is.na(bad_avg) && good_value != 0) {
        ratio <- bad_avg / good_value
      } else {
        ratio <- NA
      }
      
      # 评估区分能力
      if (!is.na(ratio)) {
        if (ratio > 1.5) power <- "强"
        else if (ratio > 1.1) power <- "中等"
        else if (ratio > 1.0) power <- "弱"
        else power <- "无"
      } else {
        power <- "无法计算"
      }
      
      results <- rbind(results, data.frame(
        Metric = metric,
        Type = "越小越好",
        Good_Value = round(good_value, 4),
        Bad_Avg_Value = round(bad_avg, 4),
        Discrimination_Ratio = round(ratio, 4),
        Discrimination_Power = power
      ))
    }
  }
  
  return(results)
}

# 生成区分能力分析表格
discrimination_analysis <- analyze_discrimination_power(
  comparison_df, 
  available_higher, 
  available_lower
)

print(discrimination_analysis)
