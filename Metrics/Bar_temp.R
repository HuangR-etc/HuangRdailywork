#函数封装版
library(dplyr)
library(stringr)
library(ape)
library(caper)
library(nlme)
library(rr2)
setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/get_obs_pred/temp")
#每次分析读取特定数据文件
filename <- "Barnum2024_temperature_clean.csv"
#------通用数据读取和分析代码----------
result_file <- "r2_comparison_results.csv"
df <- read.csv(filename)
tree <- read.tree("D:/pine/paper_redo/GTDB_tree/bac120_r226.tree")
# 读取元数据，只保留需要的列
cols_to_keep <- c("accession", "ncbi_organism_name", "gtdb_taxonomy", 
                  "gtdb_representative", "gtdb_genome_representative",
                  "ncbi_refseq_category", "ncbi_assembly_level")
gtdb_meta <- read.delim("D:/pine/paper_redo/GTDB_tree/bac120_metadata_r226.tsv", 
                        sep = "\t", 
                        header = TRUE)[, cols_to_keep]


#根据不同数据集定制的数据和树匹配代码
#-------数据集2：Barnum2024---------
# 提取数据集中的GCA编号（去掉"GCA_"前缀）
df$GCA_number <- gsub("GCA_", "", df$genbank_accession)

# 提取树tip.label中的GCF编号（去掉"RS_GCF_"前缀和".1"后缀）
tree_tip_numbers <- gsub("RS_GCF_", "", tree$tip.label)
tree_tip_numbers <- gsub("\\.\\d+$", "", tree_tip_numbers)  # 移除末尾的版本号

# 创建匹配的数据框
tree_tip_df <- data.frame(
  tip_label = tree$tip.label,
  GCF_number = tree_tip_numbers
)

# 将数据集与树tip进行匹配
matched_data <- merge(df, tree_tip_df, 
                      by.x = "GCA_number", by.y = "GCF_number", 
                      all.x = FALSE, all.y = FALSE)

# 重命名列以符合后续代码要求

final_data <- matched_data
final_data$GCF <- final_data$tip_label
rownames(final_data) <- final_data$GCF
#计算整体数据的系统发育信号(观察值列)
# response_values <- final_data$reported_temperature_optimum
# names(response_values) <- final_data$GCF
# tree_subset <- prune_tree_and_data(tree, final_data)$tree
# lambda_result <- phytools::phylosig(tree_subset, response_values, method = "lambda", test = FALSE)
# lambda <- as.numeric(lambda_result$lambda)
# 检查匹配结果
cat("原始数据集物种数:", nrow(df), "\n")
cat("匹配成功的物种数:", nrow(final_data), "\n")
cat("匹配的物种:", final_data$genbank_accession, "\n")
#------函数------
# 修剪树和数据匹配的通用函数
prune_tree_and_data <- function(tree, final_data) {
  # 修剪树
  tree_subset <- keep.tip(tree, final_data$GCF)
  tree_subset$node.label <- NULL  # 移除内部节点
  
  # 确保数据顺序与树一致
  final_data <- final_data[tree_subset$tip.label, ]
  
  # 创建系统发育比较对象
  comp_data <- comparative.data(phy = tree_subset, 
                                data = final_data, 
                                names.col = "GCF", 
                                vcv = TRUE, 
                                warn.dropped = TRUE)
  
  # 检查匹配结果
  cat("原始数据集物种数:", nrow(final_data), "\n")
  cat("修剪后的树包含", length(tree_subset$tip.label), "个物种\n")
  if (length(comp_data$dropped$tips) > 0) {
    cat("被丢弃的物种:", comp_data$dropped$tips, "\n")
  }
  if (length(comp_data$dropped$unmatched.rows) > 0) {
    cat("被丢弃的数据:", comp_data$dropped$unmatched.rows, "\n")
  }
  
  return(list(tree = tree_subset, data = final_data, comp_data = comp_data))
}

# 修改后的函数：计算所有指标
calculate_all_metrics <- function(final_data, tree_subset, comp_data, 
                                  response_var = "reported_temperature_optimum", 
                                  predictor_var = "temperature_optimum") {
  
  # 提取观测值和预测值
  y_true <- final_data[[response_var]]
  y_pred <- final_data[[predictor_var]]
  n <- length(y_true)
  
  # 1. 基础误差指标
  residuals <- y_true - y_pred
  
  # OLS回归
  raw_model <- lm(as.formula(paste(response_var, "~", predictor_var)), data = final_data)
  ols_r2 <- summary(raw_model)$r.squared
  
  # 调整R² (p=1个预测变量)
  p <- 1
  adjusted_r2 <- 1 - (1 - ols_r2) * (n - 1) / (n - p - 1)
  
  # 误差指标
  mse <- mean(residuals^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(residuals))
  median_ae <- median(abs(residuals))
  
  # 相关性指标
  pearson_corr <- cor(y_true, y_pred, method = "pearson")
  spearman_corr <- cor(y_true, y_pred, method = "spearman")
  
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
  response_values <- y_true
  names(response_values) <- final_data$GCF
  lambda_result <- phytools::phylosig(tree_subset, response_values, method = "lambda", test = FALSE)
  lambda <- as.numeric(lambda_result$lambda)
  
  # PIC相关指标
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
  } else {
    pic_pearson <- NA
    pic_spearman <- NA
    pic_rmse <- NA
    pic_mae <- NA
    pic_r2 <- NA
  }
  
  # 系统发育白化指标
  calculate_whitened_metrics <- function(y, yhat, tip_labels, tree, lambda) {
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
  }
  
  whitened_results <- calculate_whitened_metrics(
    y = y_true,
    yhat = y_pred,
    tip_labels = final_data$GCF,
    tree = tree_subset,
    lambda = lambda
  )
  
  # 3. GLS相关指标
  final_data$weights <- diag(ape::vcv(tree_subset))
  gls_model <- gls(as.formula(paste(response_var, "~", predictor_var)), 
                   data = final_data,
                   correlation = corBrownian(1, form = ~GCF, phy = tree_subset),
                   weights = varFixed(~weights),
                   method = "ML")
  
  # Ives提出的R2
  resid_r2 <- R2_resid(gls_model)
  lik_r2 <- R2_lik(gls_model)
  
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
    
    # GLS指标
    resid_r2 = resid_r2,
    lik_r2 = lik_r2
  ))
}
# # 测试确认无误
# result <- calculate_all_metrics(final_data, tree_subset, comp_data,
#                            response_var = "reported_temperature_optimum",
#                            predictor_var = "temperature_optimum")

#--------抽样准备，通用函数---------
#------新增：创建随机划分结果文件头----------
# 更新结果文件头
random_split_file <- "random_split_results_1118.csv"
if (!file.exists(random_split_file)) {
  random_split_header <- data.frame(
    filename = character(),
    iteration = integer(),
    random_seed = integer(),
    train_size = integer(),
    test_size = integer(),
    
    # 训练集指标
    ols_r2_train = numeric(),
    adjusted_r2_train = numeric(),
    mse_train = numeric(),
    rmse_train = numeric(),
    mae_train = numeric(),
    median_ae_train = numeric(),
    pearson_corr_train = numeric(),
    spearman_corr_train = numeric(),
    mape_train = numeric(),
    smape_train = numeric(),
    cn_smape_train = numeric(),
    lambda_train = numeric(),
    pic_r2_train = numeric(),
    pic_pearson_train = numeric(),
    pic_spearman_train = numeric(),
    pic_rmse_train = numeric(),
    pic_mae_train = numeric(),
    whitened_r2_train = numeric(),
    whitened_rmse_train = numeric(),
    whitened_mae_train = numeric(),
    resid_r2_train = numeric(),
    lik_r2_train = numeric(),
    
    # 测试集指标
    ols_r2_test = numeric(),
    adjusted_r2_test = numeric(),
    mse_test = numeric(),
    rmse_test = numeric(),
    mae_test = numeric(),
    median_ae_test = numeric(),
    pearson_corr_test = numeric(),
    spearman_corr_test = numeric(),
    mape_test = numeric(),
    smape_test = numeric(),
    cn_smape_test = numeric(),
    lambda_test = numeric(),
    pic_r2_test = numeric(),
    pic_pearson_test = numeric(),
    pic_spearman_test = numeric(),
    pic_rmse_test = numeric(),
    pic_mae_test = numeric(),
    whitened_r2_test = numeric(),
    whitened_rmse_test = numeric(),
    whitened_mae_test = numeric(),
    resid_r2_test = numeric(),
    lik_r2_test = numeric(),
    
    stringsAsFactors = FALSE
  )
  write.csv(random_split_header, random_split_file, row.names = FALSE)
}
# 修改随机划分分析函数
perform_random_split_analysis <- function(final_data, tree, n_iterations = 50, train_ratio = 0.5) {
  
  results <- list()
  successful_iterations <- 0
  iteration_count <- 0
  
  while (successful_iterations < n_iterations && iteration_count < n_iterations * 2) {
    iteration_count <- iteration_count + 1
    cat("正在进行第", iteration_count, "次尝试，成功次数:", successful_iterations, "/", n_iterations, "...\n")
    
    # 设置随机种子
    current_seed <- as.integer(Sys.time()) + iteration_count
    set.seed(current_seed)
    
    # 使用tryCatch进行错误处理
    result <- tryCatch({
      # 随机划分数据
      n_total <- nrow(final_data)
      
      # 特殊处理：当train_ratio为0.5时，确保两个数据集大小完全相等
      if (train_ratio == 0.5) {
        n_half <- floor(n_total / 2)
        n_train <- n_half
        n_test <- n_half
        n_used <- n_train + n_test
        
        used_indices <- sample(1:n_total, n_used, replace = FALSE)
        train_indices <- sample(used_indices, n_train, replace = FALSE)
        test_indices <- setdiff(used_indices, train_indices)
        
        cat("比例设置为0.5，使用对等划分: 训练集=", n_train, "测试集=", n_test)
        if (n_used < n_total) {
          cat("，舍弃了", n_total - n_used, "个样本\n")
        } else {
          cat("\n")
        }
      } else {
        n_train <- round(n_total * train_ratio)
        n_test <- n_total - n_train
        n_used <- n_total
        
        train_indices <- sample(1:n_total, n_train, replace = FALSE)
        test_indices <- setdiff(1:n_total, train_indices)
      }
      
      # 创建训练集和测试集
      train_data <- final_data[train_indices, ]
      test_data <- final_data[test_indices, ]
      
      # 对训练集进行计算
      prune_train <- prune_tree_and_data(tree, train_data)
      tree_train <- prune_train$tree
      train_data_ordered <- prune_train$data
      comp_data_train <- prune_train$comp_data
      
      metrics_train <- calculate_all_metrics(train_data_ordered, tree_train, comp_data_train)
      
      # 对测试集进行计算
      prune_test <- prune_tree_and_data(tree, test_data)
      tree_test <- prune_test$tree
      test_data_ordered <- prune_test$data
      comp_data_test <- prune_test$comp_data
      
      metrics_test <- calculate_all_metrics(test_data_ordered, tree_test, comp_data_test)
      
      # 构建结果行
      result_row <- data.frame(
        filename = filename,
        iteration = successful_iterations + 1,
        random_seed = current_seed,
        train_size = n_train,
        test_size = n_test,
        
        # 训练集指标
        ols_r2_train = metrics_train$ols_r2,
        adjusted_r2_train = metrics_train$adjusted_r2,
        mse_train = metrics_train$mse,
        rmse_train = metrics_train$rmse,
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
        
        # 测试集指标
        ols_r2_test = metrics_test$ols_r2,
        adjusted_r2_test = metrics_test$adjusted_r2,
        mse_test = metrics_test$mse,
        rmse_test = metrics_test$rmse,
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
        lik_r2_test = metrics_test$lik_r2
      )
      
      # 返回成功的结果
      list(success = TRUE, result = result_row)
      
    }, error = function(e) {
      cat("第", iteration_count, "次尝试失败，错误信息:", e$message, "\n")
      list(success = FALSE, result = NULL)
    })
    
    # 如果成功，保存结果
    if (result$success) {
      successful_iterations <- successful_iterations + 1
      results[[successful_iterations]] <- result$result
      
      # 每次成功迭代后立即写入文件
      write.table(result$result, random_split_file, 
                  sep = ",", col.names = FALSE, row.names = FALSE, 
                  append = TRUE)
      
      cat("第", successful_iterations, "次成功划分完成\n")
    }
  }
  
  if (successful_iterations < n_iterations) {
    cat("警告: 只完成了", successful_iterations, "次成功迭代，目标为", n_iterations, "次\n")
  }
  
  return(do.call(rbind, results))
}
#--------划分1：对整体数据集进行82划分-----------
prune <- prune_tree_and_data(tree, final_data)
tree_subset <- prune$tree
final_data <- prune$data
comp_data <- prune$comp_data
# 计算完整数据的R2
# full_r2 <- calculate_all_r2(final_data_ordered, tree_subset, comp_data)
cat("开始50次随机划分分析...\n")
random_results <- perform_random_split_analysis(final_data, tree, n_iterations = 100,train_ratio = 0.8)

cat("所有分析完成！结果已保存至:", random_split_file, "\n")
#--------划分2：筛选出原训练集，进行82划分-----------
# library(dplyr)
# train_data <- final_data %>%
#   filter(partition_temperature == "train")
# prune <- prune_tree_and_data(tree, train_data)
# tree_subset <- prune$tree
# train_data <- prune$data
# comp_data <- prune$comp_data
# # 计算完整数据的R2
# # full_r2 <- calculate_all_r2(final_data_ordered, tree_subset, comp_data)
# cat("开始50次随机划分分析...\n")
# random_results <- perform_random_split_analysis(train_data, tree, n_iterations = 100,train_ratio = 0.8)
# 
# cat("所有分析完成！结果已保存至:", random_split_file, "\n")
#--------划分3：筛选出原测试集，进行55划分-----------
library(dplyr)
test_data <- final_data %>%
  filter(partition_temperature == "test")
prune <- prune_tree_and_data(tree, test_data)
tree_subset <- prune$tree
test_data <- prune$data
comp_data <- prune$comp_data
# 计算完整数据的R2
# full_r2 <- calculate_all_r2(final_data_ordered, tree_subset, comp_data)
cat("开始50次随机划分分析...\n")
random_results <- perform_random_split_analysis(test_data, tree, n_iterations = 100,train_ratio = 0.5)

cat("所有分析完成！结果已保存至:", random_split_file, "\n")
#--------划分4：对整体数据集进行55划分-----------
prune <- prune_tree_and_data(tree, final_data)
tree_subset <- prune$tree
final_data <- prune$data
comp_data <- prune$comp_data
# 计算完整数据的R2
# full_r2 <- calculate_all_r2(final_data_ordered, tree_subset, comp_data)
cat("开始50次随机划分分析...\n")
random_results <- perform_random_split_analysis(final_data, tree, n_iterations = 100,train_ratio = 0.5)

cat("所有分析完成！结果已保存至:", random_split_file, "\n")
#-------1118版本结果分析--------
result <- read.csv("random_split_results_1118.csv")
result1 <- result[1:100,]
# result1 <- result[101:200,]
# result1 <- result[201:300,]

# 定义各类指标
r2_metrics <- c("ols", "adjusted")
error_metrics <- c("mse", "rmse", "mae", "median_ae")
corr_metrics <- c("pearson", "spearman")
percentage_metrics <- c("mape","smape", "cn_smape")

# 计算delta列（test - train）
delta_columns <- c()

# 1. R²类指标
for(metric in r2_metrics) {
  train_col <- paste0(metric, "_r2_train")
  test_col <- paste0(metric, "_r2_test")
  delta_col <- paste0("delta_", metric, "_r2")
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    result1[[delta_col]] <- result1[[test_col]] - result1[[train_col]]
    delta_columns <- c(delta_columns, delta_col)
  }
}

# 2. 误差类指标
for(metric in error_metrics) {
  train_col <- paste0(metric, "_train")
  test_col <- paste0(metric, "_test")
  delta_col <- paste0("delta_", metric)
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    result1[[delta_col]] <- result1[[test_col]] - result1[[train_col]]
    delta_columns <- c(delta_columns, delta_col)
  }
}

# 3. 相关性指标 - 修复列名生成逻辑
for(metric in corr_metrics) {
  # 检查两种可能的列名格式
  train_col1 <- paste0(metric, "_corr_train")
  test_col1 <- paste0(metric, "_corr_test")
  train_col2 <- paste0(metric, "_train")
  test_col2 <- paste0(metric, "_test")
  
  if(train_col1 %in% colnames(result1) && test_col1 %in% colnames(result1)) {
    # 使用第一种格式：metric_corr_train/test
    delta_col <- paste0("delta_", metric, "_corr")
    result1[[delta_col]] <- result1[[test_col1]] - result1[[train_col1]]
    delta_columns <- c(delta_columns, delta_col)
  } else if(train_col2 %in% colnames(result1) && test_col2 %in% colnames(result1)) {
    # 使用第二种格式：metric_train/test
    delta_col <- paste0("delta_", metric)
    result1[[delta_col]] <- result1[[test_col2]] - result1[[train_col2]]
    delta_columns <- c(delta_columns, delta_col)
  }
}

# 4. 百分比误差指标
for(metric in percentage_metrics) {
  train_col <- paste0(metric, "_train")
  test_col <- paste0(metric, "_test")
  delta_col <- paste0("delta_", metric)
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    result1[[delta_col]] <- result1[[test_col]] - result1[[train_col]]
    delta_columns <- c(delta_columns, delta_col)
  }
}

# 新增指标1：相对变化率 (delta/train)
relative_columns <- c()

# R²类指标相对变化率
for(metric in r2_metrics) {
  train_col <- paste0(metric, "_r2_train")
  test_col <- paste0(metric, "_r2_test")
  relative_col <- paste0("relative_", metric, "_r2")
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    denominator <- ifelse(abs(result1[[train_col]]) < 1e-10, sign(result1[[train_col]]) * 1e-10, result1[[train_col]])
    result1[[relative_col]] <- (result1[[test_col]] - result1[[train_col]]) / abs(denominator)
    relative_columns <- c(relative_columns, relative_col)
  }
}

# 误差类指标相对变化率
for(metric in error_metrics) {
  train_col <- paste0(metric, "_train")
  test_col <- paste0(metric, "_test")
  relative_col <- paste0("relative_", metric)
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    denominator <- ifelse(abs(result1[[train_col]]) < 1e-10, sign(result1[[train_col]]) * 1e-10, result1[[train_col]])
    result1[[relative_col]] <- (result1[[test_col]] - result1[[train_col]]) / abs(denominator)
    relative_columns <- c(relative_columns, relative_col)
  }
}

# 相关性指标相对变化率 - 修复列名生成逻辑
for(metric in corr_metrics) {
  # 检查两种可能的列名格式
  train_col1 <- paste0(metric, "_corr_train")
  test_col1 <- paste0(metric, "_corr_test")
  train_col2 <- paste0(metric, "_train")
  test_col2 <- paste0(metric, "_test")
  
  if(train_col1 %in% colnames(result1) && test_col1 %in% colnames(result1)) {
    # 使用第一种格式
    relative_col <- paste0("relative_", metric, "_corr")
    denominator <- ifelse(abs(result1[[train_col1]]) < 1e-10, sign(result1[[train_col1]]) * 1e-10, result1[[train_col1]])
    result1[[relative_col]] <- (result1[[test_col1]] - result1[[train_col1]]) / abs(denominator)
    relative_columns <- c(relative_columns, relative_col)
  } else if(train_col2 %in% colnames(result1) && test_col2 %in% colnames(result1)) {
    # 使用第二种格式
    relative_col <- paste0("relative_", metric)
    denominator <- ifelse(abs(result1[[train_col2]]) < 1e-10, sign(result1[[train_col2]]) * 1e-10, result1[[train_col2]])
    result1[[relative_col]] <- (result1[[test_col2]] - result1[[train_col2]]) / abs(denominator)
    relative_columns <- c(relative_columns, relative_col)
  }
}

# 百分比误差指标相对变化率
for(metric in percentage_metrics) {
  train_col <- paste0(metric, "_train")
  test_col <- paste0(metric, "_test")
  relative_col <- paste0("relative_", metric)
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    denominator <- ifelse(abs(result1[[train_col]]) < 1e-10, sign(result1[[train_col]]) * 1e-10, result1[[train_col]])
    result1[[relative_col]] <- (result1[[test_col]] - result1[[train_col]]) / abs(denominator)
    relative_columns <- c(relative_columns, relative_col)
  }
}

# 新增指标2：标准化变化 (delta/sd)
standardized_columns <- c()

# R²类指标标准化变化
for(metric in r2_metrics) {
  train_col <- paste0(metric, "_r2_train")
  test_col <- paste0(metric, "_r2_test")
  standardized_col <- paste0("standardized_", metric, "_r2")
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    combined_values <- c(result1[[train_col]], result1[[test_col]])
    sd_value <- sd(combined_values, na.rm = TRUE)
    if(sd_value < 1e-10) sd_value <- 1e-10
    result1[[standardized_col]] <- (result1[[test_col]] - result1[[train_col]]) / sd_value
    standardized_columns <- c(standardized_columns, standardized_col)
  }
}

# 误差类指标标准化变化
for(metric in error_metrics) {
  train_col <- paste0(metric, "_train")
  test_col <- paste0(metric, "_test")
  standardized_col <- paste0("standardized_", metric)
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    combined_values <- c(result1[[train_col]], result1[[test_col]])
    sd_value <- sd(combined_values, na.rm = TRUE)
    if(sd_value < 1e-10) sd_value <- 1e-10
    result1[[standardized_col]] <- (result1[[test_col]] - result1[[train_col]]) / sd_value
    standardized_columns <- c(standardized_columns, standardized_col)
  }
}

# 相关性指标标准化变化 - 修复列名生成逻辑
for(metric in corr_metrics) {
  # 检查两种可能的列名格式
  train_col1 <- paste0(metric, "_corr_train")
  test_col1 <- paste0(metric, "_corr_test")
  train_col2 <- paste0(metric, "_train")
  test_col2 <- paste0(metric, "_test")
  
  if(train_col1 %in% colnames(result1) && test_col1 %in% colnames(result1)) {
    # 使用第一种格式
    standardized_col <- paste0("standardized_", metric, "_corr")
    combined_values <- c(result1[[train_col1]], result1[[test_col1]])
    sd_value <- sd(combined_values, na.rm = TRUE)
    if(sd_value < 1e-10) sd_value <- 1e-10
    result1[[standardized_col]] <- (result1[[test_col1]] - result1[[train_col1]]) / sd_value
    standardized_columns <- c(standardized_columns, standardized_col)
  } else if(train_col2 %in% colnames(result1) && test_col2 %in% colnames(result1)) {
    # 使用第二种格式
    standardized_col <- paste0("standardized_", metric)
    combined_values <- c(result1[[train_col2]], result1[[test_col2]])
    sd_value <- sd(combined_values, na.rm = TRUE)
    if(sd_value < 1e-10) sd_value <- 1e-10
    result1[[standardized_col]] <- (result1[[test_col2]] - result1[[train_col2]]) / sd_value
    standardized_columns <- c(standardized_columns, standardized_col)
  }
}

# 百分比误差指标标准化变化
for(metric in percentage_metrics) {
  train_col <- paste0(metric, "_train")
  test_col <- paste0(metric, "_test")
  standardized_col <- paste0("standardized_", metric)
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    combined_values <- c(result1[[train_col]], result1[[test_col]])
    sd_value <- sd(combined_values, na.rm = TRUE)
    if(sd_value < 1e-10) sd_value <- 1e-10
    result1[[standardized_col]] <- (result1[[test_col]] - result1[[train_col]]) / sd_value
    standardized_columns <- c(standardized_columns, standardized_col)
  }
}

# 计算所有指标的统计量
calculate_stats <- function(columns, type_name) {
  stats_df <- data.frame(
    Metric = character(),
    Type = character(),
    Mean = numeric(),
    Variance = numeric(),
    StdDev = numeric(),
    AbsMean = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(col in columns) {
    values <- result1[[col]]
    values <- values[!is.na(values)]  # 移除NA值
    
    if(length(values) > 0) {
      # 确定指标类型
      if(grepl("_r2$", col)) {
        metric_type <- "R-squared"
      } else if(grepl("mse|rmse|mae|median_ae", col)) {
        metric_type <- "Error"
      } else if(grepl("pearson|spearman", col)) {
        metric_type <- "Correlation"
      } else if(grepl("mape|smape|cn_smape", col)) {
        metric_type <- "Percentage Error"
      } else {
        metric_type <- "Other"
      }
      
      stats_df <- rbind(stats_df, data.frame(
        Metric = col,
        Type = metric_type,
        Mean = mean(values),
        Variance = var(values),
        StdDev = sd(values),
        AbsMean = mean(abs(values)),
        Source = type_name
      ))
    }
  }
  return(stats_df)
}

# 计算三种指标的统计量
delta_stats <- calculate_stats(delta_columns, "Delta")
relative_stats <- calculate_stats(relative_columns, "Relative")
standardized_stats <- calculate_stats(standardized_columns, "Standardized")

# 合并所有统计量
all_stats <- rbind(delta_stats, relative_stats, standardized_stats)

print("所有指标的统计量：")
print(all_stats)

# 准备绘图数据（长格式）
library(tidyr)
library(ggplot2)
library(dplyr)

# 创建绘图函数
create_plots <- function(columns, source_name, title_suffix) {
  # 创建长格式数据
  long_data <- result1[columns] %>%
    pivot_longer(cols = everything(), 
                 names_to = "Metric", 
                 values_to = "Value") %>%
    mutate(Metric_clean = gsub(paste0("^", tolower(source_name), "_"), "", Metric))
  
  # 添加类型信息
  if(source_name == "Delta") {
    stats_subset <- delta_stats
  } else if(source_name == "Relative") {
    stats_subset <- relative_stats
  } else {
    stats_subset <- standardized_stats
  }
  
  long_data <- long_data %>%
    left_join(stats_subset[, c("Metric", "Type")], by = "Metric")
  
  # 箱线图
  p1 <- ggplot(long_data, aes(x = reorder(Metric_clean, Value), y = Value, fill = Type)) +
    geom_boxplot(alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste(source_name, title_suffix),
         x = "指标", 
         y = source_name,
         fill = "指标类型") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),        # 修改x轴文字大小
          axis.text.y = element_text(size = 15),                              # 修改y轴文字大小
          axis.title = element_text(size = 14),                               # 修改坐标轴标题大小
          plot.title = element_text(size = 16, face = "bold"),                # 修改标题大小
          legend.text = element_text(size = 10),                              # 修改图例文字大小
          legend.title = element_text(size = 12),                             # 修改图例标题大小
          legend.position = "bottom") +
    coord_flip()
  
  # 方差条形图
  stats_source <- stats_subset
  p2 <- ggplot(stats_source, aes(x = reorder(Metric, Variance), y = Variance, fill = Type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    coord_flip() +
    labs(title = paste(source_name, "指标的方差比较"),
         x = "指标", 
         y = "方差",
         fill = "指标类型") +
    theme_minimal()+
    theme(axis.text.x = element_text(size = 12),        # 修改x轴文字大小
          axis.text.y = element_text(size = 15),        # 修改y轴文字大小
          axis.title = element_text(size = 14),         # 修改坐标轴标题大小
          plot.title = element_text(size = 16, face = "bold"),  # 修改标题大小
          legend.text = element_text(size = 10),        # 修改图例文字大小
          legend.title = element_text(size = 12))       # 修改图例标题大小
  
  return(list(boxplot = p1, variance_plot = p2))
}

# 为三种指标创建图形
delta_plots <- create_plots(delta_columns, "Delta", "训练-测试差异分布 (Test - Train)")
relative_plots <- create_plots(relative_columns, "Relative", "相对变化分布 ((Test-Train)/Train)")
standardized_plots <- create_plots(standardized_columns, "Standardized", "标准化变化分布 ((Test-Train)/SD)")

# 显示图形
print(delta_plots$boxplot)
print(delta_plots$variance_plot)

print(relative_plots$boxplot)
print(relative_plots$variance_plot)

print(standardized_plots$boxplot)
print(standardized_plots$variance_plot)

# #-------新抽样函数：五等分和筛选lambda-------
# perform_quintile_split_analysis <- function(final_data, tree, n_iterations = 50, output_file = "quintile_analysis_results.csv") {
#   
#   # 初始化结果列表
#   results <- list()
#   successful_iterations <- 0
#   total_attempts <- 0
#   max_total_attempts <- n_iterations * 20  # 防止无限循环
#   
#   # 创建输出文件并写入表头
#   header <- data.frame(
#     iteration = integer(),
#     random_seed = integer(),
#     quintile_pair = character(),
#     dataset_type = character(),
#     dataset_size = integer(),
#     lambda = numeric(),
#     ols_r2 = numeric(),
#     pic_r2 = numeric(),
#     phylo_r2 = numeric(),
#     caper_r2 = numeric(),
#     cor_r2 = numeric(),
#     var_r2 = numeric(),
#     resid_r2 = numeric(),
#     lik_r2 = numeric(),
#     stringsAsFactors = FALSE
#   )
#   
#   write.csv(header, output_file, row.names = FALSE)
#   cat("创建输出文件:", output_file, "\n")
#   
#   while (successful_iterations < n_iterations && total_attempts < max_total_attempts) {
#     total_attempts <- total_attempts + 1
#     
#     # 设置随机种子
#     current_seed <- as.integer(Sys.time()) + total_attempts
#     set.seed(current_seed)
#     
#     cat("正在进行第", total_attempts, "次尝试，成功次数:", successful_iterations, "/", n_iterations, "...\n")
#     
#     # 使用tryCatch进行错误处理
#     result <- tryCatch({
#       # 计算五等分的大小
#       n_total <- nrow(final_data)
#       quintile_size <- floor(n_total / 5)  # 每份的大小
#       n_used <- quintile_size * 5  # 实际使用的样本数
#       
#       cat("总样本量:", n_total, "，每等分大小:", quintile_size, "，使用样本数:", n_used)
#       if (n_used < n_total) {
#         cat("，舍弃了", n_total - n_used, "个样本\n")
#       } else {
#         cat("\n")
#       }
#       
#       # 随机选择要使用的样本
#       used_indices <- sample(1:n_total, n_used, replace = FALSE)
#       
#       # 创建五个等分数据集
#       quintiles <- list()
#       for (i in 1:5) {
#         start_idx <- (i-1) * quintile_size + 1
#         end_idx <- i * quintile_size
#         quintile_indices <- used_indices[start_idx:end_idx]
#         quintiles[[i]] <- final_data[quintile_indices, ]
#       }
#       
#       # 计算每个等分的lambda值
#       lambda_values <- numeric(5)
#       valid_quintiles <- list()
#       valid_indices <- integer(0)
#       
#       for (i in 1:5) {
#         quintile_data <- quintiles[[i]]
#         
#         # 对等分数据进行修剪和计算lambda
#         prune_result <- prune_tree_and_data(tree, quintile_data)
#         tree_quintile <- prune_result$tree
#         quintile_data_ordered <- prune_result$data
#         
#         # 计算系统发育信号lambda
#         r2_result <- calculate_all_r2(quintile_data_ordered, tree_quintile, prune_result$comp_data)
#         lambda_values[i] <- r2_result$lambda
#         
#         cat("  等分", i, "的lambda =", lambda_values[i], "\n")
#         
#         # 检查lambda是否大于0.95
#         if (lambda_values[i] > 0.95) {
#           valid_quintiles[[length(valid_quintiles) + 1]] <- list(
#             data = quintile_data_ordered,
#             tree = tree_quintile,
#             comp_data = prune_result$comp_data,
#             lambda = lambda_values[i],
#             index = i
#           )
#           valid_indices <- c(valid_indices, i)
#         }
#       }
#       
#       # 检查是否有至少两个等分满足条件
#       if (length(valid_quintiles) >= 2) {
#         cat("找到", length(valid_quintiles), "个lambda>0.95的等分，使用前两个进行分析\n")
#         
#         # 选择前两个有效等分
#         train_quintile <- valid_quintiles[[1]]
#         test_quintile <- valid_quintiles[[2]]
#         
#         # 对train等分进行计算
#         r2_train <- calculate_all_r2(train_quintile$data, train_quintile$tree, train_quintile$comp_data)
#         
#         # 对test等分进行计算
#         r2_test <- calculate_all_r2(test_quintile$data, test_quintile$tree, test_quintile$comp_data)
#         
#         # 构建结果行 - train数据集
#         result_row_train <- data.frame(
#           iteration = successful_iterations + 1,
#           random_seed = current_seed,
#           quintile_pair = paste(train_quintile$index, "vs", test_quintile$index),
#           dataset_type = "train",
#           dataset_size = nrow(train_quintile$data),
#           lambda = train_quintile$lambda,
#           ols_r2 = r2_train$ols_r2,
#           pic_r2 = r2_train$pic_r2,
#           phylo_r2 = r2_train$phylo_r2,
#           caper_r2 = r2_train$caper_r2,
#           cor_r2 = r2_train$cor_r2,
#           var_r2 = r2_train$var_r2,
#           resid_r2 = r2_train$resid_r2,
#           lik_r2 = r2_train$lik_r2,
#           stringsAsFactors = FALSE
#         )
#         
#         # 构建结果行 - test数据集
#         result_row_test <- data.frame(
#           iteration = successful_iterations + 1,
#           random_seed = current_seed,
#           quintile_pair = paste(train_quintile$index, "vs", test_quintile$index),
#           dataset_type = "test",
#           dataset_size = nrow(test_quintile$data),
#           lambda = test_quintile$lambda,
#           ols_r2 = r2_test$ols_r2,
#           pic_r2 = r2_test$pic_r2,
#           phylo_r2 = r2_test$phylo_r2,
#           caper_r2 = r2_test$caper_r2,
#           cor_r2 = r2_test$cor_r2,
#           var_r2 = r2_test$var_r2,
#           resid_r2 = r2_test$resid_r2,
#           lik_r2 = r2_test$lik_r2,
#           stringsAsFactors = FALSE
#         )
#         
#         # 返回成功的结果
#         list(success = TRUE, 
#              result_train = result_row_train, 
#              result_test = result_row_test)
#         
#       } else {
#         cat("只有", length(valid_quintiles), "个等分满足lambda>0.95的条件，需要重新抽样\n")
#         list(success = FALSE, result_train = NULL, result_test = NULL)
#       }
#       
#     }, error = function(e) {
#       # 捕获错误，返回失败
#       cat("第", total_attempts, "次尝试失败，错误信息:", e$message, "\n")
#       list(success = FALSE, result_train = NULL, result_test = NULL)
#     })
#     
#     # 如果成功，保存结果
#     if (result$success) {
#       successful_iterations <- successful_iterations + 1
#       
#       # 保存train数据集结果
#       results[[length(results) + 1]] <- result$result_train
#       # 保存test数据集结果
#       results[[length(results) + 1]] <- result$result_test
#       
#       # 每次成功迭代后立即写入文件
#       write.table(result$result_train, output_file, 
#                   sep = ",", col.names = FALSE, row.names = FALSE, 
#                   append = TRUE)
#       write.table(result$result_test, output_file, 
#                   sep = ",", col.names = FALSE, row.names = FALSE, 
#                   append = TRUE)
#       
#       cat("第", successful_iterations, "次成功分析完成！使用的等分对:", 
#           result$result_train$quintile_pair, 
#           "，数据集大小: train=", result$result_train$dataset_size, 
#           ", test=", result$result_test$dataset_size, "\n")
#     }
#   }
#   
#   # 检查是否达到目标迭代次数
#   if (successful_iterations < n_iterations) {
#     cat("警告: 只完成了", successful_iterations, "次成功迭代，目标为", n_iterations, "次\n")
#     cat("总尝试次数:", total_attempts, "\n")
#   } else {
#     cat("成功完成所有", n_iterations, "次迭代分析！\n")
#   }
#   
#   # 返回所有成功的结果
#   return(do.call(rbind, results))
# }
# 
# # 使用示例
# # quintile_results <- perform_quintile_split_analysis(final_data, tree, n_iterations = 50)
# # 使用示例（与之前相同）
# prune <- prune_tree_and_data(tree, final_data)
# tree_subset <- prune$tree
# final_data_processed <- prune$data
# comp_data <- prune$comp_data
# quintile_results <- perform_quintile_split_analysis(final_data, tree_subset, n_iterations = 50)
# 
# cat("五等分分析完成！\n")

#--------结果分析----------
result <- read.csv("random_split_results_1115.csv")
# result1 <- result[1:50,]
# result1 <- result[51:100,]
# result1 <- result[101:150,]
# result1 <- result[1:100,]
# result1 <- result[101:200,]
result1 <- result[201:300,]
# result1 <- result[301:400,]
# 计算delta列（test - train）
delta_columns <- c()  # 用于存储新列名

# 遍历所有R2指标，计算delta
r2_methods <- c("ols", "pic", "phylo", "caper", "cor", "var", "resid", "lik")

for(method in r2_methods) {
  train_col <- paste0(method, "_r2_train")
  test_col <- paste0(method, "_r2_test")
  delta_col <- paste0("delta_", method, "_r2")
  
  # 确保列存在后再计算
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    result1[[delta_col]] <- result1[[test_col]] - result1[[train_col]]
    delta_columns <- c(delta_columns, delta_col)
  }
}

# 计算delta列的均值和方差
delta_stats <- data.frame(
  Method = gsub("delta_|_r2", "", delta_columns),
  Mean = sapply(result1[delta_columns], mean, na.rm = TRUE),
  Variance = sapply(result1[delta_columns], var, na.rm = TRUE)
)

print("Delta列的均值和方差：")
print(delta_stats)

# 准备绘图数据（长格式）
library(tidyr)
library(ggplot2)

delta_long <- result1[delta_columns] %>%
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "Delta") %>%
  mutate(Method = gsub("delta_|_r2", "", Method))

# 绘制散点图
ggplot(delta_long, aes(x = Method, y = Delta, color = Method)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  labs(title = "各方法R²差异分布 (Test - Train)",
       x = "方法", 
       y = "ΔR² (Test - Train)",
       caption = "黑色菱形表示均值，虚线y=0作为参考") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# 如果想要更直观地显示方差差异，可以添加箱线图版本
ggplot(delta_long, aes(x = Method, y = Delta, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "各方法R²差异分布 (Test - Train) - 箱线图",
       x = "方法", 
       y = "ΔR² (Test - Train)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#--------另一种形式尝试----------
# 准备绘图数据（长格式）
library(tidyr)
library(ggplot2)
library(dplyr)

delta_long <- result1[delta_columns] %>%
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "Delta") %>%
  mutate(Method = gsub("delta_|_r2", "", Method))

# 为每个方法创建秩次散点图
rank_plot_data <- delta_long %>%
  group_by(Method) %>%
  arrange(Delta) %>%  # 按Delta值排序
  mutate(Rank = row_number()) %>%  # 创建秩次（排序位置）
  ungroup()

# 绘制秩次散点图
ggplot(rank_plot_data, aes(x = Rank, y = Delta, color = Method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ Method, scales = "free_x", nrow = 2) +
  labs(title = "各方法R²差异分布 (Test - Train) - 秩次散点图",
       x = "秩次（按数值大小排序）", 
       y = "ΔR² (Test - Train)",
       caption = "点按数值从小到大排列，x轴表示排序位置") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# 如果希望所有方法在同一张图上，可以用这个版本
ggplot(rank_plot_data, aes(x = Rank, y = Delta, color = Method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "各方法R²差异分布 (Test - Train) - 秩次散点图",
       x = "秩次（按数值大小排序）", 
       y = "ΔR² (Test - Train)",
       caption = "点按数值从小到大排列，x轴表示排序位置") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# 另一种风格：添加连接线，更清晰地显示分布形状
ggplot(rank_plot_data, aes(x = Rank, y = Delta, color = Method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(alpha = 0.5) +  # 添加连接线
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ Method, scales = "free_x", nrow = 2) +
  labs(title = "各方法R²差异分布 (Test - Train) - 带连接线的秩次图",
       x = "秩次（按数值大小排序）", 
       y = "ΔR² (Test - Train)") +
  theme_minimal() +
  theme(legend.position = "none")
#----------相关性分析-------------
# 分析训练集内部各R2指标与lambda_train的关系
train_columns <- c("ols_r2_train", "pic_r2_train", "phylo_r2_train", 
                   "caper_r2_train", "cor_r2_train", "var_r2_train", 
                   "resid_r2_train", "lik_r2_train")

train_cor_results <- data.frame(
  Method = character(),
  Correlation = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for(col in train_columns) {
  if(col %in% colnames(result1) && "lambda_train" %in% colnames(result1)) {
    cor_test <- cor.test(result1[[col]], result1$lambda_train, 
                         use = "complete.obs", method = "pearson")
    
    train_cor_results <- rbind(train_cor_results, data.frame(
      Method = gsub("_r2_train", "", col),
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# 分析测试集内部各R2指标与lambda_test的关系
test_columns <- c("ols_r2_test", "pic_r2_test", "phylo_r2_test", 
                  "caper_r2_test", "cor_r2_test", "var_r2_test", 
                  "resid_r2_test", "lik_r2_test")

test_cor_results <- data.frame(
  Method = character(),
  Correlation = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for(col in test_columns) {
  if(col %in% colnames(result1) && "lambda_test" %in% colnames(result1)) {
    cor_test <- cor.test(result1[[col]], result1$lambda_test, 
                         use = "complete.obs", method = "pearson")
    
    test_cor_results <- rbind(test_cor_results, data.frame(
      Method = gsub("_r2_test", "", col),
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# 分析delta_r2与lambda_test的关系（测试性能变化与lambda的关系）
delta_columns <- c("delta_ols_r2", "delta_pic_r2", "delta_phylo_r2", 
                   "delta_caper_r2", "delta_cor_r2", "delta_var_r2", 
                   "delta_resid_r2", "delta_lik_r2")

delta_cor_results <- data.frame(
  Method = character(),
  Correlation = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

for(col in delta_columns) {
  if(col %in% colnames(result1) && "lambda_test" %in% colnames(result1)) {
    cor_test <- cor.test(result1[[col]], result1$lambda_test, 
                         use = "complete.obs", method = "pearson")
    
    delta_cor_results <- rbind(delta_cor_results, data.frame(
      Method = gsub("delta_|_r2", "", col),
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# 输出结果
print("训练集R2与lambda_train的相关性：")
print(train_cor_results)

print("测试集R2与lambda_test的相关性：")
print(test_cor_results)

print("Delta R2与lambda_test的相关性：")
print(delta_cor_results)

# 可视化相关性结果
library(ggplot2)
library(dplyr)
library(tidyr)

# 合并所有结果
all_results <- bind_rows(
  train_cor_results %>% mutate(Set = "Train"),
  test_cor_results %>% mutate(Set = "Test"),
  delta_cor_results %>% mutate(Set = "Delta")
)

# 创建相关性热图
ggplot(all_results, aes(x = Method, y = Set, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f\n(p=%.3f)", Correlation, P_value)), 
            size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  labs(title = "各方法R2指标与lambda的相关性分析",
       x = "方法", y = "数据集") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 分别绘制训练集和测试集的相关性散点图
# 训练集散点图 - 修改后版本
train_plot_data <- result1 %>%
  dplyr::select(all_of(train_columns), lambda_train) %>%
  pivot_longer(cols = -lambda_train, names_to = "Method", values_to = "R2") %>%
  mutate(Method = gsub("_r2_train", "", Method))

# 合并相关性结果到绘图数据
train_plot_data_with_cor <- train_plot_data %>%
  left_join(train_cor_results, by = "Method") %>%
  mutate(annotation = sprintf("r = %.3f\np = %.3f", Correlation, P_value))

# 为每个方法计算合适的标注位置
annotation_positions <- train_plot_data_with_cor %>%
  group_by(Method) %>%
  summarise(
    x_pos = quantile(lambda_train, 0.7, na.rm = TRUE),
    y_pos = quantile(R2, 0.1, na.rm = TRUE),
    Correlation = first(Correlation),
    P_value = first(P_value),
    annotation = first(annotation)
  )

ggplot(train_plot_data_with_cor, aes(x = lambda_train, y = R2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_text(data = annotation_positions, 
            aes(x = x_pos, y = y_pos, label = annotation),
            size = 3, color = "darkred", fontface = "bold") +
  facet_wrap(~ Method, scales = "free_y", ncol = 4) +
  labs(title = "训练集：各方法R2与lambda的关系",
       x = "Lambda (训练集)", y = "R2") +
  theme_bw()

# 测试集散点图 - 修改后版本
test_plot_data <- result1 %>%
  dplyr::select(all_of(test_columns), lambda_test) %>%
  pivot_longer(cols = -lambda_test, names_to = "Method", values_to = "R2") %>%
  mutate(Method = gsub("_r2_test", "", Method))

# 合并相关性结果到绘图数据
test_plot_data_with_cor <- test_plot_data %>%
  left_join(test_cor_results, by = "Method") %>%
  mutate(annotation = sprintf("r = %.3f\np = %.3f", Correlation, P_value))

# 为每个方法计算合适的标注位置
annotation_positions_test <- test_plot_data_with_cor %>%
  group_by(Method) %>%
  summarise(
    x_pos = quantile(lambda_test, 0.7, na.rm = TRUE),
    y_pos = quantile(R2, 0.1, na.rm = TRUE),
    Correlation = first(Correlation),
    P_value = first(P_value),
    annotation = first(annotation)
  )

ggplot(test_plot_data_with_cor, aes(x = lambda_test, y = R2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_text(data = annotation_positions_test, 
            aes(x = x_pos, y = y_pos, label = annotation),
            size = 3, color = "darkblue", fontface = "bold") +
  facet_wrap(~ Method, scales = "free_y", ncol = 4) +
  labs(title = "测试集：各方法R2与lambda的关系",
       x = "Lambda (测试集)", y = "R2") +
  theme_bw()

# Delta R2散点图 - 新增
delta_plot_data <- result1 %>%
  dplyr::select(all_of(delta_columns), lambda_test) %>%
  pivot_longer(cols = -lambda_test, names_to = "Method", values_to = "Delta_R2") %>%
  mutate(Method = gsub("delta_|_r2", "", Method))

# 合并相关性结果到绘图数据
delta_plot_data_with_cor <- delta_plot_data %>%
  left_join(delta_cor_results, by = "Method") %>%
  mutate(annotation = sprintf("r = %.3f\np = %.3f", Correlation, P_value))

# 为每个方法计算合适的标注位置
annotation_positions_delta <- delta_plot_data_with_cor %>%
  group_by(Method) %>%
  summarise(
    x_pos = quantile(lambda_test, 0.7, na.rm = TRUE),
    y_pos = quantile(Delta_R2, 0.1, na.rm = TRUE),
    Correlation = first(Correlation),
    P_value = first(P_value),
    annotation = first(annotation)
  )

ggplot(delta_plot_data_with_cor, aes(x = lambda_test, y = Delta_R2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "green") +
  geom_text(data = annotation_positions_delta, 
            aes(x = x_pos, y = y_pos, label = annotation),
            size = 3, color = "darkgreen", fontface = "bold") +
  facet_wrap(~ Method, scales = "free_y", ncol = 4) +
  labs(title = "Delta R2与lambda_test的关系",
       x = "Lambda (测试集)", y = "Delta R2") +
  theme_bw()

# 显著性筛选：p < 0.05的相关性
significant_correlations <- all_results %>%
  filter(P_value < 0.05)

print("显著的相关性结果 (p < 0.05)：")
print(significant_correlations)
#------相对差异--------
result <- read.csv("random_split_results_1115.csv")
result1 <- result[1:100,]
# result1 <- result[101:200,]
# result1 <- result[201:300,]
# 计算相对差异列 [deltaR2(Test - Train)]/R2_train
relative_delta_columns <- c()  # 用于存储新列名

# 遍历所有R2指标，计算相对差异
r2_methods <- c("ols", "pic", "phylo", "caper", "cor", "var", "resid", "lik")

for(method in r2_methods) {
  train_col <- paste0(method, "_r2_train")
  test_col <- paste0(method, "_r2_test")
  relative_delta_col <- paste0("relative_delta_", method, "_r2")
  
  # 确保列存在后再计算
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    # 计算相对差异：(test - train) / train
    result1[[relative_delta_col]] <- (result1[[test_col]] - result1[[train_col]]) / result1[[train_col]]
    relative_delta_columns <- c(relative_delta_columns, relative_delta_col)
  }
}

# 计算相对差异列的均值和方差
relative_delta_stats <- data.frame(
  Method = gsub("relative_delta_|_r2", "", relative_delta_columns),
  Mean = sapply(result1[relative_delta_columns], mean, na.rm = TRUE),
  Variance = sapply(result1[relative_delta_columns], var, na.rm = TRUE)
)

print("相对差异列 [(Test - Train)/Train] 的均值和方差：")
print(relative_delta_stats)

# 准备绘图数据（长格式）
library(tidyr)
library(ggplot2)

relative_delta_long <- result1[relative_delta_columns] %>%
  pivot_longer(cols = everything(), 
               names_to = "Method", 
               values_to = "Relative_Delta") %>%
  mutate(Method = gsub("relative_delta_|_r2", "", Method))

# 绘制散点图
ggplot(relative_delta_long, aes(x = Method, y = Relative_Delta, color = Method)) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  labs(title = "各方法R²相对差异分布 [(Test - Train)/Train]",
       x = "方法", 
       y = "相对ΔR² [(Test - Train)/Train]",
       caption = "黑色菱形表示均值，虚线y=0作为参考") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# 箱线图版本
ggplot(relative_delta_long, aes(x = Method, y = Relative_Delta, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "各方法R²相对差异分布 [(Test - Train)/Train] - 箱线图",
       x = "方法", 
       y = "相对ΔR² [(Test - Train)/Train]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 可选：如果相对差异值范围很大，可以使用对数尺度
ggplot(relative_delta_long, aes(x = Method, y = Relative_Delta, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_y_continuous(trans = "pseudo_log", 
                     breaks = c(-10, -1, 0, 1, 10, 100)) +
  labs(title = "各方法R²相对差异分布 [(Test - Train)/Train] - 对数尺度",
       x = "方法", 
       y = "相对ΔR² [(Test - Train)/Train] (伪对数尺度)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))