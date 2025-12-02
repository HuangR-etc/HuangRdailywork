#函数封装版
library(dplyr)
library(stringr)
library(ape)
library(caper)
library(nlme)
library(rr2)
setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/get_obs_pred/cimen")
#每次分析读取特定数据文件
filename <- "ogt_Cimen2020_bac_ran.csv"
#------通用数据读取和分析代码----------
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
#-------数据集1：Cimen2020-------
df$search_name <- gsub("_", " ", df$species)

# 执行匹配
matched_df <- df %>%
  left_join(gtdb_meta, by = c("search_name" = "ncbi_organism_name"))

# 为每个物种选择最佳的GCF编号
matched_df <- matched_df %>%
  group_by(species) %>%
  arrange(
    # 优先选择GTDB代表基因组
    desc(gtdb_representative == "t"),
    # 其次选择NCBI参考基因组
    desc(ncbi_refseq_category == "reference genome"),
    # 然后选择组装级别高的
    desc(factor(ncbi_assembly_level, 
                levels = c("Complete Genome", "Chromosome", "Scaffold", "Contig"))),
    # 最后选择第一个匹配的
    row_number()
  ) %>%
  slice(1) %>%
  ungroup()

# 重命名accession列为GCF
names(matched_df)[names(matched_df) == "accession"] <- "GCF"

matched_df$in_tree <- matched_df$GCF %in% tree$tip.label
# 查看不在树中的物种
if(any(!matched_df$in_tree)) {
  missing_in_tree <- matched_df[!matched_df$in_tree, ]
  cat("以下物种不在发育树中:\n")
  print(missing_in_tree[, c("species", "GCF")])
} else {
  cat("所有物种都在发育树中!\n")
}

library(caper)

# 移除不在树中的物种
final_data <- matched_df[matched_df$in_tree, ]
# 将final_data转换为普通数据框（而不是tibble）
final_data <- as.data.frame(final_data)

# 设置行名为GCF编号（这样comparative.data可以正确匹配）
rownames(final_data) <- final_data$GCF
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
                                  response_var = "OGT", 
                                  predictor_var = "predicted_OGT") {
  
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
#                            response_var = "OGT",
#                            predictor_var = "predicted_OGT")

#--------抽样准备，通用函数---------
#------新增：创建随机划分结果文件头----------
# 更新结果文件头
random_split_file <- "random_split_results_bac_ran.csv"
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
random_results <- perform_random_split_analysis(final_data, tree, n_iterations = 50,train_ratio = 0.8)

cat("所有分析完成！结果已保存至:", random_split_file, "\n")
#-------1118版本结果分析--------
result <- read.csv("random_split_results_bac_ran.csv")
result1 <- result[1:50,]
# result1 <- result[101:200,]
# result1 <- result[201:300,]

# 定义各类指标
r2_metrics <- c("ols", "adjusted")
error_metrics <- c("mse", "rmse", "mae",  "pic_rmse", "pic_mae", "whitened_rm
se", "whitened_mae")
corr_metrics <- c("pearson", "spearman", "pic_pearson", "pic_spearman")
percentage_metrics <- c("smape", "cn_smape","mape")


# 计算新的delta指标：delta = (test - train) / [(test + train)/2]
delta_columns <- c()

# 1. R²类指标
for(metric in r2_metrics) {
  train_col <- paste0(metric, "_r2_train")
  test_col <- paste0(metric, "_r2_test")
  delta_col <- paste0("delta_", metric, "_r2")
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    # 计算均值分母，避免除以零
    mean_value <- (result1[[test_col]] + result1[[train_col]]) / 2
    denominator <- ifelse(abs(mean_value) < 1e-10, sign(mean_value) * 1e-10, mean_value)
    result1[[delta_col]] <- (result1[[test_col]] - result1[[train_col]]) / denominator
    delta_columns <- c(delta_columns, delta_col)
  }
}

# 2. 误差类指标
for(metric in error_metrics) {
  train_col <- paste0(metric, "_train")
  test_col <- paste0(metric, "_test")
  delta_col <- paste0("delta_", metric)
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    # 计算均值分母，避免除以零
    mean_value <- (result1[[test_col]] + result1[[train_col]]) / 2
    denominator <- ifelse(abs(mean_value) < 1e-10, sign(mean_value) * 1e-10, mean_value)
    result1[[delta_col]] <- (result1[[test_col]] - result1[[train_col]]) / denominator
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
  
  if(train_col1 %in% colnames(result1) && test_col1 %in% colnames(result1)) {
    # 使用第一种格式
    delta_col <- paste0("delta_", metric, "_corr")
    mean_value <- (result1[[test_col1]] + result1[[train_col1]]) / 2
    denominator <- ifelse(abs(mean_value) < 1e-10, sign(mean_value) * 1e-10, mean_value)
    result1[[delta_col]] <- (result1[[test_col1]] - result1[[train_col1]]) / denominator
    delta_columns <- c(delta_columns, delta_col)
  } else if(train_col2 %in% colnames(result1) && test_col2 %in% colnames(result1)) {
    # 使用第二种格式
    delta_col <- paste0("delta_", metric)
    mean_value <- (result1[[test_col2]] + result1[[train_col2]]) / 2
    denominator <- ifelse(abs(mean_value) < 1e-10, sign(mean_value) * 1e-10, mean_value)
    result1[[delta_col]] <- (result1[[test_col2]] - result1[[train_col2]]) / denominator
    delta_columns <- c(delta_columns, delta_col)
  }
}

# 4. 百分比误差指标
for(metric in percentage_metrics) {
  train_col <- paste0(metric, "_train")
  test_col <- paste0(metric, "_test")
  delta_col <- paste0("delta_", metric)
  
  if(train_col %in% colnames(result1) && test_col %in% colnames(result1)) {
    # 计算均值分母，避免除以零
    mean_value <- (result1[[test_col]] + result1[[train_col]]) / 2
    denominator <- ifelse(abs(mean_value) < 1e-10, sign(mean_value) * 1e-10, mean_value)
    result1[[delta_col]] <- (result1[[test_col]] - result1[[train_col]]) / denominator
    delta_columns <- c(delta_columns, delta_col)
  }
}

# 计算统计量
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

# 计算delta指标的统计量
delta_stats <- calculate_stats(delta_columns, "Delta")

print("Delta指标的统计量：")
print(delta_stats)

# 准备绘图数据
library(tidyr)
library(ggplot2)
library(dplyr)

# 创建长格式数据
long_data <- result1[delta_columns] %>%
  pivot_longer(cols = everything(), 
               names_to = "Metric", 
               values_to = "Value") %>%
  mutate(Metric_clean = gsub("^delta_", "", Metric))

# 添加类型信息
long_data <- long_data %>%
  left_join(delta_stats[, c("Metric", "Type")], by = "Metric")


# 方差条形图
# 创建符合学术期刊标准的方差条形图
p2 <- ggplot(delta_stats, aes(x = reorder(Metric, Variance), y = Variance, fill = Type)) +
  # 绘制条形图，设置透明度为0.8
  geom_bar(stat = "identity", alpha = 0.8, width = 0.7) +
  
  # 使用Set2配色方案，适合学术出版（颜色对比度适中）
  scale_fill_manual(
    values = c(
      "R-squared" = "#C73E1D",     
      "Error" = "#2E86AB",           
      "Correlation" = "#F18F01",   
      "Percentage Error" = "#A23B72"
    )
  ) +
  
  # 坐标轴翻转，使长指标名称更易读
  coord_flip() +
  
  # 设置图形标题和坐标轴标签
  labs(title = "Variance Comparison of Delta Metrics",
       x = "Metric", 
       y = "Variance",
       fill = "Metric Type") +
  
  # 使用经典黑白主题作为基础
  theme_bw() +
  
  # 详细的主题定制
  theme(
    # 图形标题样式：加粗、居中、适当大小
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, 
                              margin = margin(b = 15)),
    
    # 坐标轴标题样式：适当大小、加粗
    axis.title = element_text(size = 22, face = "bold"),
    
    # X轴文本样式（翻转后为底部横轴）
    axis.text.x = element_text(size = 18, color = "black", 
                               margin = margin(t = 5)),
    
    # Y轴文本样式（翻转后为左侧纵轴）
    axis.text.y = element_text(size = 20, color = "black", 
                               margin = margin(r = 5)),
    
    # 坐标轴线样式：设置颜色和粗细
    axis.line = element_line(color = "black", linewidth = 0.5),
    
    # 面板网格线：只保留主要网格线，颜色调浅
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),  # 去除次要网格线
    
    # 面板边框：去除默认边框，使用坐标轴线
    panel.border = element_blank(),
    
    # 图例样式：设置位置、背景等
    legend.position = "top",  # 图例放在顶部
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 17),
    legend.key = element_rect(fill = "white"),  # 图例键背景
    legend.background = element_rect(fill = "transparent"),  # 透明背景
    
    # 整体边距设置
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  # 在底部和左侧添加坐标轴线（替代被移除的面板边框）
  # 注意：由于coord_flip，实际需要添加在翻转后的位置
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, 
           color = "black", linewidth = 0.5) +  # 左侧轴线
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, 
           color = "black", linewidth = 0.5)    # 底部轴线


# 显示图形

print(p2)
# 1. 计算当前数据集的方差排名
delta_stats$Variance_Rank <- rank(delta_stats$Variance, ties.method = "min")

# 2. 按方差从小到大排序（方差越小越好）
delta_stats_sorted <- delta_stats[order(delta_stats$Variance), ]

# 3. 创建排名结果数据框
ranking_results <- data.frame(
  Metric = delta_stats_sorted$Metric,
  Type = delta_stats_sorted$Type,
  Variance = delta_stats_sorted$Variance,
  Rank = 1:nrow(delta_stats_sorted),
  Dataset = "Dataset_1"  # 根据当前数据集修改
)

# 4. 打印排名结果
print("方差排名结果：")
print(ranking_results)

# 5. 保存排名结果到文件（便于后续整合分析）
# 生成文件名（包含数据集信息）
dataset_id <- "bac_ran"  # 根据当前数据集修改
ranking_filename <- paste0("variance_ranking_", dataset_id, ".csv")
write.csv(ranking_results, ranking_filename, row.names = FALSE)
print(paste("排名结果已保存到:", ranking_filename))