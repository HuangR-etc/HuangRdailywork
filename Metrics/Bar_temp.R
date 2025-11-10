#函数封装版
library(dplyr)
library(stringr)
library(ape)
library(caper)
library(nlme)
library(rr2)
setwd("/home/huangr/projects/dailywork/Metrics")
#每次分析读取特定数据文件
filename <- "/home/huangr/r2_data/Barnum2024_temperature_clean.csv"
#------通用数据读取和分析代码----------
result_file <- "r2_comparison_results.csv"
df <- read.csv(filename)
tree <- read.tree("/home/huangr/r2_data/bac120_r226.tree")
# 读取元数据，只保留需要的列
cols_to_keep <- c("accession", "ncbi_organism_name", "gtdb_taxonomy", 
                  "gtdb_representative", "gtdb_genome_representative",
                  "ncbi_refseq_category", "ncbi_assembly_level")
gtdb_meta <- read.delim("/home/huangr/r2_data/bac120_metadata_r226.tsv", 
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

# 计算各种R2的函数
calculate_all_r2 <- function(final_data, tree_subset, comp_data, 
                             response_var = "reported_temperature_optimum", 
                             predictor_var = "temperature_optimum") {
  
  # 1. 原始OLS回归R2
  raw_model <- lm(as.formula(paste(response_var, "~", predictor_var)), data = final_data)
  ols_r2 <- summary(raw_model)$r.squared
  
  # 2. 系统发育独立对比 (PIC) R2
  pic_obs <- pic(final_data[[response_var]], tree_subset)
  pic_pred <- pic(final_data[[predictor_var]], tree_subset)
  pic_df <- data.frame(pic_obs = pic_obs, pic_pred = pic_pred)
  pic_model <- lm(pic_obs ~ pic_pred - 1)  # 无截距项
  pic_r2 <- summary(pic_model)$r.squared
  
  # 3. caper::pgls R2
  pgls_model <- pgls(as.formula(paste(response_var, "~", predictor_var)), 
                     data = comp_data, lambda = "ML")
  caper_r2 <- summary(pgls_model)$r.squared
  lambda <- as.numeric(summary(pgls_model)$param[2])
  # 4.相关系数平方：R2=r^2
  observed <- pgls_model$y
  predicted <- fitted(pgls_model)
  cor_r2 <- as.numeric(cor(observed, predicted)^2)
  # 5.残差平方和和总平方和
  rss <- sum((observed - predicted)^2)  # 残差平方和
  tss <- sum((observed - mean(observed))^2)  # 总平方和
  var_r2 <- 1 - (rss / tss)
  
  # 6. 系统发育白化R2
  calculate_phylo_r2 <- function(y, yhat, tip_labels, tree, model = "lambda", lambda = 1) {
    stopifnot(length(y) == length(yhat))
    names(y) <- names(yhat) <- tip_labels
    
    # 构造协方差矩阵
    if (model == "BM") {
      V <- ape::vcv(tree)
    } else if (model == "lambda") {
      # 使用lambda调整的协方差矩阵
      if (is.null(lambda)) {
        stop("当model='lambda'时，必须提供lambda参数")
      }
      V <- adjust_vcv_by_lambda(tree, lambda)
    } else {
      stop("目前支持BM和lambda模型")
    }
    
    # 确保顺序一致
    V <- V[tip_labels, tip_labels, drop = FALSE]
    
    # Cholesky 白化
    eps <- 1e-8
    W <- solve(V + diag(eps, nrow(V)))
    L <- chol(W)  # W = LᵀL
    
    # 白化变换
    y_star <- as.numeric(L %*% y)
    yhat_star <- as.numeric(L %*% yhat)
    
    # 在白化后空间计算均值
    ybar_star <- mean(y_star)
    
    RSS <- sum((y_star - yhat_star)^2)
    TSS <- sum((y_star - ybar_star)^2)
    
    R2 <- 1 - RSS/TSS
    return(R2)
  }
  # 辅助函数：根据lambda调整协方差矩阵
  adjust_vcv_by_lambda <- function(tree, lambda) {
    # 获取原始BM模型的协方差矩阵
    V_bm <- ape::vcv(tree)
    
    # 对角线元素（种间方差）
    diag_vals <- diag(V_bm)
    
    # 调整非对角线元素（种间协方差）
    V_lambda <- V_bm
    n <- nrow(V_bm)
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        V_lambda[i, j] <- V_lambda[j, i] <- V_bm[i, j] * lambda
      }
    }
    
    # 确保矩阵是数值稳定的
    eps <- 1e-8
    V_lambda <- V_lambda + diag(eps, n)
    
    return(V_lambda)
  }
  phylo_r2 <- calculate_phylo_r2(
    y = final_data$reported_temperature_optimum,
    yhat = final_data$temperature_optimum,
    tip_labels = final_data$GCF,
    tree = tree_subset,
    model = "lambda",
    lambda = lambda
  )
  
  # 7. Ives R2 
  final_data$weights <- diag(ape::vcv(tree_subset))
  # 使用权重重新拟合模型
  gls_model <- gls(as.formula(paste(response_var, "~", predictor_var)), 
                   data = final_data,
                   correlation = corBrownian(1, form = ~GCF, phy = tree_subset),
                   weights = varFixed(~weights),
                   method = "ML")
  
  # 计算Ives提出的R2
  resid_r2 <- R2_resid(gls_model)
  lik_r2 <- R2_lik(gls_model)
  #计算时间实在太长，而且原理上未必有用，先不算了
  # pred_r2 <- R2_pred(gls_model)
  
  # 返回所有R2值
  return(list(
    ols_r2 = ols_r2,
    pic_r2 = pic_r2,
    phylo_r2 = phylo_r2,
    caper_r2 = caper_r2,
    cor_r2 = cor_r2,
    var_r2 = var_r2,
    resid_r2 = resid_r2,
    lik_r2 = lik_r2,
    lambda = lambda
  ))
}
# 测试确认无误
# result <- calculate_all_r2(final_data, tree_subset, comp_data, 
#                            response_var = "reported_temperature_optimum", 
#                            predictor_var = "temperature_optimum")

#--------抽样准备，通用函数---------
#------新增：创建随机划分结果文件头----------
random_split_file <- "random_split_results.csv"
if (!file.exists(random_split_file)) {
  random_split_header <- data.frame(
    filename = character(),
    iteration = integer(),
    random_seed = integer(),
    train_size = integer(),
    test_size = integer(),
    
    ols_r2_train = numeric(),
    pic_r2_train = numeric(),
    phylo_r2_train = numeric(),
    caper_r2_train = numeric(),
    cor_r2_train = numeric(),
    var_r2_train = numeric(),
    resid_r2_train = numeric(),
    lik_r2_train = numeric(),
    lambda_train = numeric(),
    
    ols_r2_test = numeric(),
    pic_r2_test = numeric(),
    phylo_r2_test = numeric(),
    caper_r2_test = numeric(),
    cor_r2_test = numeric(),
    var_r2_test = numeric(),
    resid_r2_test = numeric(),
    lik_r2_test = numeric(),
    lambda_test = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(random_split_header, random_split_file, row.names = FALSE)
}
perform_random_split_analysis <- function(final_data, tree, n_iterations = 50, train_ratio = 0.8) {
  
  results <- list()
  
  for (i in 1:n_iterations) {
    cat("正在进行第", i, "次随机划分...\n")
    
    # 设置随机种子
    current_seed <- as.integer(Sys.time()) + i
    set.seed(current_seed)
    
    # 随机划分数据
    n_total <- nrow(final_data)
    n_train <- round(n_total * train_ratio)
    train_indices <- sample(1:n_total, n_train, replace = FALSE)
    
    train_data <- final_data[train_indices, ]
    test_data <- final_data[-train_indices, ]
    
    # 对训练集进行计算
    prune_train <- prune_tree_and_data(tree, train_data)
    tree_train <- prune_train$tree
    train_data_ordered <- prune_train$data
    comp_data_train <- prune_train$comp_data
    
    r2_train <- calculate_all_r2(train_data_ordered, tree_train, comp_data_train)
    
    # 对测试集进行计算
    prune_test <- prune_tree_and_data(tree, test_data)
    tree_test <- prune_test$tree
    test_data_ordered <- prune_test$data
    comp_data_test <- prune_test$comp_data
    
    r2_test <- calculate_all_r2(test_data_ordered, tree_test, comp_data_test)
    
    # 保存结果
    result_row <- data.frame(
      filename = filename,
      iteration = i,
      random_seed = current_seed,
      train_size = n_train,
      test_size = n_total - n_train,
      
      ols_r2_train = r2_train$ols_r2,
      pic_r2_train = r2_train$pic_r2,
      phylo_r2_train = r2_train$phylo_r2,
      caper_r2_train = r2_train$caper_r2,
      cor_r2_train = r2_train$cor_r2,
      var_r2_train = r2_train$var_r2,
      resid_r2_train = r2_train$resid_r2,
      lik_r2_train = r2_train$lik_r2,
      lambda_train = r2_train$lambda,
      
      ols_r2_test = r2_test$ols_r2,
      pic_r2_test = r2_test$pic_r2,
      phylo_r2_test = r2_test$phylo_r2,
      caper_r2_test = r2_test$caper_r2,
      cor_r2_test = r2_test$cor_r2,
      var_r2_test = r2_test$var_r2,
      resid_r2_test = r2_test$resid_r2,
      lik_r2_test = r2_test$lik_r2,
      lambda_test = r2_test$lambda
    )
    
    results[[i]] <- result_row
    
    # 每次迭代后立即写入文件（追加模式）
    write.table(result_row, random_split_file, 
                sep = ",", col.names = FALSE, row.names = FALSE, 
                append = TRUE)
    
    cat("第", i, "次划分完成，训练集大小:", n_train, "测试集大小:", n_total - n_train, "\n")
  }
  
  # 返回所有结果
  return(do.call(rbind, results))
}

#--------划分1：对整体数据集进行划分-----------
prune <- prune_tree_and_data(tree, final_data)
tree_subset <- prune$tree
final_data <- prune$data
comp_data <- prune$comp_data
# 计算完整数据的R2
# full_r2 <- calculate_all_r2(final_data_ordered, tree_subset, comp_data)
cat("开始50次随机划分分析...\n")
random_results <- perform_random_split_analysis(final_data, tree, n_iterations = 100)

cat("所有分析完成！结果已保存至:", random_split_file, "\n")
#--------划分2：筛选出原训练集，进行划分-----------
library(dplyr)
train_data <- final_data %>%
  filter(partition_temperature == "train")
prune <- prune_tree_and_data(tree, train_data)
tree_subset <- prune$tree
train_data <- prune$data
comp_data <- prune$comp_data
# 计算完整数据的R2
# full_r2 <- calculate_all_r2(final_data_ordered, tree_subset, comp_data)
cat("开始50次随机划分分析...\n")
random_results <- perform_random_split_analysis(train_data, tree, n_iterations = 100)

cat("所有分析完成！结果已保存至:", random_split_file, "\n")
#--------划分3：筛选出原测试集，进行划分-----------
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
#--------结果分析----------
result <- read.csv("random_split_results.csv")
# result1 <- result[1:50,]
# result1 <- result[51:100,]
# result1 <- result[101:150,]
result1 <- result[1:100,]
# result1 <- result[101:200,]
# result1 <- result[201:300,]

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