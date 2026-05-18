#运行调试版
library(dplyr)
library(stringr)
library(ape)
library(caper)
# library(reticulate)
# use_condaenv("r-ggtree", required = TRUE)
setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/get_obs_pred")
#每次分析读取特定数据文件
filename <- "Barnum2024_temperature_clean.csv"
#------通用数据读取和分析代码----------
result_file <- "r2_comparison_results.csv"
if (!file.exists(result_file)) {
  result_header <- data.frame(
    filename = character(),
    raw_r2 = numeric(),
    pic_r2 = numeric(),
    phylo_r2 = numeric(),
    pgls_r2_residual = numeric(),
    pgls_r2_phylo = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(result_header, result_file, row.names = FALSE)
}

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

# 检查匹配结果
cat("原始数据集物种数:", nrow(df), "\n")
cat("匹配成功的物种数:", nrow(final_data), "\n")
cat("匹配的物种:", final_data$genbank_accession, "\n")
#---------通用，修剪树并和数据组合--------
tree_subset <- keep.tip(tree, final_data$GCF)
#移除内部节点
tree_subset$node.label <- NULL
# 检查修剪后的树
cat("修剪后的树包含", length(tree_subset$tip.label), "个物种\n")

# 创建系统发育比较对象
comp_data <- comparative.data(phy = tree_subset, 
                              data = final_data, 
                              names.col = "GCF", 
                              vcv = TRUE, 
                              warn.dropped = TRUE)

# 检查是否有物种被丢弃
cat("被丢弃的物种:", comp_data$dropped$tips, "\n")
cat("被丢弃的数据:", comp_data$dropped$unmatched.rows, "\n")
# 检查是否有物种被丢弃
print(comp_data$dropped)

#重要步骤：设定行名，和树统一顺序
rownames(final_data) <- final_data$GCF
final_data <- final_data[tree_subset$tip.label, ]
#-------测试：计算几种不同的pgls R2----------
# caper::pgls
pgls_model_opt <- pgls(reported_temperature_optimum ~ temperature_optimum, data = comp_data, 
                       lambda = "ML")
caper_r2 <- summary(pgls_model_opt)$r.square
lambda <- as.numeric(summary(pgls_model)$param[2])
# 相关系数平方：R2=r^2
observed <- pgls_model_opt$y
predicted <- fitted(pgls_model_opt)
cor_r2 <- as.numeric(cor(observed, predicted)^2)
# 残差平方和和总平方和
rss <- sum((observed - predicted)^2)  # 残差平方和
tss <- sum((observed - mean(observed))^2)  # 总平方和
var_r2 <- 1 - (rss / tss)

cat("完整数据集pgls R2报告:\n")
#模型X对Y的变异解释程度
cat("1. caper系统发育R²:", round(caper_r2, 4), "\n")
#单纯看回归线对模型散点的拟合程度
cat("2. 相关性R²：", round(cor_r2, 4), "\n")
cat("3. 1-残差平方和/总平方和R²:", round(var_r2, 4), "\n")
#---------ols、pic、white、Ives---------
# 计算原始数据OLS回归的R²
raw_model <- lm(reported_temperature_optimum ~ temperature_optimum, data = final_data)
ols_r2 <- summary(raw_model)$r.squared

# 计算系统发育独立对比 (PIC)
final_data <- final_data[tree_subset$tip.label, ]
pic_obs <- pic(final_data$reported_temperature_optimum, tree_subset)
pic_pred <- pic(final_data$temperature_optimum, tree_subset)
# 创建PIC数据框
pic_df <- data.frame(
  pic_obs = pic_obs,
  pic_pred = pic_pred
)
# PIC回归模型（强制通过原点）
pic_model <- lm(pic_obs~pic_pred - 1)  # -1 表示无截距项
summary(pic_model)
# 计算R²值
pic_r2 <- summary(pic_model)$r.squared

#白化R2
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

cat("系统发育R2 (白化方法):", round(phylo_r2, 4), "\n")


#Ives R2
library(phylolm)
library(rr2)

# 使用phylolm拟合PGLS模型（与rr2兼容）
#r2_resid需要超度量树，如果不是，则需要指定权重
weights <- diag(ape::vcv(tree_subset))
# 使用权重重新拟合模型
gls_model <- gls(reported_temperature_optimum ~ temperature_optimum, 
                 data = final_data,
                 correlation = corBrownian(1, form = ~GCF, phy = tree_subset),
                 weights = varFixed(~weights),
                 method = "ML")

# 计算Ives提出的三种R²
resid_r2 <- R2_resid(gls_model)
# 最省心可以直接算
lik_r2 <- R2_lik(gls_model)
#用留一法算的，计算时间太长了
# r2_pred <- R2_pred(gls_model)


cat("完整数据集R2报告:\n")
cat("1. ols R²:", round(ols_r2, 4), "\n")
cat("2. pic R²：", round(pic_r2, 4), "\n")
cat("3. white R²:", round(phylo_r2, 4), "\n")
cat("4. R².resid:", round(resid_r2, 4), "\n")
cat("5. R².lik:", round(lik_r2, 4), "\n")
