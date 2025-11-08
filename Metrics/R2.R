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

#-------数据集3：Sauer2019---------
# 创建搜索名称：将下划线替换为空格（匹配GTDB格式）
df$search_name <- gsub("_", " ", df$species)

# 标准化GTDB元数据中的物种名（小写并移除多余空格）
gtdb_meta$clean_ncbi_name <- tolower(trimws(gsub("\\s+", " ", gtdb_meta$ncbi_organism_name)))

matched_df <- df %>%
  left_join(gtdb_meta, by = c("search_name" = "clean_ncbi_name"))
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

# 检查是否在树中
matched_df$in_tree <- matched_df$GCF %in% tree$tip.label

# 查看匹配结果
cat("成功匹配的物种数:", sum(!is.na(matched_df$GCF)), "/", nrow(df), "\n")
cat("在树中的物种数:", sum(matched_df$in_tree, na.rm = TRUE), "\n")

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
final_data <- final_data %>%
  dplyr::select(
    species,         # 物种名
    OGT,             # 实测OGT值
    predicted_OGT,   # 预测OGT值
    GCF,             # 基因组编号
  )
final_data <- final_data[complete.cases(final_data[, c("OGT", "predicted_OGT")]), ]
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

#-------完成数据和树的组合，开始模型----------
# PGLS:自动估计lambda的模型
pgls_model_opt <- pgls(OGT ~ predicted_OGT, data = comp_data, 
                       lambda = "ML")
summary(pgls_model_opt)

plot(final_data$predicted_OGT,final_data$OGT)
final_data <- final_data[tree_subset$tip.label, ]
# 计算原始数据OLS回归的R²
raw_model <- lm(OGT ~ predicted_OGT, data = final_data)
raw_r_squared <- summary(raw_model)$r.squared
# 计算系统发育独立对比 (PIC)
pic_OGT <- pic(final_data$OGT, tree_subset)
pic_predicted <- pic(final_data$predicted_OGT, tree_subset)

# 创建PIC数据框
pic_df <- data.frame(
  pic_OGT = pic_OGT,
  pic_predicted = pic_predicted
)

# 添加PIC回归模型（强制通过原点）
pic_model <- lm(pic_OGT ~ pic_predicted - 1)  # -1 表示无截距项
summary(pic_model)

# 计算R²值
r_squared <- summary(pic_model)$r.squared
cat("PIC模型的R-squared:", round(r_squared, 4), "\n")

#----------手动计算R2尝试------------

# 1. 提取实际值和预测值
actual_y <- pic_df$pic_OGT
predicted_y <- predict(pic_model)  # 模型预测值

# 2. 计算残差平方和 (SS_residual)
residuals <- actual_y - predicted_y
SS_residual <- sum(residuals^2)

# 3. 计算总平方和 (SS_total)
# 注意：对于有截距的模型，用 mean(actual_y)；对于无截距模型，用0
mean_y <- 0  # 因为PIC模型是无截距的(y ~ x - 1)
SS_total <- sum((actual_y - mean_y)^2)

# 4. 计算R²
R_squared_manual <- 1 - (SS_residual / SS_total)
# r_squared <- summary(pic_model)$r.squared
#二者确实相等

library(nlme)
calculate_r2_resid <- function(final_data, tree_subset) {
  # 创建系统发育距离矩阵
  phylo_dist <- cophenetic(tree_subset)
  
  # 构建数据框
  comp_data <- comparative.data(phy = tree_subset, 
                                data = final_data, 
                                names.col = "GCF", 
                                vcv = TRUE, 
                                warn.dropped = TRUE)
  
  # 使用gls函数拟合PGLS模型（考虑系统发育结构）
  pgls_model <- gls(OGT ~ predicted_OGT, 
                    data = final_data,
                    correlation = corBrownian(1, form = ~GCF, phy = tree_subset),
                    method = "ML")
  pgls_model_opt <- pgls(OGT ~ predicted_OGT, data = comp_data, 
                         lambda = "ML")
  actual <- final_data$OGT
  # 计算R²_resid
  ss_residual <- sum(residuals(pgls_model)^2)
  ss_total <- sum((actual - mean(actual))^2)
  r2_resid_1 <- 1 - (ss_residual / ss_total)
  ss_residual <- sum(residuals(pgls_model_opt)^2)
  r2_resid_2 <- 1 - (ss_residual / ss_total)
  return(list(r2_1 = r2_resid_1,r2_2 = r2_resid_2, model1 = pgls_model,model2 = pgls_model_opt))
}

# 应用计算
r2_resid_result <- calculate_r2_resid(final_data, tree_subset)

#--------在画图部分之前加入系统发育R2的计算----------
#增加：计算系统发育信号lambda值
calculate_phylogenetic_signal <- function(data, tree, tip_labels_col = "GCF") {
  # 确保数据顺序与树梢一致
  data_ordered <- data[tree$tip.label, ]
  
  # 提取观测值和预测值向量
  obs_ogt <- data_ordered$OGT
  pred_ogt <- data_ordered$predicted_OGT
  names(obs_ogt) <- names(pred_ogt) <- data_ordered[[tip_labels_col]]
  
  # 计算观测OGT的Pagel's lambda
  obs_lambda_result <- phylosig(tree, obs_ogt, method = "lambda", test = FALSE)
  obs_lambda <- ifelse(is.numeric(obs_lambda_result$lambda), 
                       obs_lambda_result$lambda, NA)
  
  # 计算预测OGT的Pagel's lambda
  pred_lambda_result <- phylosig(tree, pred_ogt, method = "lambda", test = FALSE)
  pred_lambda <- ifelse(is.numeric(pred_lambda_result$lambda), 
                        pred_lambda_result$lambda, NA)
  
  return(list(
    obs_lambda = obs_lambda,
    pred_lambda = pred_lambda
  ))
}
# 使用白化方法计算系统发育R2
calculate_phylo_r2 <- function(y, yhat, tip_labels, tree, model = "BM") {
  stopifnot(length(y) == length(yhat))
  names(y) <- names(yhat) <- tip_labels
  
  # 1) 构造协方差 V
  if (model == "BM") {
    V <- ape::vcv(tree) # 布朗运动协方差
  } else {
    stop("目前只实现了BM模型，可根据需要扩展其他模型")
  }
  
  # 确保顺序一致
  V <- V[tip_labels, tip_labels, drop = FALSE]
  
  # 2) 取 W = V^(-1) 并 Cholesky 白化
  # 数值稳定，对角线上加入极小抖动避免奇异
  eps <- 1e-8
  W <- solve(V + diag(eps, nrow(V)))
  L <- chol(W) # W = L^T L
  
  y_star <- as.numeric(L %*% y)
  yhat_star <- as.numeric(L %*% yhat)
  
  # 3) 计算 GLS 加权均值（等价地：用白化后的普通均值即可）
  ybar_star <- mean(y_star)
  
  RSS <- sum((y_star - yhat_star)^2)
  TSS <- sum((y_star - ybar_star)^2)
  
  R2 <- 1 - RSS/TSS
  return(R2)
}

# 应用计算系统发育R2
# 确保数据顺序与树梢标签一致
final_data_ordered <- final_data[tree_subset$tip.label, ]

phylo_r2 <- calculate_phylo_r2(
  y = final_data_ordered$OGT,
  yhat = final_data_ordered$predicted_OGT,
  tip_labels = final_data_ordered$GCF,
  tree = tree_subset,
  model = "BM"
)

cat("系统发育R2 (白化方法):", round(phylo_r2, 4), "\n")

#--------计算PGLS模型的R2----------
pgls_summary <- summary(pgls_model_opt)

# 提取PGLS模型的R²值
# pgls_r2 <- 1 - (pgls_summary$sigma^2) / var(final_data_ordered$OGT)
pgls_r2 <-pgls_summary$r.squared
#--------计算系统发育信号lambda值----------
# 计算观测和预测OGT的lambda值
signal_results <- calculate_phylogenetic_signal(final_data_ordered, tree_subset)

# 从PGLS模型中提取最优lambda值
pgls_opt_lambda <- pgls_model_opt$param["lambda"]

cat("系统发育信号评估:\n")
cat("观测OGT的lambda值:", round(signal_results$obs_lambda, 4), "\n")
cat("预测OGT的lambda值:", round(signal_results$pred_lambda, 4), "\n")
cat("PGLS最优lambda值:", round(pgls_opt_lambda, 4), "\n")

# 进行显著性检验
cat("\n系统发育信号显著性检验:\n")
obs_signal_test <- phylosig(tree_subset, setNames(final_data_ordered$OGT, final_data_ordered$GCF), 
                            method = "lambda", test = TRUE)
pred_signal_test <- phylosig(tree_subset, setNames(final_data_ordered$predicted_OGT, final_data_ordered$GCF), 
                             method = "lambda", test = TRUE)

cat("观测OGT lambda显著性: p =", format.pval(obs_signal_test$P, digits = 3), "\n")
cat("预测OGT lambda显著性: p =", format.pval(pred_signal_test$P, digits = 3), "\n")

#--------收集p值-------
# 1. 原始OLS模型的p值
raw_model_summary <- summary(raw_model)
raw_p_value <- raw_model_summary$coefficients["predicted_OGT", "Pr(>|t|)"]
# 2. PIC模型的p值
pic_model_summary <- summary(pic_model)
pic_p_value <- pic_model_summary$coefficients["pic_predicted", "Pr(>|t|)"]

# 用formatC代替一般的format
format_p_value <- function(p) {
  if (p < 0.001) {
    # 使用 formatC 确保正确显示极小值
    return(formatC(p, format = "e", digits = 2))
  } else {
    return(sprintf("%.3f", p))
  }
}
#--------画图---------
library(ggplot2)
library(gridExtra)

# 1. 原始数据图（带OLS回归）
raw_plot <- ggplot(final_data, aes(x = predicted_OGT, y = OGT)) +
  geom_point(alpha = 0.7, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#E41A9C", size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Predicted OGT (°C)", 
       y = "Observed OGT (°C)",
       title = "A. Raw data (OLS regression)") +
  annotate("text", x = min(final_data$predicted_OGT), y = max(final_data$OGT),
           label = sprintf("p = %s\nR² = %.2f", 
                           format_p_value(raw_p_value),
                           raw_r_squared),
           hjust = 0, vjust = 1, size = 5, lineheight = 0.8) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# 2. PIC图（带强制通过原点的回归）
# 构建包含PIC模型和PGLS模型的文本标签
pic_labels <- sprintf("p = %s\nPIC R² = %.2f\nR²_resid = %.2f\nPhylo R² = %.2f", 
                      format_p_value(pic_p_value), r_squared, 
                      r2_resid_result$r2_2, phylo_r2)

pic_plot <- ggplot(pic_df, aes(x = pic_predicted, y = pic_OGT)) +
  geom_point(alpha = 0.7, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x - 1, se = TRUE, color = "#377EB8", size = 0.5) +
  geom_abline(slope = 1, linetype = "dashed", color = "black") +
  labs(x = "PIC of Predicted OGT", 
       y = "PIC of Observed OGT",
       title = "B. Phylogenetic Independent Contrasts") +
  annotate("text", x = min(pic_df$pic_predicted), y = max(pic_df$pic_OGT),
           label = pic_labels,
           hjust = 0, vjust = 1, size = 5, lineheight = 0.8) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
# 组合两个子图
combined_plot <- grid.arrange(raw_plot, pic_plot, ncol = 2)
# 运行前慎重！确认文件名！别把之前的覆盖了！
# ggsave(filename = "ogt_Cimen2020_bac_dis.png",
#        plot = combined_plot,
#        width = 10,         # 宽度（英寸）
#        height = 6,         # 高度（英寸）
#        dpi = 500,          # 分辨率（每英寸点数）
#        bg = "white")       # 背景颜色
raw_r2 <- summary(raw_model)$r.squared
pic_r2 <- summary(pic_model)$r.squared

new_result <- data.frame(
  filename = basename(filename),
  raw_r2 = raw_r2,
  pic_r2 = pic_r2,
  phylo_r2 = phylo_r2,
  pgls_r2_residual = pgls_r2_residual,
  pgls_r2_phylo = pgls_phylo_r2
)
# 追加新结果
write.table(new_result, result_file, 
            append = TRUE, 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE)

#--------读取最终结果画图----------
all_results <- read.csv(result_file)
r2_plot <- ggplot(all_results, aes(x = raw_r2, y = pic_r2)) +
  geom_point(size = 2, alpha = 0.7, color = "#4DAF4A") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text(aes(label = filename), vjust = -0.8, size = 3, check_overlap = TRUE) +
  labs(x = "Raw R² (OLS regression)",
       y = "PIC R² (Phylogenetic Independent Contrasts)",
       title = "Comparison of R² Values Across Datasets") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_line(color = "gray90"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
r2_plot

#-----------抽取子集比较-------------
# 函数：根据系统发育距离抽取子集
get_phylogenetic_subsets <- function(tree, data, n_species = 30, seed = NULL) {
  # 如果提供了种子，设置随机种子
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # 计算系统发育距离矩阵
  phylo_dist <- cophenetic(tree)
  
  # 子集A：系统发育关系近的物种
  # 随机选取一个物种
  start_species <- sample(rownames(phylo_dist), 1)
  
  # 获取与起始物种最近的n_species个物种
  dist_to_start <- phylo_dist[start_species, ]
  close_species <- names(sort(dist_to_start)[1:n_species])
  
  subset_A <- data[data$GCF %in% close_species, ]
  
  # 子集B：系统发育关系远的物种
  # 方法：从不同主要分支中随机选择
  # 首先将树分成几个主要分支
  n_clusters <- 5  # 分成5个主要分支
  clusters <- cutree(hclust(as.dist(phylo_dist)), k = n_clusters)
  
  # 从每个分支中随机选择大致相等数量的物种
  species_per_cluster <- ceiling(n_species / n_clusters)
  far_species <- c()
  
  for (cluster_id in unique(clusters)) {
    cluster_species <- names(clusters[clusters == cluster_id])
    if (length(cluster_species) > species_per_cluster) {
      selected <- sample(cluster_species, species_per_cluster)
    } else {
      selected <- cluster_species
    }
    far_species <- c(far_species, selected)
  }
  
  # 如果数量超过n_species，随机选择n_species个
  if (length(far_species) > n_species) {
    far_species <- sample(far_species, n_species)
  }
  
  subset_B <- data[data$GCF %in% far_species, ]
  
  # 对应的子树
  tree_A <- keep.tip(tree, subset_A$GCF)
  tree_B <- keep.tip(tree, subset_B$GCF)
  
  return(list(
    subset_A = list(data = subset_A, tree = tree_A, type = "近缘物种"),
    subset_B = list(data = subset_B, tree = tree_B, type = "远缘物种")
  ))
}

# 函数：计算各种R²指标和lambda值
calculate_all_r2_and_lambda <- function(data, tree) {
  results <- list()
  
  # 确保数据顺序与树梢一致
  data_ordered <- data[tree$tip.label, ]
  
  # 1. 普通OLS R²
  ols_model <- lm(OGT ~ predicted_OGT, data = data_ordered)
  ols_r2 <- summary(ols_model)$r.squared
  
  # 2. PIC R²
  pic_OGT <- pic(data_ordered$OGT, tree)
  pic_predicted <- pic(data_ordered$predicted_OGT, tree)
  pic_model <- lm(pic_OGT ~ pic_predicted - 1)
  pic_r2 <- summary(pic_model)$r.squared
  
  # 3. 系统发育R²（白化方法）
  phylo_r2 <- calculate_phylo_r2(
    y = data_ordered$OGT,
    yhat = data_ordered$predicted_OGT,
    tip_labels = data_ordered$GCF,
    tree = tree,
    model = "BM"
  )
  
  # 4. PGLS R²
  comp_data <- comparative.data(phy = tree, data = data_ordered, 
                                names.col = "GCF", vcv = TRUE)
  pgls_model <- pgls(OGT ~ predicted_OGT, data = comp_data, lambda = "ML")
  pgls_predicted <- predict(pgls_model)
  pgls_actual <- data_ordered$OGT
  pgls_r2 <- summary(pgls_model)$r.squared
  
  # 5. 计算系统发育信号lambda值
  lambda_results <- calculate_phylogenetic_signal(data_ordered, tree)
  
  return(list(
    ols_r2 = ols_r2,
    pic_r2 = pic_r2,
    phylo_r2 = phylo_r2,
    pgls_r2 = pgls_r2,
    obs_lambda = lambda_results$obs_lambda,
    pred_lambda = lambda_results$pred_lambda,
    pgls_opt_lambda = pgls_model$param["lambda"],
    n_species = nrow(data_ordered)
  ))
}

# 函数：结果汇总和可视化
summarize_results_with_lambda <- function(results) {
  library(dplyr)
  library(tidyr)
  
  # 转换为数据框
  results_df <- do.call(rbind, lapply(names(results), function(name) {
    result <- results[[name]]
    data.frame(
      dataset = name,
      type = result$type,
      ols_r2 = result$ols_r2,
      pic_r2 = result$pic_r2,
      phylo_r2 = result$phylo_r2,
      pgls_r2 = result$pgls_r2,
      obs_lambda = result$obs_lambda,
      pred_lambda = result$pred_lambda,
      pgls_opt_lambda = result$pgls_opt_lambda,
      n_species = result$n_species,
      seed = result$seed,  # 新增种子列
      replicate = result$replicate,
      stringsAsFactors = FALSE
    )
  }))
  
  # 按子集类型汇总
  summary <- results_df %>%
    filter(dataset != "full") %>%
    group_by(type) %>%
    summarise(
      # R²统计
      mean_ols_r2 = mean(ols_r2, na.rm = TRUE),
      sd_ols_r2 = sd(ols_r2, na.rm = TRUE),
      mean_pic_r2 = mean(pic_r2, na.rm = TRUE),
      sd_pic_r2 = sd(pic_r2, na.rm = TRUE),
      mean_phylo_r2 = mean(phylo_r2, na.rm = TRUE),
      sd_phylo_r2 = sd(phylo_r2, na.rm = TRUE),
      mean_pgls_r2 = mean(pgls_r2, na.rm = TRUE),
      sd_pgls_r2 = sd(pgls_r2, na.rm = TRUE),
      
      # Lambda统计
      mean_obs_lambda = mean(obs_lambda, na.rm = TRUE),
      sd_obs_lambda = sd(obs_lambda, na.rm = TRUE),
      mean_pred_lambda = mean(pred_lambda, na.rm = TRUE),
      sd_pred_lambda = sd(pred_lambda, na.rm = TRUE),
      mean_pgls_lambda = mean(pgls_opt_lambda, na.rm = TRUE),
      sd_pgls_lambda = sd(pgls_opt_lambda, na.rm = TRUE),
      
      n = n()
    )
  
  return(list(detailed = results_df, summary = summary))
}

# 函数：执行比较分析
run_phylogenetic_comparison_with_lambda <- function(tree, data, n_species = 30, n_repeats = 10) {
  
  all_results <- list()
  
  # 生成固定的种子序列，确保可复现性
  seeds <- sample(1:100000, n_repeats)
  
  for (i in 1:n_repeats) {
    cat("正在进行第", i, "次重复抽样，种子:", seeds[i], "...\n")
    
    tryCatch({
      # 使用特定的种子进行抽样
      subsets <- get_phylogenetic_subsets(tree, data, n_species, seed = seeds[i])
      
      for (subset_name in c("subset_A", "subset_B")) {
        subset_data <- subsets[[subset_name]]$data
        subset_tree <- subsets[[subset_name]]$tree
        
        # 检查数据量
        if (nrow(subset_data) < 10) {
          message("子集", subset_name, "第", i, "次重复数据量不足，种子:", seeds[i])
          next
        }
        
        # 使用新的函数计算R²和lambda
        result <- calculate_all_r2_and_lambda(subset_data, subset_tree)
        result$replicate <- i
        result$type <- subsets[[subset_name]]$type
        result$seed <- seeds[i]  # 记录使用的种子
        
        key <- paste0(subset_name, "_rep", i, "_seed", seeds[i])
        all_results[[key]] <- result
      }
    }, error = function(e) {
      message("第", i, "次重复失败，种子:", seeds[i], "错误信息: ", e$message)
      # 即使失败也记录种子信息
      all_results[[paste0("failed_rep", i, "_seed", seeds[i])]] <- list(
        replicate = i,
        seed = seeds[i],
        error = e$message,
        type = "失败"
      )
    })
  }
  
  return(all_results)
}

# 运行分析
comparison_results_with_lambda <- run_phylogenetic_comparison_with_lambda(
  tree_subset, final_data, n_species = 30, n_repeats = 500
)
# 提取结果和抽样数据
# all_results <- comparison_results_with_subsets$results
# all_subsets <- comparison_results_with_subsets$subsets
# 汇总结果
final_summary_with_lambda <- summarize_results_with_lambda(comparison_results_with_lambda)
print(final_summary_with_lambda$summary)
# write.csv(final_summary_with_lambda$detailed,"300_sample+lambda.csv")

#---------统计检验：比较近缘vs远缘物种的R²差异-------------
# 对每种R²方法进行t检验
# 筛选数据：只保留pgls_opt_lambda > 0.95的结果
filtered_data <- final_summary_with_lambda$detailed
original_n <- nrow(filtered_data)
# 筛选条件
filtered_data <- filtered_data %>%
  filter(
    pgls_opt_lambda > 0.9999,
    is.finite(ols_r2),
    is.finite(pic_r2), 
    is.finite(phylo_r2),
    is.finite(pgls_r2),
    is.finite(obs_lambda),
    is.finite(pred_lambda),
    is.finite(pgls_opt_lambda)
  )
cat("\n=== 去除重复子集 ===\n")
before_deduplication <- nrow(filtered_data)

# 对近缘和远缘子集分别去重
filtered_data_dedup <- filtered_data %>%
  group_by(type) %>%9
distinct(pic_r2, .keep_all = TRUE) %>%
  ungroup()

after_deduplication <- nrow(filtered_data_dedup)
removed_duplicates <- before_deduplication - after_deduplication

cat("去重前样本数:", before_deduplication, "\n")
cat("去重后样本数:", after_deduplication, "\n")
cat("移除的重复子集数:", removed_duplicates, "\n")

# 使用去重后的数据
filtered_data <- as.data.frame(filtered_data_dedup)
# 统计筛选后数据集的情况
cat("=== 数据筛选结果 ===\n")
cat("原始数据集样本数:", original_n, "\n")
cat("筛选后数据集样本数:", nrow(filtered_data), "\n\n")

# 统计近缘和远缘子集的个数
subset_counts <- table(filtered_data$type)
cat("=== 筛选后各子集样本数 ===\n")
cat("近缘物种子集数:", subset_counts["近缘物种"], "\n")
cat("远缘物种子集数:", subset_counts["远缘物种"], "\n\n")

# 进行统计检验
cat("=== 统计检验结果 (pgls_opt_lambda = 1) ===\n")

# 定义要检验的R²指标
r2_metrics <- c("ols_r2", "pic_r2", "phylo_r2", "pgls_r2")
results_list <- list()

for (metric in r2_metrics) {
  cat("\n---", metric, "的t检验结果 ---\n")
  
  # 提取近缘和远缘子集的数据
  close_data <- filtered_data[filtered_data$type == "近缘物种", metric]
  far_data <- filtered_data[filtered_data$type == "远缘物种", metric]
  
  # 基本统计信息
  cat("近缘子集 (n =", length(close_data), "): 均值 =", round(mean(close_data, na.rm = TRUE), 4), 
      ", 标准差 =", round(sd(close_data, na.rm = TRUE), 4), "\n")
  cat("远缘子集 (n =", length(far_data), "): 均值 =", round(mean(far_data, na.rm = TRUE), 4),
      ", 标准差 =", round(sd(far_data, na.rm = TRUE), 4), "\n")
  
  # 进行t检验
  t_test <- t.test(close_data, far_data)
  
  # 输出结果
  cat("t值 =", round(t_test$statistic, 4), "\n")
  cat("自由度 =", round(t_test$parameter, 2), "\n")
  cat("p值 =", format.pval(t_test$p.value, digits = 4), "\n")
}

#-------可视化----------
# 可视化比较不同抽样方法的R2均值
library(ggplot2)
library(reshape2)

# 准备绘图数据
plot_data <- final_summary$detailed %>%
  filter(dataset != "full") %>%
  dplyr::select(type, ols_r2, pic_r2, phylo_r2, pgls_r2) %>%
  melt(id.vars = "type", variable.name = "method", value.name = "r2")

# 绘制比较图
ggplot(plot_data, aes(x = type, y = r2, fill = method)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "不同R²计算方法在系统发育近缘vs远缘物种子集中的比较",
       x = "物种子集类型", 
       y = "R²值",
       fill = "计算方法") +
  scale_fill_brewer(palette = "Set1",
                    labels = c("OLS R²", "PIC R²", "系统发育R²", "PGLS R²")) +
  theme_bw() +
  theme(legend.position = "bottom")

# 抽样物种点可视化
visualize_original_sampling_clean <- function(tree, data, n_species = 30, n_clusters = 5) {
  # 使用原始抽样方法进行一次抽样
  subsets <- get_phylogenetic_subsets(tree, data, n_species)
  
  # 计算簇群信息（用于可视化）
  phylo_dist <- cophenetic(tree)
  clusters <- cutree(hclust(as.dist(phylo_dist)), k = n_clusters)
  
  # 创建簇群颜色映射
  cluster_colors <- rainbow(length(unique(clusters)))
  names(cluster_colors) <- unique(clusters)
  
  # 为每个物种分配颜色
  tip_colors <- cluster_colors[as.character(clusters[tree$tip.label])]
  
  # 只标记选中的物种，未选中的物种不显示点
  selected_tips <- tree$tip.label %in% c(subsets$subset_A$data$GCF, subsets$subset_B$data$GCF)
  
  # 设置点的形状和大小：近缘物种用实心圆，远缘物种用三角形，未选中的不显示
  tip_pch <- rep(NA, length(tree$tip.label))  # 默认不显示
  tip_pch[tree$tip.label %in% subsets$subset_A$data$GCF] <- 16  # 实心圆
  tip_pch[tree$tip.label %in% subsets$subset_B$data$GCF] <- 17  # 三角形
  
  tip_cex <- rep(0, length(tree$tip.label))  # 默认大小为0（不显示）
  tip_cex[selected_tips] <- 1.2  # 选中的物种点大小为1.2
  
  # 绘制系统发育树
  par(mar = c(1, 1, 3, 1))
  plot(tree, show.tip.label = FALSE, tip.color = tip_colors, 
       edge.color = "gray70", edge.width = 0.6)
  
  # 只添加选中的物种点（未选中的不会显示）
  tiplabels(pch = tip_pch, cex = tip_cex, col = "black")
  
  title("原始抽样方法示意图\n(实心点:近缘子集, 三角:远缘子集)")
  
  # 添加图例
  legend("topright", 
         legend = c("近缘子集", "远缘子集"),
         pch = c(16, 17),
         cex = 0.8, bty = "n")
  
  # 显示簇群信息
  cluster_info <- table(clusters[subsets$subset_A$data$GCF])
  cat("近缘子集在各簇群的分布:\n")
  print(cluster_info)
  
  cluster_info_far <- table(clusters[subsets$subset_B$data$GCF])
  cat("\n远缘子集在各簇群的分布:\n")
  print(cluster_info_far)
  
  # 返回抽样信息用于后续分析
  return(list(
    subsets = subsets,
    clusters = clusters,
    cluster_colors = cluster_colors
  ))
}
# 可视化随机一次原始抽样
original_sampling_clean <- visualize_original_sampling_clean(tree_subset, final_data, n_species = 30, n_clusters = 5)

#----------转化为超度量树，并绘制热图------------
library(ape)

# 使用chronopl函数进行转换
ultrametric_tree <- chronopl(tree_subset, lambda = 1)

y_obs = final_data_ordered$OGT
names(y_obs) <- final_data_ordered$GCF
y_pred = final_data_ordered$predicted_OGT
names(y_pred) <- final_data_ordered$GCF
residuals <- y_obs - y_pred

# 创建数据框
residual_df <- data.frame(
  species = names(residuals),
  residual = as.numeric(residuals)
)

# 为residual_df增加一列随机抽样数据,表示独立同分布热图
residual_df$random <- rnorm(nrow(residual_df), mean = 0, sd = sd(residual_df$residual))

# 确保残差数据与树中的物种顺序匹配
tip_order <- ultrametric_tree$tip.label
residuals_ordered <- residuals[match(tip_order, names(residuals))]
random_ordered <- residual_df$random[match(tip_order, residual_df$species)]

# 设置颜色梯度
col_palette <- colorRampPalette(c("blue", "white", "red"))(100)
residual_range <- range(c(residuals_ordered, random_ordered), na.rm = TRUE)
residual_colors <- col_palette[cut(residuals_ordered, 100, labels = FALSE)]
random_colors <- col_palette[cut(random_ordered, 100, labels = FALSE)]

# 设置图形布局：树 + 热图1 + 热图2
layout(matrix(1:3, 1, 3), widths = c(10, 1, 1))

# 绘制系统发育树
par(mar = c(2, 0, 2, 0))
plot(ultrametric_tree, show.tip.label = FALSE, cex = 0.8)

# 绘制残差热图 - 修正后的版本
par(mar = c(2, 0, 2, 1))
# 创建一个空图，设置合适的坐标轴范围
plot(0, 0, type = "n", 
     xlim = c(0, 1), 
     ylim = c(0.5, length(tip_order) + 0.5),
     axes = FALSE, 
     xlab = "", 
     ylab = "",
     main = "Residuals")

# 添加热图条 - 确保顺序与树匹配
for(i in 1:length(tip_order)) {
  rect(0, i-0.5, 1, i+0.5, col = residual_colors[i], border = NA)
}

# 绘制随机热图 - 修正后的版本
par(mar = c(2, 0, 2, 1))
# 创建一个空图，设置合适的坐标轴范围
plot(0, 0, type = "n", 
     xlim = c(0, 1), 
     ylim = c(0.5, length(tip_order) + 0.5),
     axes = FALSE, 
     xlab = "", 
     ylab = "",
     main = "Random")

# 添加热图条 - 确保顺序与树匹配
for(i in 1:length(tip_order)) {
  rect(0, i-0.5, 1, i+0.5, col = random_colors[i], border = NA)
}
#----------绘制ols R2和pic R2差距最大时的散点图------------
# 1. 将指定列转换为小数表示（避免科学计数法）
library(dplyr)
detailed_data <- final_summary_with_lambda$detailed
# 筛选条件
detailed_data <- detailed_data %>%
  filter(
    is.finite(pgls_r2),
  )
# 指定需要转换的列
numeric_cols <- c("ols_r2", "pic_r2", "phylo_r2", "pgls_r2", 
                  "obs_lambda", "pred_lambda", "pgls_opt_lambda")

# 将这些列转换为数值型并保留足够的小数位数
detailed_data[numeric_cols] <- lapply(detailed_data[numeric_cols], function(x) {
  as.numeric(format(round(x, 6), scientific = FALSE))
})

# 2. 新增差值列
detailed_data$ols_pic_diff <- detailed_data$ols_r2 - detailed_data$pic_r2

# 3. 找出四种极端情况
# 近缘物种中差值最大（正数绝对值最大）
close_max_diff <- detailed_data %>%
  filter(type == "近缘物种") %>%
  filter(ols_pic_diff == max(ols_pic_diff, na.rm = TRUE)) %>%
  slice(1)  # 如果有多个相同值，取第一个

# 近缘物种中差值最小（负数绝对值最大）
close_min_diff <- detailed_data %>%
  filter(type == "近缘物种") %>%
  filter(ols_pic_diff == min(ols_pic_diff, na.rm = TRUE)) %>%
  slice(1)

# 远缘物种中差值最大（正数绝对值最大）
far_max_diff <- detailed_data %>%
  filter(type == "远缘物种") %>%
  filter(ols_pic_diff == max(ols_pic_diff, na.rm = TRUE)) %>%
  slice(1)

# 远缘物种中差值最小（负数绝对值最大）
far_min_diff <- detailed_data %>%
  filter(type == "远缘物种") %>%
  filter(ols_pic_diff == min(ols_pic_diff, na.rm = TRUE)) %>%
  slice(1)

# 4. 重新抽样获取这四个数据集
extreme_cases <- list(
  close_max = close_max_diff,
  close_min = close_min_diff,
  far_max = far_max_diff,
  far_min = far_min_diff
)

# 函数：重新抽样并获取指定子集
get_specific_subset <- function(tree, data, seed, subset_type) {
  subsets <- get_phylogenetic_subsets(tree, data, n_species = 30, seed = seed)
  
  if (subset_type == "近缘物种") {
    return(subsets$subset_A)
  } else {
    return(subsets$subset_B)
  }
}

# 获取四个极端情况的数据集
extreme_datasets <- list()

for (case_name in names(extreme_cases)) {
  case_data <- extreme_cases[[case_name]]
  seed <- case_data$seed
  subset_type <- case_data$type
  
  cat("正在获取案例:", case_name, "种子:", seed, "类型:", subset_type, "\n")
  
  dataset <- get_specific_subset(tree_subset, final_data, seed, subset_type)
  extreme_datasets[[case_name]] <- list(
    data = dataset$data,
    tree = dataset$tree,
    type = dataset$type,
    seed = seed,
    case_name = case_name,
    ols_pic_diff = case_data$ols_pic_diff
  )
}

# 5. 优化后的绘制散点图函数
plot_extreme_cases_optimized <- function(extreme_datasets, 
                                         point_size = 2,           # 点大小
                                         title_size = 12,          # 标题字体大小
                                         axis_title_size = 11,     # 坐标轴标题字体大小
                                         axis_text_size = 10,      # 坐标轴标签字体大小
                                         case_title_size = 14) {   # 案例标题字体大小
  
  library(ggplot2)
  library(patchwork)
  
  # 存储每个案例的图形
  case_plots <- list()
  
  # 遍历四个极端案例
  for (case_name in names(extreme_datasets)) {
    dataset_info <- extreme_datasets[[case_name]]
    data <- dataset_info$data
    tree <- dataset_info$tree
    type <- dataset_info$type
    seed <- dataset_info$seed
    diff_value <- dataset_info$ols_pic_diff
    
    # 确保数据顺序与树梢一致
    data_ordered <- data[tree$tip.label, ]
    
    # 计算PIC值
    pic_OGT <- pic(data_ordered$OGT, tree)
    pic_predicted <- pic(data_ordered$predicted_OGT, tree)
    
    # 创建数据框用于绘图
    orig_data <- data.frame(
      OGT = data_ordered$OGT,
      predicted_OGT = data_ordered$predicted_OGT
    )
    
    pic_data <- data.frame(
      pic_OGT = pic_OGT,
      pic_predicted = pic_predicted
    )
    
    # 计算线性回归模型用于显示R²
    lm_orig <- lm(OGT ~ predicted_OGT, data = orig_data)
    lm_pic <- lm(pic_OGT ~ pic_predicted, data = pic_data)
    
    r2_orig <- round(summary(lm_orig)$r.squared, 3)
    r2_pic <- round(summary(lm_pic)$r.squared, 3)
    
    # 自定义主题：白色背景，黑色边框，无网格线
    custom_theme <- function(base_size = 10) {
      theme(
        # 设置白色背景
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        # 去除网格线
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # 设置坐标轴
        axis.line = element_line(colour = "black", linewidth = 0.5),
        # 设置文本大小
        axis.title = element_text(size = axis_title_size, colour = "black"),
        axis.text = element_text(size = axis_text_size, colour = "black"),
        plot.title = element_text(size = title_size, hjust = 0.5, face = "bold"),
        # 设置图例（如果有的话）
        legend.position = "none"
      )
    }
    
    # 绘制原始数据散点图
    p_orig <- ggplot(orig_data, aes(x = predicted_OGT, y = OGT)) +
      # 添加散点
      geom_point(color = "blue", alpha = 0.7, size = point_size) +
      # 添加线性拟合线
      geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid", 
                  fill = "lightpink", alpha = 0.3) +
      # 添加R²标注
      annotate("text", x = min(orig_data$predicted_OGT), 
               y = max(orig_data$OGT), 
               label = paste0("R² = ", r2_orig),
               hjust = 0, vjust = 1, size = 4, color = "darkred") +
      # 设置坐标轴标签和标题
      labs(title = paste0(type, " - 原始数据"),
           x = "origin_pred_OGT(°C)", y = "origin_obs_OGT(°C)") +
      # 应用自定义主题
      custom_theme() +
      # 添加全包围边框
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
    
    # 绘制PIC数据散点图
    p_pic <- ggplot(pic_data, aes(x = pic_predicted, y = pic_OGT)) +
      # 添加散点
      geom_point(color = "darkgreen", alpha = 0.7, size = point_size) +
      # 添加线性拟合线
      geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid",
                  fill = "lightpink", alpha = 0.3) +
      # 添加通过原点的对角线（PIC理论线）
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                  color = "gray", linewidth = 0.8) +
      # 添加R²标注
      annotate("text", x = min(pic_data$pic_predicted), 
               y = max(pic_data$pic_OGT), 
               label = paste0("R² = ", r2_pic),
               hjust = 0, vjust = 1, size = 4, color = "darkred") +
      # 设置坐标轴标签和标题
      labs(title = "PIC数据",
           x = "PIC_pred_OGT(°C)", y = "PIC_obs_OGT(°C)") +
      # 应用自定义主题
      custom_theme() +
      # 添加全包围边框
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
    
    # 合并两个子图，并添加案例信息
    combined_plot <- p_orig + p_pic + 
      plot_annotation(
        title = paste0("案例: ", case_name, 
                       " (种子: ", seed, 
                       ", 差值: ", round(diff_value, 4), ")"),
        theme = theme(plot.title = element_text(size = case_title_size, 
                                                hjust = 0.5, 
                                                face = "bold"))
      )
    
    case_plots[[case_name]] <- combined_plot
  }
  
  # 6. 将所有八个子图组合成一个整体图形
  
  # 使用patchwork将四个案例垂直排列
  # 每个案例包含两个水平排列的子图
  final_plot <- (case_plots[["close_max"]] | case_plots[["close_min"]]) / 
    (case_plots[["far_max"]] | case_plots[["far_min"]]) +
    plot_layout(heights = c(1, 1), widths = c(1, 1)) +  # 设置行高和列宽比例
    plot_annotation(
      title = "极端案例系统发育比较分析",
      theme = theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold", margin = margin(b = 10)),
        plot.subtitle = element_text(size = 14, hjust = 0.5, lineheight = 1.2)
      )
    )
  
  return(list(
    individual_plots = case_plots,  # 保留单个案例的图形
    combined_plot = final_plot      # 返回2x2网格组合图形
  ))
  
}

# 7. 使用优化函数绘制图形
# 可根据实际情况调整参数
optimized_plots <- plot_extreme_cases_optimized(
  extreme_datasets,
  point_size = 2.5,        # 点大小
  title_size = 11,         # 子标题字体大小
  axis_title_size = 10,    # 坐标轴标题字体大小
  axis_text_size = 9,      # 坐标轴标签字体大小
  case_title_size = 12     # 案例标题字体大小
)

# 8. 显示组合图形
cat("显示组合图形...\n")
print(optimized_plots$combined_plot)


# 9. 输出极端案例的详细信息
cat("=== 极端案例详细信息 ===\n")
for (case_name in names(extreme_cases)) {
  case <- extreme_cases[[case_name]]
  cat("案例:", case_name, "\n")
  cat("类型:", case$type, "\n")
  cat("种子:", case$seed, "\n")
  cat("重复次数:", case$replicate, "\n")
  cat("OLS R²:", round(case$ols_r2, 4), "\n")
  cat("PIC R²:", round(case$pic_r2, 4), "\n")
  cat("差值 (OLS-PIC):", round(case$ols_pic_diff, 4), "\n")
  cat("物种数量:", case$n_species, "\n")
  cat("------------------------\n")
}

#----转换final_summary数据框-------
library(tidyr)
library(dplyr)

# 方法1：使用pivot_wider
wide_df <- final_summary_with_lambda$detailed %>%
  pivot_wider(
    id_cols = c(replicate, seed, n_species),  # 保持不变的标识列
    names_from = type,                         # 根据type列的值创建新列名
    values_from = c(ols_r2, pic_r2, phylo_r2, pgls_r2, obs_lambda, pred_lambda, pgls_opt_lambda)
  )

# 查看结果
head(wide_df)
write.csv(wide_df,"500sampleResult.csv")
