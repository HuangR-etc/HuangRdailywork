library(dplyr)
library(stringr)
library(ape)

setwd("/home/huangr/projects/dailywork")
#每次分析读取特定数据文件
filename <- "ogt_Cimen2020_bac_dis.csv"

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
tree <- read.tree("GTDB_tree/bac120_r226.tree")
# 读取元数据，只保留需要的列
cols_to_keep <- c("accession", "ncbi_organism_name", "gtdb_taxonomy", 
                  "gtdb_representative", "gtdb_genome_representative",
                  "ncbi_refseq_category", "ncbi_assembly_level")
gtdb_meta <- read.delim("GTDB_tree/bac120_metadata_r226.tsv", 
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
# 运行PGLS模型
pgls_model <- pgls(OGT ~ predicted_OGT, data = comp_data)
summary(pgls_model)
# 让PGLS自动估计lambda
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
# 方法1: 使用残差计算R2
pgls_predicted <- predict(pgls_model_opt)
pgls_actual <- final_data_ordered$OGT

# 计算PGLS的R2 (残差方法)
pgls_residuals <- pgls_actual - pgls_predicted
pgls_r2_residual <- 1 - sum(pgls_residuals^2) / sum((pgls_actual - mean(pgls_actual))^2)

# 方法2: 使用白化方法计算PGLS的R2
pgls_phylo_r2 <- calculate_phylo_r2(
  y = final_data_ordered$OGT,
  yhat = pgls_predicted,
  tip_labels = final_data_ordered$GCF,
  tree = tree_subset,
  model = "BM"
)

cat("PGLS R2 (残差方法):", round(pgls_r2_residual, 4), "\n")
cat("PGLS R2 (白化方法):", round(pgls_phylo_r2, 4), "\n")

pgls_summary <- summary(pgls_model_opt)

# 提取PGLS模型的R²值
# 注意：caper包的pgls函数可能不直接提供R²，我们需要计算
pgls_r2 <- 1 - (pgls_summary$sigma^2) / var(final_data_ordered$OGT)
pgls_r2
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


#-----------抽取子集试运行-------------
# 函数：根据系统发育距离抽取子集
get_phylogenetic_subsets <- function(tree, data, n_species = 30) {
  # 计算系统发育距离矩阵
  phylo_dist <- cophenetic(tree)
  
  # 子集A：系统发育关系近的物种
  # 找到距离最近的一对物种，然后扩展
  min_pair <- which(phylo_dist == min(phylo_dist[phylo_dist > 0]), arr.ind = TRUE)[1, ]
  start_species <- rownames(phylo_dist)[min_pair[1]]
  
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

# 函数：计算各种R²指标
calculate_all_r2 <- function(data, tree) {
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
  pgls_r2 <- 1 - sum((pgls_actual - pgls_predicted)^2) / 
    sum((pgls_actual - mean(pgls_actual))^2)
  
  return(list(
    ols_r2 = ols_r2,
    pic_r2 = pic_r2,
    phylo_r2 = phylo_r2,
    pgls_r2 = pgls_r2,
    n_species = nrow(data_ordered)
  ))
}

# 执行比较分析
run_phylogenetic_comparison <- function(tree, data, n_species = 30, n_repeats = 10) {
  all_results <- list()
  
  # 完整数据集的基准结果
  full_result <- calculate_all_r2(data, tree)
  full_result$type <- "完整数据集"
  all_results[["full"]] <- full_result
  
  # 多次重复抽样以减少随机性
  for (i in 1:n_repeats) {
    subsets <- get_phylogenetic_subsets(tree, data, n_species)
    
    for (subset_name in c("subset_A", "subset_B")) {
      subset_data <- subsets[[subset_name]]$data
      subset_tree <- subsets[[subset_name]]$tree
      subset_type <- subsets[[subset_name]]$type
      
      # 只有当有足够数据时才计算
      if (nrow(subset_data) >= 10) {
        result <- calculate_all_r2(subset_data, subset_tree)
        result$type <- subset_type
        result$replicate <- i
        
        key <- paste0(subset_name, "_rep", i)
        all_results[[key]] <- result
      }
    }
  }
  
  return(all_results)
}

# 运行分析
comparison_results <- run_phylogenetic_comparison(tree_subset, final_data, 
                                                  n_species = 30, n_repeats = 10)

# 结果汇总和可视化
summarize_results <- function(results) {
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
      n_species = result$n_species,
      stringsAsFactors = FALSE
    )
  }))
  
  # 按子集类型汇总
  summary <- results_df %>%
    filter(dataset != "full") %>%
    group_by(type) %>%
    summarise(
      mean_ols_r2 = mean(ols_r2, na.rm = TRUE),
      sd_ols_r2 = sd(ols_r2, na.rm = TRUE),
      mean_pic_r2 = mean(pic_r2, na.rm = TRUE),
      sd_pic_r2 = sd(pic_r2, na.rm = TRUE),
      mean_phylo_r2 = mean(phylo_r2, na.rm = TRUE),
      sd_phylo_r2 = sd(phylo_r2, na.rm = TRUE),
      mean_pgls_r2 = mean(pgls_r2, na.rm = TRUE),
      sd_pgls_r2 = sd(pgls_r2, na.rm = TRUE),
      n = n()
    )
  
  return(list(detailed = results_df, summary = summary))
}

# 汇总结果
final_summary <- summarize_results(comparison_results)
print(final_summary$summary)

# 可视化比较
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

# 统计检验：比较近缘vs远缘物种的R²差异
# 对每种R²方法进行t检验
r2_methods <- c("ols_r2", "pic_r2", "phylo_r2", "pgls_r2")
for (method in r2_methods) {
  near_data <- final_summary$detailed %>% 
    filter(type == "近缘物种") %>% pull(method)
  far_data <- final_summary$detailed %>% 
    filter(type == "远缘物种") %>% pull(method)
  
  t_test <- t.test(near_data, far_data)
  cat(sprintf("\n%s方法的t检验结果:\n", method))
  cat(sprintf("近缘物种均值: %.3f, 远缘物种均值: %.3f\n", 
              mean(near_data), mean(far_data)))
  cat(sprintf("t = %.3f, p = %.3f\n", t_test$statistic, t_test$p.value))
}

# 可视化两个子集的系统发育关系
plot_phylogenetic_structure_ape <- function(tree, subset_A, subset_B) {
  # 设置颜色
  tip_colors <- ifelse(tree$tip.label %in% subset_A$data$GCF, "red",
                       ifelse(tree$tip.label %in% subset_B$data$GCF, "blue", "gray"))
  
  # 绘制系统发育树
  par(mar = c(2, 0, 2, 0))
  plot(tree, show.tip.label = TRUE, cex = 0.6, tip.color = tip_colors)
  title("物种在系统发育树中的分布")
  
  # 添加图例
  legend("topright", 
         legend = c("近缘子集", "远缘子集", "其他"),
         fill = c("red", "blue", "gray"),
         cex = 0.8, bty = "n")
}

# 绘制
plot_phylogenetic_structure_ape(tree_subset, subsets$subset_A, subsets$subset_B)