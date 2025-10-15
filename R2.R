library(dplyr)
library(stringr)
library(ape)

setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/get_obs_pred")
filename <- "ogt_Cimen2020_bac_ran.csv"

df <- read.csv(filename)
tree <- read.tree("D:/pine/paper_redo/GTDB_tree/bac120_r226.tree")
# 读取元数据，只保留需要的列
cols_to_keep <- c("accession", "ncbi_organism_name", "gtdb_taxonomy", 
                  "gtdb_representative", "gtdb_genome_representative",
                  "ncbi_refseq_category", "ncbi_assembly_level")
gtdb_meta <- read.delim("D:/pine/paper_redo/GTDB_tree/bac120_metadata_r226.tsv", 
                        sep = "\t", 
                        header = TRUE)[, cols_to_keep]
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
cat("R²_resid:", round(r2_resid_result$r2, 4), "\n")

#--------画图---------
library(ggplot2)
library(gridExtra)

# 计算原始数据OLS回归的R²
raw_model <- lm(OGT ~ predicted_OGT, data = final_data)
raw_r_squared <- summary(raw_model)$r.squared

# 1. 原始数据图（带OLS回归）
raw_plot <- ggplot(final_data, aes(x = predicted_OGT, y = OGT)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#E41A1C") +  # 红色回归线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # y=x参考线
  labs(x = "Predicted OGT (°C)", 
       y = "Observed OGT (°C)",
       title = "A. Raw data (OLS regression)") +
  annotate("text", x = min(final_data$predicted_OGT), y = max(final_data$OGT),
           label = sprintf("R² = %.2f", raw_r_squared),
           hjust = 0, vjust = 1, size = 5) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# 2. PIC图（带强制通过原点的回归）
# 构建包含两个R²值的文本标签
pic_labels <- sprintf("PIC R² = %.2f\nR²_resid = %.2f", 
                      r_squared, r2_resid_result$r2_2)

pic_plot <- ggplot(pic_df, aes(x = pic_predicted, y = pic_OGT)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x - 1, se = TRUE, color = "#377EB8") +  # 蓝色回归线
  geom_abline(slope = 1, linetype = "dashed", color = "black") +  # y=x参考线
  labs(x = "PIC of Predicted OGT", 
       y = "PIC of Observed OGT",
       title = "B. Phylogenetic Independent Contrasts") +
  annotate("text", x = min(pic_df$pic_predicted), y = max(pic_df$pic_OGT),
           label = pic_labels,
           hjust = 0, vjust = 1, size = 5, lineheight = 0.8) +  # 添加lineheight调整行距
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# 组合两个子图
combined_plot <- grid.arrange(raw_plot, pic_plot, ncol = 2)

