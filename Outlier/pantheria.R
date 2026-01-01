setwd("D:/Lavender/paper_redo/2025_0604pic/pantheria")
# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ape)
library(geiger)
library(treeplyr)
library(caper)
library(tidyverse)
library(phylolm)
#------计算相关性函数-------
# 函数定义：phylo_correlation_analysis
# 参数：
#   tree_file: 树文件路径
#   trait_data: 数据框，包含物种名和两个性状
#   sp_col: 物种名列名
#   x_col: 第一个性状的列名
#   y_col: 第二个性状的列名
phylo_correlation_analysis <- function(tree, trait_data, sp_col, x_col, y_col) {
  # 2. 准备性状数据
  # 确保物种名列名正确
  trait_data <- as.data.frame(trait_data)
  names(trait_data)[names(trait_data) == sp_col] <- "species"
  
  # 提取需要的数据列并删除NA值
  phylo_data <- trait_data[, c("species", x_col, y_col)]
  phylo_data <- na.omit(phylo_data)
  
  # 3. 检查树与数据的匹配
  tree_names <- tree$tip.label
  
  # 查找不匹配的物种
  missing_in_data <- setdiff(tree_names, phylo_data$species)
  missing_in_tree <- setdiff(phylo_data$species, tree_names)
  
  # 4. 对树和数据进行修剪对齐
  # 修剪树：移除数据中不存在的物种
  if(length(missing_in_data) > 0) {
    tree_trimmed <- drop.tip(tree, missing_in_data)
  } else {
    tree_trimmed <- tree
  }
  
  # 修剪数据：移除树中不存在的物种
  if(length(missing_in_tree) > 0) {
    phylo_data <- phylo_data[phylo_data$species %in% tree_trimmed$tip.label, ]
  }
  
  # 最终检查是否完全匹配
  if(!all(phylo_data$species %in% tree_trimmed$tip.label) || 
     !all(tree_trimmed$tip.label %in% phylo_data$species)) {
    warning("Tree and data not perfectly aligned after trimming")
  }
  
  # 5. 提取两个性状向量
  X1 <- phylo_data[[x_col]]
  X2 <- phylo_data[[y_col]]
  names(X1) <- names(X2) <- phylo_data$species
  
  # 6. 原始数据相关性分析
  pearson_orig <- cor.test(X1, X2, method = "pearson")
  spearman_orig <- cor.test(X1, X2, method = "spearman")
  
  # 7. PIC转换和分析
  pic_X1 <- pic(X1, tree_trimmed)
  pic_X2 <- pic(X2, tree_trimmed)
  
  pearson_pic <- cor.test(pic_X1, pic_X2, method = "pearson")
  spearman_pic <- cor.test(pic_X1, pic_X2, method = "spearman")
  
  # 8. 创建包含所有结果的结构化列表
  results <- list(
    tree = tree_trimmed,
    data = phylo_data,
    n_species = nrow(phylo_data),
    
    original_data = list(
      X1 = X1,
      X2 = X2
    ),
    
    correlations_original = list(
      pearson = pearson_orig,
      spearman = spearman_orig
    ),
    
    pics = list(
      X1_pic = pic_X1,
      X2_pic = pic_X2
    ),
    
    correlations_pic = list(
      pearson = pearson_pic,
      spearman = spearman_pic
    )
  )
  
  # 9. 打印简明结果
  cat("\n===== Phylogenetic Independent Contrasts (PIC) Analysis =====\n")
  cat("Species count:", results$n_species, "\n\n")
  
  cat("Original Data Correlations:\n")
  cat("Pearson correlation (r):", round(pearson_orig$estimate, 4), 
      "p-value:", format.pval(pearson_orig$p.value, digits = 3), "\n")
  cat("Spearman correlation (rho):", round(spearman_orig$estimate, 4), 
      "p-value:", format.pval(spearman_orig$p.value, digits = 3), "\n\n")
  
  cat("PIC Correlations:\n")
  cat("Pearson correlation (r):", round(pearson_pic$estimate, 4), 
      "p-value:", format.pval(pearson_pic$p.value, digits = 3), "\n")
  cat("Spearman correlation (rho):", round(spearman_pic$estimate, 4), 
      "p-value:", format.pval(spearman_pic$p.value, digits = 3), "\n")
  cat("=============================================================\n")
  
  # 10. 返回完整结果
  return(results)
}
match_tree <- function(tree,trait_data,sp_col,x_col,y_col){
  # 确保物种名列名正确
  trait_data <- as.data.frame(trait_data)
  names(trait_data)[names(trait_data) == sp_col] <- "species"
  
  # 提取需要的数据列并删除NA值
  phylo_data <- trait_data[, c("species", x_col, y_col)]
  phylo_data <- na.omit(phylo_data)
  
  # 3. 检查树与数据的匹配
  tree_names <- tree$tip.label
  
  # 查找不匹配的物种
  missing_in_data <- setdiff(tree_names, phylo_data$species)
  missing_in_tree <- setdiff(phylo_data$species, tree_names)
  
  # 4. 对树和数据进行修剪对齐
  # 修剪树：移除数据中不存在的物种
  if(length(missing_in_data) > 0) {
    tree_trimmed <- drop.tip(tree, missing_in_data)
  } else {
    tree_trimmed <- tree
  }
  
  # 修剪数据：移除树中不存在的物种
  if(length(missing_in_tree) > 0) {
    phylo_data <- phylo_data[phylo_data$species %in% tree_trimmed$tip.label, ]
  }
  
  # 最终检查是否完全匹配
  if(!all(phylo_data$species %in% tree_trimmed$tip.label) || 
     !all(tree_trimmed$tip.label %in% phylo_data$species)) {
    warning("Tree and data not perfectly aligned after trimming")
  }
  result <- list(tree_trimmed,phylo_data)
  return(result)
}
pgls_correlation_analysis <- function(tree, trait_data, sp_col, x_col, y_col) {
  trait_data <- as.data.frame(trait_data)
  rownames(trait_data) <- trait_data$species
  names(trait_data)[names(trait_data) == sp_col] <- "species"
  
  phylo_data <- trait_data[, c("species", x_col, y_col)]
  phylo_data <- na.omit(phylo_data)
  
  # 2. 检查树与数据的匹配
  tree_names <- tree$tip.label
  missing_in_data <- setdiff(tree_names, phylo_data$species)
  missing_in_tree <- setdiff(phylo_data$species, tree_names)
  
  # 3. 对树和数据进行修剪对齐
  if(length(missing_in_data) > 0) {
    tree_trimmed <- drop.tip(tree, missing_in_data)
  } else {
    tree_trimmed <- tree
  }
  
  if(length(missing_in_tree) > 0) {
    phylo_data <- phylo_data[phylo_data$species %in% tree_trimmed$tip.label, ]
  }
  
  # 最终检查
  if(!all(phylo_data$species %in% tree_trimmed$tip.label) || 
     !all(tree_trimmed$tip.label %in% phylo_data$species)) {
    warning("树和数据在修剪后仍未完全对齐")
  }
  
  # 5. 定义模型列表
  models <- list(
    OUfixedRoot = "OUfixedRoot",
    OUrandomRoot = "OUrandomRoot",
    lambda = "lambda",
    EB = "EB",
    BM = "BM"
  )
  
  # 存储AIC值、p值和模型名称
  aic_values <- numeric(length(models))
  p_values <- numeric(length(models))
  r_values <- numeric(length(models))
  model_names <- names(models)
  
  # 循环计算每个模型的AIC和p值
  for (i in seq_along(models)) {
    model_formula <- paste(x_col,"~",y_col)
    fit <- phylolm(formula(model_formula), data = phylo_data, phy = tree_trimmed, model = models[[i]])
    aic_values[i] <- summary(fit)$aic
    p_values[i] <- summary(fit)$coefficients[2,4]
    r_values[i] <- summary(fit)$coefficients[2,1]
  }
  
  # 找到AIC最小的模型
  min_aic_index <- which.min(aic_values)
  best_model_p_value <- p_values[min_aic_index]
  best_model_r_value <- r_values[min_aic_index]
  
  # 8. 创建结果对象
  results <- list(
    tree = tree_trimmed,
    data = phylo_data,
    n_species = nrow(phylo_data),
    aic_values = as.character(aic_values[min_aic_index]),
    best_model_name = as.character(models[min_aic_index]),
    p_value = best_model_p_value,
    r_squared = best_model_r_value)
  cat("Best Fitting Model:", results$best_model_name, "\n")
  cat("p-value:", results$p_value, "\n")
  cat("r-value:", results$r_squared, "\n")
  return(results)
}

#-------配对树和数据---------
df_wr05 <- read_tsv("PanTHERIA_1-0_WR05_Aug2008.txt")
tree <- read.nexus("mammal.tre")

# 第一步：准备数据集中的物种名称 (MSW05_Species列)
df_wr05 <- df_wr05 %>%
  mutate(
    # 提取MSW05_Species列的值
    species_from_dataset = `MSW05_Species`,
    # 创建树中使用的格式: 属_种
    species_for_tree = gsub(" ", "_", `MSW05_Binomial`)
  )

# 第二步：提取树节点中的物种名称
extract_tree_species <- function(tip_label) {
  parts <- unlist(strsplit(tip_label, "_"))
  # 取前两个部分（属名+种名）
  paste(parts[1], parts[2], sep = "_")
}

# 创建树tip.label到物种名称的映射
tree_species_map <- tibble(
  original_tip = tree$tip.label,
  species_for_tree = sapply(original_tip, extract_tree_species)
)

# 第三步：数据集中的物种名称匹配
df_wr05 <- df_wr05 %>%
  left_join(tree_species_map, by = "species_for_tree") %>%
  mutate(
    in_tree = !is.na(original_tip),
    matched_tip = ifelse(in_tree, original_tip, NA)
  )

# 第四步：树中的物种名称匹配
tree_species_map <- tree_species_map %>%
  left_join(df_wr05 %>% dplyr::select(species_for_tree, in_dataset = species_from_dataset),
            by = "species_for_tree") %>%
  mutate(in_dataset = ifelse(is.na(in_dataset), FALSE, TRUE))

# 第五步：分析匹配情况
# 1. 数据集中有多少物种在树中
matched_in_data <- sum(df_wr05$in_tree, na.rm = TRUE)
total_in_data <- nrow(df_wr05)
match_percentage_data <- round(matched_in_data / total_in_data * 100, 2)

# 2. 树中有多少物种在数据集中
matched_in_tree <- sum(tree_species_map$in_dataset, na.rm = TRUE)
total_in_tree <- length(tree$tip.label)
match_percentage_tree <- round(matched_in_tree / total_in_tree * 100, 2)

# 打印匹配结果
cat("数据集匹配结果:\n")
cat(sprintf("  总物种数: %d\n", total_in_data))
cat(sprintf("  匹配物种数: %d (%.2f%%)\n", matched_in_data, match_percentage_data))
cat(sprintf("  缺失物种数: %d\n", total_in_data - matched_in_data))

cat("\n系统发育树匹配结果:\n")
cat(sprintf("  总tip数: %d\n", total_in_tree))
cat(sprintf("  匹配tip数: %d (%.2f%%)\n", matched_in_tree, match_percentage_tree))
cat(sprintf("  缺失tip数: %d\n", total_in_tree - matched_in_tree))

# 第六步：获取可以用于分析的完整数据集
df_for_analysis <- df_wr05 %>%
  filter(in_tree) %>%
  dplyr::select(MSW05_Order, MSW05_Family, MSW05_Genus, 
         MSW05_Species,
         AdultBodyMass_g = `5-1_AdultBodyMass_g`,
         PopDensity = `21-1_PopulationDensity_n/km2`,
         matched_tip)

# 第七步：获取匹配的树
matched_tips <- na.omit(df_wr05$matched_tip)
tree_pruned <- keep.tip(tree, matched_tips)

# 第八步：检查修剪后树和数据的匹配
name_check_pruned <- name.check(tree_pruned, df_for_analysis, data.names = matched_tips)
name_check_pruned

# 第十步：准备用于分析的数据
# 确保数据框的行名与树的tip.label匹配
rownames(df_for_analysis) <- df_for_analysis$matched_tip

# 按树tip.label顺序排列数据框
df_for_analysis <- df_for_analysis[tree_pruned$tip.label, ]

# 现在tree_pruned和df_for_analysis可以用于后续的系统发育分析

#-------保存和更新数据-------
write.csv(df_for_analysis,"df_wr05_clean.csv")
write.tree(tree_pruned,"pruned_tree.tre")
#--------绘图：不同目分别-----------
df_clean <- df_wr05 %>%
  dplyr::select(MSW05_Order, 
         AdultBodyMass_g = "5-1_AdultBodyMass_g",
         PopDensity = "21-1_PopulationDensity_n/km2") %>%
  drop_na() %>%
  mutate(log_BodyMass = log10(AdultBodyMass_g + 1),
         log_Density = log10(PopDensity + 1)) %>%
  filter(log_BodyMass > 0, log_Density > 0)  # 过滤无效值

# 统计每个目的物种数量
order_counts <- df_clean %>%
  count(MSW05_Order, name = "SpeciesCount")  # 只保留至少有10个物种的目

df_filtered <- df_clean %>%
  filter(MSW05_Order %in% order_counts$MSW05_Order)

# 1. 整个哺乳动物类群的分析（大尺度）
overall_plot <- ggplot(df_clean, aes(x = log_BodyMass, y = log_Density)) +
  geom_point(alpha = 0.3, color = "grey70", size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  labs(x = "Log10(Body Mass (g))", 
       y = "Log10(Population Density (ind./km²))",
       title = "All Mammals (Large Phylogenetic Scale)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# 2. 计算每个目的回归系数
order_results <- df_filtered %>%
  group_by(MSW05_Order) %>%
  summarise(
    Slope = ifelse(n() > 1, coef(lm(log_Density ~ log_BodyMass))[2], NA),
    Intercept = ifelse(n() > 1, coef(lm(log_Density ~ log_BodyMass))[1], NA),
    R_squared = ifelse(n() > 1, summary(lm(log_Density ~ log_BodyMass))$r.squared, NA),
    SpeciesCount = n()
  ) %>%
  arrange(desc(Slope))  # 按斜率排序

# 3. 绘制所有目的的关系图
orders_plot <- ggplot(df_filtered, aes(x = log_BodyMass, y = log_Density)) +
  geom_point(aes(color = MSW05_Order), alpha = 0.6, size = 1.5) +
  geom_smooth(aes(color = MSW05_Order), method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ MSW05_Order, scales = "free", ncol = 4) +  # 每行4个图
  labs(x = "Log10(Body Mass (g))", 
       y = "Log10(Population Density (ind./km²))",
       title = "Body Mass vs. Population Density by Mammalian Order") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",  # 移除图例
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.grid.minor = element_blank())

# 4. 绘制斜率分布图
slope_plot <- ggplot(order_results, aes(x = reorder(MSW05_Order, Slope), y = Slope)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(x = "Order", y = "Regression Slope",
       title = "Variation in Mass-Density Relationship Across Orders",
       subtitle = "Negative values = density decreases with body size\nPositive values = density increases with body size") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 9))

# 5. 组合图形
combined_plot <- (overall_plot | slope_plot) / orders_plot +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(title = "Phylogenetic Scale Dependence in Mammals",
                  subtitle = "Analysis of Body Mass vs. Population Density Relationships",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                                plot.subtitle = element_text(hjust = 0.5, size = 12)))

# 保存结果
ggsave("full_phylogenetic_scale_analysis.png", combined_plot, 
       width = 16, height = 20, dpi = 300)

# 显示图形
print(combined_plot)

# 输出统计结果
cat("\n回归斜率统计摘要:\n")
print(summary(order_results$Slope))

cat("\n各目回归系数:\n")
print(order_results)


#------------绘图：回归线集合---------------
# 加载包
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)  # 用于避免标签重叠

# 数据准备
df_clean <- df_wr05 %>%
  dplyr::select(MSW05_Order, 
         AdultBodyMass_g = "5-1_AdultBodyMass_g",
         PopDensity = "21-1_PopulationDensity_n/km2") %>%
  filter(!is.na(AdultBodyMass_g), !is.na(PopDensity),
         AdultBodyMass_g > 0, PopDensity > 0) %>%
  mutate(log_BodyMass = log10(AdultBodyMass_g),
         log_Density = log10(PopDensity))

# 计算每个目的回归线
order_data <- df_clean %>%
  group_by(MSW05_Order) %>%
  filter(n() >= 5) %>%  # 只包含至少有5个物种的目
  do({
    model <- lm(log_Density ~ log_BodyMass, data = .)
    min_mass <- min(.$log_BodyMass, na.rm = TRUE)
    max_mass <- max(.$log_BodyMass, na.rm = TRUE)
    
    fit_line <- data.frame(
      log_BodyMass = seq(min_mass, max_mass, length.out = 10),
      Order = unique(.$MSW05_Order)
    )
    fit_line$log_Density <- predict(model, newdata = fit_line)
    
    data.frame(
      Order = unique(.$MSW05_Order),
      slope = coef(model)[2],
      n = nrow(.),
      fit_line
    )
  }) %>%
  ungroup()

# 创建颜色映射
orders <- unique(order_data$Order)
set.seed(123)
color_mapping <- setNames(sample(hue_pal()(length(orders))), orders)

# 创建图例数据
legend_data <- data.frame(
  Order = orders,
  x = 1,
  y = seq_along(orders)
)

# 绘制主图
main_plot <- ggplot() +
  # 所有物种的灰色背景点
  geom_point(data = df_clean, 
             aes(x = log_BodyMass, y = log_Density), 
             color = "gray90", alpha = 0.3, size = 1.5) +
  
  # 每个目的散点（使用对应颜色）
  geom_point(data = df_clean %>% filter(MSW05_Order %in% orders),
             aes(x = log_BodyMass, y = log_Density, color = MSW05_Order),
             alpha = 0.6, size = 1.8) +
  
  # 每个目的回归线
  geom_line(data = order_data, 
            aes(x = log_BodyMass, y = log_Density, color = Order),
            linewidth = 1.2) +
  
  # 整体回归线
  geom_smooth(data = df_clean, 
              aes(x = log_BodyMass, y = log_Density),
              method = "lm", color = "black", linewidth = 1.5,
              linetype = "dashed", se = FALSE) +
  
  # 坐标轴和主题
  labs(x = "Log₁₀(Body Mass (g))", 
       y = "Log₁₀(Population Density (ind./km²))",
       title = "Mammalian Body Mass-Density Relationships by Order") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",  # 在主图中移除图例
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  ) +
  scale_color_manual(values = color_mapping)

# 创建单独的图例
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Order)) +
  geom_point(size = 4) +
  geom_text(aes(label = Order), 
            hjust = 0, nudge_x = 0.1, size = 4, 
            check_overlap = FALSE) +
  scale_color_manual(values = color_mapping) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  expand_limits(x = c(1, 8), y = c(0, length(orders) + 1)) +
  coord_cartesian(clip = "off")

# 组合主图和图例
library(patchwork)

combined_plot <- main_plot + 
  inset_element(legend_plot, 
                left = 0.65, bottom = 0.05, 
                right = 0.99, top = 0.95) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white")))

# 添加整体回归方程
model_all <- lm(log_Density ~ log_BodyMass, data = df_clean)
slope <- format(round(coef(model_all)[2], 2), nsmall = 2)
intercept <- format(round(coef(model_all)[1], 2), nsmall = 2)

combined_plot <- combined_plot +
  annotate("text", 
           x = min(df_clean$log_BodyMass) + 0.5, 
           y = min(df_clean$log_Density) + 0.1,
           label = paste0("Overall regression: y = ", intercept, " + ", slope, "x\n",
                          "(n = ", nrow(df_clean), " species)"),
           color = "black", size = 5, hjust = 0)

print(combined_plot)

#-------数据预处理，把所有目的数据准备好---------
library(purrr)
df_wr05 <- read.csv("df_wr05_clean.csv")[,-1]
tree_pruned <- read.tree("pruned_tree_pantheria.tre")
# 函数：process_order_data
# 目的：处理特定目的数据，转换为标准格式
# 参数：
#   data: 原始数据框 (df_wr05)
#   order_name: 目标目的名称 (e.g., "Carnivora")
process_order_data <- function(data, order_name) {
  # 检查是否有足够物种
  species_count <- data %>%
    filter(MSW05_Order == order_name) %>%
    distinct(matched_tip) %>%
    nrow()
  
  if (species_count < 2) {
    return(NULL)  # 跳过只有一个物种的目
  }
  
  # 处理特定目的数据
  order_data <- data %>%
    filter(MSW05_Order == order_name) %>%
    dplyr::select(
      species = matched_tip,
      Order = MSW05_Order,
      Family = MSW05_Family,
      Genus = MSW05_Genus,
      Species = MSW05_Species,
      BodyMass = AdultBodyMass_g,
      PopDensity = PopDensity
    ) %>%
    filter(
      !is.na(BodyMass),
      !is.na(PopDensity),
      BodyMass > 0,
      PopDensity > 0
    ) %>%
    mutate(
      BodyMass_log = log10(BodyMass),
      PopDensity_log = log10(PopDensity),
      Order = factor(Order)  # 确保目是因子类型
    ) %>%
    distinct(species, .keep_all = TRUE)  # 确保每个物种只有一条记录
  
  return(order_data)
}

# 批量处理所有目
all_orders_data <- df_wr05 %>%
  pull(MSW05_Order) %>%
  unique() %>%
  map_dfr(~ {
    result <- process_order_data(df_wr05, .x)
    if (!is.null(result)) {
      cat("成功处理目:", .x, "- 物种数:", nrow(result), "\n")
      return(result)
    } else {
      cat("跳过目的:", .x, "- 物种数不足\n")
      return(NULL)
    }
  })
#-------打印所有目名--------
all_orders <- unique(all_orders_data$Order)
unique(all_orders_data$Order)
tree_m <- match_tree(tree = tree_pruned,
                     trait_data = all_orders_data,
                     sp_col = "species",
                     x_col = "BodyMass_log",
                     y_col = "PopDensity_log")[[1]]
write.tree(tree_m,"all_mammal_order_tree.tre")
write.csv(all_orders_data,"all_mammal_order.csv")

#-------手动---------
data1 <- all_orders_data %>% filter(Order == "Lagomorpha")
data2 <- all_orders_data %>% filter(Order == "Rodentia")
df <- as.data.frame(rbind(data1,data2))
df <- cbind(df$species,df$Order,df$BodyMass_log,df$PopDensity_log)
df_sub <- df[, c(1, 2, 8, 9)]
write.csv(df_sub,"Lagomorpha_Rodentia_raw_data.csv")

tree <- match_tree(tree = tree_pruned,
                   trait_data = df,
                   sp_col = "species",
                   x_col = "BodyMass_log",
                   y_col = "PopDensity_log")[[1]]
write.tree(tree,"Lagomorpha_Rodentia_pruned_tree.tre")
library(ggplot2)
# 绘制散点图
ggplot(data, aes(x = BodyMass_log, y = PopDensity_log, color = Order)) +
  geom_point(alpha = 0.7, size = 3) +                  # 半透明点，大小=3
  scale_color_manual(                                   # 自定义颜色
    values = c("Cingulata" = "#E41A1C", "Rodentia" = "#377EB8"),
    name = "Order"
  ) +
  labs(
    x = "Log(Body Mass)", 
    y = "Log(Population Density)",
    title = "Body Mass vs Population Density by Order"
  ) +
  theme_minimal() +                                     # 简洁主题
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # 标题居中加粗
    legend.position = "bottom"                          # 图例放底部
  )

# 合并数据（您的原始代码）
data1 <- all_orders_data %>% filter(Order == "Cingulata")
data2 <- all_orders_data %>% filter(Order == "Macroscelidea")
data <- rbind(data1, data2)

# 绘制散点图
ggplot(data, aes(x = BodyMass_log, y = PopDensity_log, color = Order)) +
  geom_point(alpha = 0.7, size = 3) +                  # 半透明点，大小=3
  scale_color_manual(                                   # 自定义颜色
    values = c("Cingulata" = "#E41A1C", "Macroscelidea" = "#377EB8"),
    name = "Order"
  ) +
  labs(
    x = "Log(Body Mass)", 
    y = "Log(Population Density)",
    title = "Body Mass vs Population Density by Order"
  ) +
  theme_minimal() +                                     # 简洁主题
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # 标题居中加粗
    legend.position = "bottom"                          # 图例放底部
  )
# 运行分析
results <- pgls_correlation_analysis(
  tree = tree_pruned,
  trait_data = data,
  sp_col = "species",
  x_col = "BodyMass_log",
  y_col = "PopDensity_log"
)


# 筛选五个目的数据
selected_orders <- c("Scandentia", "Lagomorpha", "Rodentia", "Dermoptera", "Primates")
data <- all_orders_data %>% 
  filter(Order %in% selected_orders)

# 定义五个目对应的颜色（使用ColorBrewer的Set1调色板）
order_colors <- c(
  "Scandentia" = "#E41A1C",    # 红色
  "Lagomorpha" = "#377EB8",    # 蓝色
  "Rodentia" = "#4DAF4A",      # 绿色
  "Dermoptera" = "#984EA3",    # 紫色
  "Primates" = "#FF7F00"       # 橙色
)

# 绘制散点图
ggplot(data, aes(x = BodyMass_log, y = PopDensity_log, color = Order)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = order_colors) +  # 应用自定义颜色
  labs(
    x = "Log(Body Mass)", 
    y = "Log(Population Density)",
    title = "Body Mass vs Population Density by Order"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")  # 图例标题加粗
  )

#--------批量--------
order_results <- list()
error_log <- data.frame(
  Order = character(),
  Reason = character(),
  stringsAsFactors = FALSE
)
# 1. 单独分析每个目（整合pgls结果）
for (order in all_orders) {
  cat("处理目:", order, "\n")
  
  order_data <- all_orders_data %>% 
    filter(Order == order) %>%
    dplyr::select(species, BodyMass_log, PopDensity_log)
  
  # 存储本目所有分析结果
  order_result <- list()
  
  # 运行phylo_correlation_analysis
  order_result$phylo <- tryCatch({
    phylo_correlation_analysis(
      tree = tree_pruned,
      trait_data = order_data,
      sp_col = "species",
      x_col = "BodyMass_log",
      y_col = "PopDensity_log"
    )
  }, error = function(e) {
    error_log <<- rbind(error_log, data.frame(
      Order = order,
      Reason = paste("phylo分析失败:", conditionMessage(e))
    ))
    cat("! phylo分析失败:", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(order_result$phylo)) {
    next
  }
  # 运行pgls_correlation_analysis
  order_result$pgls <- tryCatch({
    pgls_correlation_analysis(
      tree = tree_pruned,
      trait_data = order_data,
      sp_col = "species",
      x_col = "BodyMass_log",
      y_col = "PopDensity_log"
    )
  }, error = function(e) {
    error_log <<- rbind(error_log, data.frame(
      Order = order,
      Reason = paste("pgls分析失败:", conditionMessage(e))
    ))
    cat("! pgls分析失败:", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(order_result$pgls)) {
    next
  }
  
  order_results[[order]] <- order_result
}


# 2. 两两组合分析
ana_orders <- unique(as.data.frame(summary(order_results))$Var1)
# 生成所有可能的目组合
order_pairs <- combn(ana_orders, 2, simplify = FALSE)

# 存储所有组合的结果
all_pairs <- list()

for (pair in order_pairs) {
  order1 <- as.character(pair[1])
  order2 <- as.character(pair[2])
  order1_pgls <- order_results[[order1]]$pgls
  order2_pgls <- order_results[[order2]]$pgls
  cat("\n分析组合:", order1, "+", order2, "\n")
  
  # 获取两个目的数据
  data1 <- all_orders_data %>% 
    filter(Order == order1) %>%
    dplyr::select(species, BodyMass_log, PopDensity_log)
  
  data2 <- all_orders_data %>% 
    filter(Order == order2) %>%
    dplyr::select(species, BodyMass_log, PopDensity_log)
  
  combined_data <- rbind(data1, data2)
  
  
  # 运行融合分析
  combined_result_pic <- tryCatch({
    phylo_correlation_analysis(
      tree = tree_pruned,
      trait_data = combined_data,
      sp_col = "species",
      x_col = "BodyMass_log",
      y_col = "PopDensity_log"
    )
  }, error = function(e) {
    cat("!! 分析组合", order1, "+", order2, "时出错:", e$message, "\n")
    NULL
  })
  combined_result_pgls <- tryCatch({
    pgls_correlation_analysis(
      tree = tree_pruned,
      trait_data = combined_data,
      sp_col = "species",
      x_col = "BodyMass_log",
      y_col = "PopDensity_log"
    )
  }, error = function(e) {
    cat("!! 分析组合", order1, "+", order2, "时出错:", e$message, "\n")
    NULL
  })
  
  # 跳过分析失败的情况
  if (is.null(combined_result_pic)) {
    next
  }
  if (is.null(combined_result_pgls)) {
    next
  }
  
  # 获取单个目的所有相关性结果
  order1_results <- order_results[[order1]]$phylo
  order2_results <- order_results[[order2]]$phylo
  
  # 提取结果信息
  order1_pearson <- order1_results$correlations_original$pearson
  order2_pearson <- order2_results$correlations_original$pearson
  
  safe_extract <- function(obj, path, default = NA) {
    tryCatch({
      reduce(path, ~ .x[[.y]], .init = obj)
    }, error = function(e) default)
  }
  # 创建组合数据框
  pair_summary <- data.frame(
    Order1 = order1,
    Order2 = order2,
    
    Order1_num=order1_results$n_species,
    Order2_num=order2_results$n_species,
    
    # Order1 的各种相关性
    Order1_orig_pearson_cor = order1_pearson$estimate,
    Order1_orig_pearson_p = order1_pearson$p.value,
    Order1_orig_spearman_cor = order1_results$correlations_original$spearman$estimate,
    Order1_orig_spearman_p = order1_results$correlations_original$spearman$p.value,
    Order1_pic_pearson_cor = order1_results$correlations_pic$pearson$estimate,
    Order1_pic_pearson_p = order1_results$correlations_pic$pearson$p.value,
    Order1_pic_spearman_cor = order1_results$correlations_pic$spearman$estimate,
    Order1_pic_spearman_p = order1_results$correlations_pic$spearman$p.value,
    Order1_pgls_cor = ifelse(!is.null(order1_pgls), order1_pgls$r_squared, NA),
    Order1_pgls_p = ifelse(!is.null(order1_pgls), order1_pgls$p_value, NA),
    
    # Order2 的各种相关性
    Order2_orig_pearson_cor = order2_pearson$estimate,
    Order2_orig_pearson_p = order2_pearson$p.value,
    Order2_orig_spearman_cor = order2_results$correlations_original$spearman$estimate,
    Order2_orig_spearman_p = order2_results$correlations_original$spearman$p.value,
    Order2_pic_pearson_cor = order2_results$correlations_pic$pearson$estimate,
    Order2_pic_pearson_p = order2_results$correlations_pic$pearson$p.value,
    Order2_pic_spearman_cor = order2_results$correlations_pic$spearman$estimate,
    Order2_pic_spearman_p = order2_results$correlations_pic$spearman$p.value,
    Order2_pgls_cor = ifelse(!is.null(order2_pgls), order2_pgls$r_squared, NA),
    Order2_pgls_p = ifelse(!is.null(order2_pgls), order2_pgls$p_value, NA),
    
    # 组合后的各种相关性
    Combined_orig_pearson_cor = combined_result_pic$correlations_original$pearson$estimate,
    Combined_orig_pearson_p = combined_result_pic$correlations_original$pearson$p.value,
    Combined_orig_spearman_cor = combined_result_pic$correlations_original$spearman$estimate,
    Combined_orig_spearman_p = combined_result_pic$correlations_original$spearman$p.value,
    Combined_pic_pearson_cor = combined_result_pic$correlations_pic$pearson$estimate,
    Combined_pic_pearson_p = combined_result_pic$correlations_pic$pearson$p.value,
    Combined_pic_spearman_cor = combined_result_pic$correlations_pic$spearman$estimate,
    Combined_pic_spearman_p = combined_result_pic$correlations_pic$spearman$p.value,
    Combined_pgls_cor = ifelse(!is.null(combined_result_pgls), combined_result_pgls$r_squared, NA),
    Combined_pgls_p = ifelse(!is.null(combined_result_pgls), combined_result_pgls$p_value, NA),
    
    stringsAsFactors = FALSE
  )
  
  # 添加到结果列表
  all_pairs[[paste(order1, order2, sep = "+")]] <- pair_summary
}

# 将所有组合的结果合并为一个数据框
all_pairs_df <- bind_rows(all_pairs)

writexl::write_xlsx(all_pairs_df,"all_pairs_result_withNspecies.xlsx")

#-----------分析结果------------
all_pairs_df <- readxl::read_excel("all_pairs_result.xlsx")
# 现在从所有组合中筛选不一致的组合
inconsistent_df <- all_pairs_df %>%
  filter(
    # 情况1：两个单目都不显著，但组合后显著
    (
      (Order1_orig_pearson_p >= 0.05 & Order2_orig_pearson_p >= 0.05) & 
        Combined_orig_pearson_p < 0.05
    ) | 
      # 情况2：两个单目都显著且同向，但组合后变得不显著或方向反转
      (
        (Order1_orig_pearson_p < 0.05 & Order2_orig_pearson_p < 0.05) &
          (sign(Order1_orig_pearson_cor) == sign(Order2_orig_pearson_cor)) &
          (
            Combined_orig_pearson_p >= 0.05 | 
              sign(Combined_orig_pearson_cor) != sign(Order1_orig_pearson_cor)
          )
      )
  ) %>%
  rowwise() %>%
  mutate(
    # 增加一个标记列，说明不一致的原因
    inconsistent_reason = case_when(
      Order1_orig_pearson_p >= 0.05 & Order2_orig_pearson_p >= 0.05 & Combined_orig_pearson_p < 0.05 ~ 
        "两个目的原始相关性都不显著，但合并后显著",
      Order1_orig_pearson_p < 0.05 & Order2_orig_pearson_p < 0.05 & sign(Order1_orig_pearson_cor) == sign(Order2_orig_pearson_cor) & Combined_orig_pearson_p >= 0.05 ~ 
        "两个目的原始相关性都显著且方向一致，但合并后不再显著",
      Order1_orig_pearson_p < 0.05 & Order2_orig_pearson_p < 0.05 & sign(Order1_orig_pearson_cor) == sign(Order2_orig_pearson_cor) & sign(Combined_orig_pearson_cor) != sign(Order1_orig_pearson_cor) ~ 
        "两个目的原始相关性都显著且方向一致，但合并后方向反转"
    )
  ) %>%
  ungroup()

# 输出不一致组合的数量
cat("\n=== 共发现", nrow(inconsistent_df), "个不一致组合 ===\n")


# 筛选两个Order各自至少有一种相关性显著且方向一致的情况
consistent_pairs_df <- all_pairs_df %>%
  mutate(
    # 检查Order1是否有任何显著的正/负相关性
    Order1_pos = (Order1_pic_pearson_p < 0.05 & Order1_pic_pearson_cor > 0),
    Order1_neg =(Order1_pic_pearson_p < 0.05 & Order1_pic_pearson_cor < 0) ,
    
    # 检查Order2是否有任何显著的正/负相关性
    Order2_pos = (Order2_pic_pearson_p < 0.05 & Order2_pic_pearson_cor > 0),
    Order2_neg =(Order2_pic_pearson_p < 0.05 & Order2_pic_pearson_cor < 0)
  ) %>%
  filter(
    # 两个Order都有显著的正相关性，或都有显著的负相关性
    (Order1_pos & Order2_pos) | (Order1_neg & Order2_neg)
  ) %>%
  dplyr::select(-Order1_pos, -Order1_neg, -Order2_pos, -Order2_neg)  # 移除临时列

# 输出筛选结果数量
cat("\n=== 共发现", nrow(consistent_pairs_df), 
    "个组合：两个目的各自至少有一种相关性显著且方向一致 ===\n")

# 打印前几个不一致组合（如果需要）
if (nrow(inconsistent_df) > 0) {
  print(head(inconsistent_df))
}


#------图：pic散点图-------
# 加载必要的包
library(tidyverse)
library(ggpubr)
library(patchwork)

# 获取所有目的列表
all_orders <- unique(all_orders_data$Order)

# 1. 创建存储所有PIC数据的数据框
all_pic_data <- data.frame()

# 2. 遍历每个目，读取PIC数据并计算回归统计量
order_results_df <- data.frame()  # 存储回归结果

for (order in all_orders) {
  # 检查该目是否有PIC数据
  if (!is.null(order_results[[order]]) && 
      !is.null(order_results[[order]]$phylo$pics)) {
    # 读取PIC数据文件
    if (1) {
      pic_x1 <- order_ <- ults[[order]]$phylo$pics$X1_pic
      pic_x2 <- order_results[[order]]$phylo$pics$X2_pic
      
      # 创建当前目的的数据框
      order_df <- data.frame(
        Order = order,
        pic_x1 = pic_x1,
        pic_x2 = pic_x2
      )
      
      # 添加到总数据框
      all_pic_data <- bind_rows(all_pic_data, order_df)
      
      # 计算回归统计量（强制通过原点）
      if (length(pic_x1) >= 3) {  # 至少3个对比点
        # 提取相关系数结果
        cor_results <- order_results[[order]]$phylo$correlations_pic$pearson
        
        # 创建回归结果行
        results_row <- data.frame(
          Order = order,
          Slope = lm(pic_x2 ~ pic_x1 + 0)$coefficients["pic_x1"],
          R_squared = cor_results$estimate^2,  # R² = 相关系数的平方
          P_value = cor_results$p.value,
          ContrastsCount = length(pic_x1)  # 对比点数量
        )
        order_results_df <- bind_rows(order_results_df, results_row)
      }
    }
  }
}

# 3. 准备排序后的回归结果数据
order_results_df <- order_results_df %>%
  drop_na() %>%
  arrange(desc(Slope)) %>%
  mutate(Slope = round(Slope, 4),
         R_squared = round(R_squared, 4))

# 4. 绘制分面散点图（每个目）
# 确保使用相同的颜色方案
palette <- scales::hue_pal()(length(all_orders))

pic_facets_plot <- ggplot(all_pic_data, aes(x = pic_x1, y = pic_x2)) +
  geom_point(aes(color = Order), alpha = 0.7, size = 2) +
  geom_smooth(aes(color = Order), method = "lm", formula = y ~ x + 0, 
              se = FALSE, linewidth = 1) +
  facet_wrap(~ Order, scales = "free", ncol = 4) +
  scale_color_manual(values = setNames(palette, all_orders)) +
  labs(x = "PIC of Body Mass", 
       y = "PIC of Population Density",
       title = "Phylogenetic Independent Contrasts by Order") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"))

# 5. 绘制斜率条形图（修正版）
# 先创建一个有序的因子变量
order_results_df <- order_results_df %>%
  mutate(Order_fct = factor(Order, levels = Order[order(Slope)]))

slope_plot <- ggplot(order_results_df, aes(x = Order_fct, y = Slope)) +
  geom_bar(stat = "identity", aes(fill = Slope > 0), alpha = 0.8) +
  geom_text(aes(label = paste0("Slope=", round(Slope, 3), 
                               "\nR²=", round(R_squared, 3), 
                               "\np=", round(P_value, 3))), 
            hjust = -0.05, size = 3, angle = 90) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "firebrick")) +
  coord_flip(ylim = c(min(order_results_df$Slope) * 1.5, 
                      max(order_results_df$Slope) * 1.2)) +
  labs(x = NULL, y = "Regression Slope (Origin Constrained)",
       title = "Slope Variation in Mass-Density PIC Relationships") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(hjust = 0.5))
# 6. 组合图形
combined_plot <- pic_facets_plot / slope_plot +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "Phylogenetically Independent Contrasts Analysis",
    subtitle = "Body Mass vs. Population Density Relationships Across Mammalian Orders",
    caption = paste("Total orders:", nrow(order_results_df), 
                    "| Total contrasts:", nrow(all_pic_data)),
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5))
  )

print(combined_plot)

# 8. 输出统计结果
cat("\n=== 各目PIC回归结果 ===\n")
print(order_results_df)

cat("\n斜率统计摘要:\n")
summary(order_results_df$Slope)

#------手动成对pic散点图(分别）------
# 选择目
selected_orders <- c("Cingulata","Didelphimorphia")

# 筛选数据
pic_data_selected <- all_pic_data %>% 
  filter(Order %in% selected_orders)

# 定义颜色方案
order_colors <- c(
  # "Scandentia" = "#E41A1C",    # 红色
  "Cingulata" = "#377EB8",    # 蓝色
  "Didelphimorphia" = "#E41A1C"     # 绿色
  # "Dermoptera" = "#984EA3",    # 紫色
  # "Primates" = "#FF7F00"       # 橙色
)

# 绘制聚合散点图
aggregate_plot <- ggplot(pic_data_selected, aes(x = pic_x1, y = pic_x2, color = Order)) +
  geom_point(alpha = 0.7, size = 3, shape = 16) +
  #平滑回归线
  # geom_smooth(method = "lm", formula = y ~ x + 0, se = TRUE,
  #             size = 1.2, alpha = 0.2, aes(fill = Order)) +
  scale_color_manual(values = order_colors) +
  scale_fill_manual(values = order_colors) +
  labs(
    x = "PIC of Body Mass", 
    y = "PIC of Population Density",
    title = "Phylogenetically Independent Contrasts (PIC) Analysis",
    subtitle = "Body Mass vs Population Density Relationships"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    axis.title = element_text(face = "bold", size = 13),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )
# +
#   # 添加回归方程和统计值到图中
#   stat_regline_equation(
#     formula = y ~ x + 0,
#     aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#     size = 4, color = "black",
#     label.x.npc = "left",
#     label.y.npc = "top",
#     vjust = 0.5, hjust = -0.1
#   )

# 显示图形
print(aggregate_plot)


#------手动成对pic，先合并数据再算pic---------
# 美化散点图
library(ggplot2)

# 筛选数据并计算PIC
combined_data <- all_orders_data %>%
  filter(Order %in% c("Lagomorpha","Rodentia")) %>%
  mutate(Order = as.factor(Order))
# results1 <- phylo_correlation_analysis(
#   tree = tree,
#   trait_data = combined_data,
#   sp_col = "species",
#   x_col = "BodyMass_log",
#   y_col = "PopDensity_log"
# )
# results2 <- pgls_correlation_analysis(
#   tree = tree,
#   trait_data = combined_data,
#   sp_col = "species",
#   x_col = "BodyMass_log",
#   y_col = "PopDensity_log"
# )
# 
combined_data <- read.csv("Lagomorpha_Rodentia_raw_data.csv")[,-1]
tree <- read.tree("Lagomorpha_Rodentia_pruned_tree.tre")

plot(tree, label.offset = 0.1)  # 调整标签间距
nodelabels(cex = 0.8)          # 显示节点编号（可选）
tiplabels(tree$tip.label,       # 添加物种名标签
          adj = c(-0.1, 0.5),  # 右对齐（向左偏移）
          cex = 0.8,           # 字体大小
          font = 3)            # 斜体字体
# 计算PIC
combined_pic <- data.frame(
  BodyMass_pic = ape::pic(combined_data$BodyMass_log, tree),
  Density_pic = ape::pic(combined_data$PopDensity_log, tree)
)
# 进行线性回归分析
pearson_orig <- cor.test(combined_pic$Density_pic,combined_pic$BodyMass_pic, method = "pearson")
pearson_orig$estimate
pearson_orig$p.value

data1 <- all_orders_data %>% filter(Order == "Macroscelidea")
data2 <- all_orders_data %>% filter(Order == "Cetacea")
data <- rbind(data1,data2)

plot(combined_pic$BodyMass_pic,combined_pic$Density_pic)
# 绘制散点图并添加统计信息
ggplot(combined_pic, aes(x = BodyMass_pic, y = Density_pic)) +
  geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Phylogenetically Independent Contrasts Analysis",
    subtitle = "Lagomorpha_Rodentia_Primates",
    x = "pic log(Body Mass)", 
    y = "pic log(Population Density)"
  ) +
  theme_bw() +  # 白底主题
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  annotate("text", x = min(combined_pic$BodyMass_pic) * 0.9, 
           y = max(combined_pic$Density_pic) * 0.9,
           label = paste("pics number", nrow(combined_pic)), 
           size = 3.5, color = "black", hjust = 0)

#------另一篇论文---------
#-------配对树和数据---------
trophic <- readxl::read_xlsx("TrophicLevel_Appendix1.xlsx")
tree <- read.nexus("mammal.tre")

# 加载必要的包
library(ape)
library(dplyr)
library(stringr)

# 步骤1: 从树中提取规范化的物种名称
extract_tree_species <- function(tip_label) {
  # 分割字符串并提取前两部分（属名+种名）
  parts <- unlist(strsplit(tip_label, "_"))
  paste(parts[1], parts[2], sep = "_")
}

# 创建树tip.label的映射表
tree_species_map <- data.frame(
  original_tip = tree$tip.label,
  normalized_name = sapply(tree$tip.label, extract_tree_species),
  stringsAsFactors = FALSE
)

# 步骤2: 准备数据框中的物种名称
trophic_data <- trophic %>%
  mutate(
    # 确保数据中的物种名称格式一致
    normalized_name = as.character(Taxon),
    
    # 处理可能的特殊情况（如亚种标记）
    normalized_name = str_remove(normalized_name, "_ssp.*$")  # 移除亚种标记如"_sspA"
  ) 

# 步骤3: 匹配数据集和系统发育树
# 数据集中的物种名称匹配
trophic_data <- trophic_data %>%
  left_join(tree_species_map, by = "normalized_name") %>%
  mutate(
    in_tree = !is.na(original_tip),
    matched_tip = ifelse(in_tree, original_tip, NA)
  )

# 系统发育树中的物种名称匹配
tree_species_map <- tree_species_map %>%
  left_join(trophic_data %>% dplyr::select(normalized_name, in_data = Taxon),
            by = "normalized_name") %>%
  mutate(in_data = ifelse(is.na(in_data), FALSE, TRUE))

# 步骤4: 分析匹配情况
# 数据集匹配统计
matched_in_data <- sum(trophic_data$in_tree, na.rm = TRUE)
total_in_data <- nrow(trophic_data)
match_percentage_data <- round(matched_in_data / total_in_data * 100, 2)

# 系统发育树匹配统计
matched_in_tree <- sum(tree_species_map$in_data, na.rm = TRUE)
total_in_tree <- length(tree$tip.label)
match_percentage_tree <- round(matched_in_tree / total_in_tree * 100, 2)

# 步骤5: 提取匹配的物种
matched_trophic <- trophic_data %>% filter(in_tree)
matched_tips <- na.omit(trophic_data$matched_tip)

# 步骤6: 修剪系统发育树
tree_pruned <- keep.tip(tree, matched_tips)

# 步骤7: 创建最终分析数据集
# 确保数据框的行名与树的tip.label匹配
rownames(matched_trophic) <- matched_trophic$matched_tip

# 按树tip.label顺序排列数据框
matched_trophic <- matched_trophic[tree_pruned$tip.label, ]

# 步骤8: 检查修剪后树和数据的匹配
name_check_pruned <- name.check(tree_pruned, matched_trophic, data.names = rownames(matched_trophic))

writexl::write_xlsx(matched_trophic,"matched_trophic.xlsx")
write.tree(tree_pruned,"tree_pruned_trophic.tre")

#------相关性分析------
trophic <- readxl::read_excel("matched_trophic.xlsx")
tree <- read.tree("tree_pruned_trophic.tre")
colnames(trophic)[2] <- "Mass_log"
colnames(trophic)[3] <- "Trophic_level"
colnames(trophic)[16] <- "species"
results_combined_pic <- phylo_correlation_analysis(
  tree = tree,
  trait_data = trophic,
  sp_col = "species",
  x_col = "Mass_log",
  y_col = "Trophic_level"
)
results_combined_pgls <- pgls_correlation_analysis(
  tree = tree,
  trait_data = trophic,
  sp_col = "species",
  x_col = "Mass_log",
  y_col = "Trophic_level"
)

# 计算PIC
combined_pic <- data.frame(
  Mass_log = ape::pic(trophic$Mass_log, tree),
  Trophic_level = ape::pic(trophic$Trophic_level, tree)
)

# 进行线性回归分析
pearson_orig <- cor.test(combined_pic$Mass_log,combined_pic$Trophic_level, method = "pearson")
pearson_orig$estimate
pearson_orig$p.value

# 绘制散点图并添加统计信息
ggplot(combined_pic, aes(x = Mass_log, y = Trophic_level)) +
  geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
  # geom_smooth(method = "lm", se = FALSE, color = "firebrick", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Phylogenetically Independent Contrasts Analysis",
    subtitle = "all Data",
    x = "pic log(Mass_log)", 
    y = "pic Trophic_level"
  ) +
  theme_bw() +  # 白底主题
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 11),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  annotate("text", x = min(combined_pic$Mass_log) * 0.9, 
           y = max(combined_pic$Trophic_level) * 0.9,
           label = paste("pics number", nrow(combined_pic)), 
           size = 3.5, color = "black", hjust = 0)

data1 <- trophic %>% filter(Environment == "Marine")
data2 <- trophic %>% filter(Environment == "Terrestrial")
results_Marine_pic <- phylo_correlation_analysis(
  tree = tree,
  trait_data = data1,
  sp_col = "species",
  x_col = "Mass_log",
  y_col = "Trophic_level"
)
results_Marine_pgls <- pgls_correlation_analysis(
  tree = tree,
  trait_data = data1,
  sp_col = "species",
  x_col = "Mass_log",
  y_col = "Trophic_level"
)
res <- match_tree(
  tree = tree,
  trait_data = combined_data,
  sp_col = "species",
  x_col = "BodyMass_log",
  y_col = "PopDensity_log"
)

# 计算PIC
combined_pic <- data.frame(
  Mass_log = ape::pic(trophic$Mass_log, tree),
  Trophic_level = ape::pic(trophic$Trophic_level, tree)
)

# 进行线性回归分析
pearson_orig <- cor.test(combined_pic$Mass_log,combined_pic$Trophic_level, method = "pearson")
pearson_orig$estimate
pearson_orig$p.value


results_Terrestrial_pic <- phylo_correlation_analysis(
  tree = tree,
  trait_data = data2,
  sp_col = "species",
  x_col = "Mass_log",
  y_col = "Trophic_level"
)
results_Terrestrial_pgls <- pgls_correlation_analysis(
  tree = tree,
  trait_data = data2,
  sp_col = "species",
  x_col = "Mass_log",
  y_col = "Trophic_level"
)

res <- match_tree(
  tree = tree,
  trait_data = combined_data,
  sp_col = "species",
  x_col = "BodyMass_log",
  y_col = "PopDensity_log"
)

# 计算PIC
combined_pic <- data.frame(
  Mass_log = ape::pic(trophic$Mass_log, tree),
  Trophic_level = ape::pic(trophic$Trophic_level, tree)
)

# 进行线性回归分析
pearson_orig <- cor.test(combined_pic$Mass_log,combined_pic$Trophic_level, method = "pearson")
pearson_orig$estimate
pearson_orig$p.value


#---------手动识别须鲸和齿鲸----------
# 创建新列并初始化为NA
trophic$Whale_Group <- NA

# 识别须鲸 (Mysticeti)
mysticeti_genera <- c("Balaena", "Balaenoptera", "Megaptera", "Eschrichtius", "Eubalaena", "Caperea")
trophic$Whale_Group[grep(paste(mysticeti_genera, collapse = "|"), trophic$Taxon)] <- "Mysticeti"

# 识别齿鲸 (Odontoceti)
odontoceti_genera <- c("Kogia", "Berardius", "Mesoplodon", "Hyperoodon", "Ziphius", "Pontoporia", 
                       "Phocoena", "Delphinapterus", "Monodon", "Orcinus", "Globicephala", "Feresa", 
                       "Sotalia", "Tursiops", "Delphinus", "Lagenorhynchus", "Cephalorhynchus")
trophic$Whale_Group[grep(paste(odontoceti_genera, collapse = "|"), trophic$Taxon)] <- "Odontoceti"

data1 <- trophic %>% filter(Whale_Group == "Mysticeti")
data2 <- trophic %>% filter(Whale_Group == "Odontoceti")
data3 <- trophic %>% filter(Environment == "Terrestrial")
data_combine <- rbind(data2,data3)

results_Marine_pic <- phylo_correlation_analysis(
  tree = tree,
  trait_data = data_combine,
  sp_col = "species",
  x_col = "Mass_log",
  y_col = "Trophic_level"
)
results_Marine_pgls <- pgls_correlation_analysis(
  tree = tree,
  trait_data = data_combine,
  sp_col = "species",
  x_col = "Mass_log",
  y_col = "Trophic_level"
)
# 2. 准备性状数据
# 确保物种名列名正确
trait_data <- as.data.frame(trait_data)
names(trait_data)[names(trait_data) == sp_col] <- "species"

# 提取需要的数据列并删除NA值
phylo_data <- trait_data[, c("species", x_col, y_col)]
phylo_data <- na.omit(phylo_data)

# 3. 检查树与数据的匹配
tree_names <- tree$tip.label

# 查找不匹配的物种
missing_in_data <- setdiff(tree_names, phylo_data$species)
missing_in_tree <- setdiff(phylo_data$species, tree_names)

# 4. 对树和数据进行修剪对齐
# 修剪树：移除数据中不存在的物种
if(length(missing_in_data) > 0) {
  tree_trimmed <- drop.tip(tree, missing_in_data)
} else {
  tree_trimmed <- tree
}

# 修剪数据：移除树中不存在的物种
if(length(missing_in_tree) > 0) {
  phylo_data <- phylo_data[phylo_data$species %in% tree_trimmed$tip.label, ]
}

# 最终检查是否完全匹配
if(!all(phylo_data$species %in% tree_trimmed$tip.label) || 
   !all(tree_trimmed$tip.label %in% phylo_data$species)) {
  warning("Tree and data not perfectly aligned after trimming")
}

# 5. 提取两个性状向量
X1 <- phylo_data[[x_col]]
X2 <- phylo_data[[y_col]]
names(X1) <- names(X2) <- phylo_data$species

# 6. 原始数据相关性分析
pearson_orig <- cor.test(X1, X2, method = "pearson")
spearman_orig <- cor.test(X1, X2, method = "spearman")

# 7. PIC转换和分析
pic_X1 <- pic(X1, tree_trimmed)
pic_X2 <- pic(X2, tree_trimmed)

pearson_pic <- cor.test(pic_X1, pic_X2, method = "pearson")
spearman_pic <- cor.test(pic_X1, pic_X2, method = "spearman")

# 8. 创建包含所有结果的结构化列表
results <- list(
  tree = tree_trimmed,
  data = phylo_data,
  n_species = nrow(phylo_data),
  
  original_data = list(
    X1 = X1,
    X2 = X2
  ),
  
  correlations_original = list(
    pearson = pearson_orig,
    spearman = spearman_orig
  ),
  
  pics = list(
    X1_pic = pic_X1,
    X2_pic = pic_X2
  ),
  
  correlations_pic = list(
    pearson = pearson_pic,
    spearman = spearman_pic
  )
)

idx <- match(tree$tip.label, combined_data$species)

# 根据索引重新排序数据框
combined_data_reordered <- combined_data[idx, ]
rownames(combined_data_reordered) <- rep(1:324)

filtered <- combined_data_reordered[212:295,]
data <- filtered
rownames(data) <- rep(1:84)
data1 <- data[1:48,]
data2 <- data[49:84,]


# data <- all_orders_data %>%
#   filter(Order %in% c("Proboscidea","Hyracoidea","Macroscelidea","Afrosoricida")) %>%
#   mutate(Order = as.factor(Order))

data <- all_orders_data %>%
  filter(Order %in% c("Hyracoidea","Macroscelidea","Afrosoricida")) %>%
  mutate(Order = as.factor(Order))

tree <- match_tree(tree_pruned,data,"species",x_col = "BodyMass_log",
           y_col = "PopDensity_log")[[1]]
X1 <- data$BodyMass_log
X2 <- data$PopDensity_log
names(X1) <- names(X2) <- data$species
X1_pic <- pic(X1,tree)
X2_pic <- pic(X2,tree)
plot(X1_pic,X2_pic)+title("Proboscidea_Hyracoidea_Macroscelidea_Afrosoricida")
plot(X1_pic,X2_pic)+title("Hyracoidea_Macroscelidea_Afrosoricida")
cor(X1_pic,X2_pic,method = "spearman")
results1 <- phylo_correlation_analysis(
  tree = tree,
  trait_data = data,
  sp_col = "species",
  x_col = "BodyMass_log",
  y_col = "PopDensity_log"
)
results2 <- pgls_correlation_analysis(
  tree = tree_pruned,
  trait_data = data,
  sp_col = "species",
  x_col = "BodyMass_log",
  y_col = "PopDensity_log"
)
