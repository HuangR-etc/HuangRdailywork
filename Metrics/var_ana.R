# 整合所有数据集的排名结果
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("D:/pine/paper_redo/2025_09_obs_vs_pred/get_obs_pred/var_ana")

# 定义指标类型和颜色（与原始分析保持一致）
r2_metrics <- c("ols", "adjusted")
error_metrics <- c("mse", "rmse", "mae", "pic_rmse", "pic_mae", "whitened_rmse", "whitened_mae")
corr_metrics <- c("pearson", "spearman", "pic_pearson", "pic_spearman")
percentage_metrics <- c("smape", "cn_smape","mape")

# 定义颜色方案
type_colors <- c(
  "R-squared" = "#C73E1D",     
  "Error" = "#2E86AB",           
  "Correlation" = "#F18F01",   
  "Percentage Error" = "#A23B72"
)

# 1. 自动读取目录下所有以"variance_ranking_"开头的CSV文件
ranking_files <- list.files(pattern = "^variance_ranking_.*\\.csv$")

if(length(ranking_files) == 0) {
  stop("没有找到以'variance_ranking_'开头的CSV文件")
}

print(paste("找到", length(ranking_files), "个排名文件:"))
print(ranking_files)

# 2. 整合排名，处理指标缺失的情况
all_rankings <- data.frame()

for(file in ranking_files) {
  df <- read.csv(file)
  # 添加文件名作为数据集标识
  df$Dataset <- gsub("\\.csv$", "", file)
  all_rankings <- rbind(all_rankings, df)
}

# 3. 计算平均排名（自动处理缺失指标）
final_ranking <- all_rankings %>%
  group_by(Metric, Type) %>%
  summarise(
    Num_Datasets = n(),  # 该指标在多少个数据集中出现
    Avg_Rank = mean(Rank, na.rm = TRUE),
    Min_Rank = min(Rank, na.rm = TRUE),
    Max_Rank = max(Rank, na.rm = TRUE),
    Std_Rank = sd(Rank, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # 按平均排名排序（排名越小越好）
  arrange(Avg_Rank)

# 4. 添加排名稳定性评估
final_ranking$Stability_Score <- 1 / (1 + final_ranking$Std_Rank)

print("最终平均排名结果：")
print(final_ranking, width = Inf)

# 5. 保存最终结果
write.csv(final_ranking, "final_variance_ranking.csv", row.names = FALSE)

# 6. 创建排名热图，显示每个指标在不同数据集中的排名
rank_matrix <- all_rankings %>%
  dplyr::select(Metric, Dataset, Rank) %>%
  pivot_wider(names_from = Dataset, values_from = Rank)

print("各指标在不同数据集中的排名矩阵：")
print(rank_matrix)

# 7. 创建排名点图（改进版）- 按排名排序，显示标准差误差线
p_ranking_dot <- ggplot(final_ranking, 
                        aes(x = reorder(Metric, Avg_Rank), y = Avg_Rank, color = Type)) +
  
  # 先添加误差线（不偏移）
  geom_errorbar(aes(ymin = pmax(Avg_Rank - Std_Rank, 0),  # 确保没有负值
                    ymax = Avg_Rank + Std_Rank), 
                width = 0.3,     # 误差线宽度
                size = 0.8,      # 误差线粗细
                alpha = 0.7) +   # 透明度
  
  # 添加点，点大小与稳定性成正比（标准差的倒数）
  geom_point(aes(size = 1/(Std_Rank + 0.1)),  # 加0.1避免除零错误
             alpha = 0.9,     # 点透明度
             stroke = 1.5) +  # 点边框粗细
  
  # 添加排名数值标签 - 向上偏移避免重叠
  geom_text(aes(label = sprintf("%.2f", Avg_Rank)), 
            hjust = -0.5,     # 水平对齐（负值表示向右偏移）
            vjust = 0.5,      # 垂直对齐（0.5表示居中，正值向下，负值向上）
            size = 5,         # 字体大小
            fontface = "bold", # 粗体
            nudge_x = 0.3,    # 垂直偏移量（正值向上，负值向下）
            check_overlap = FALSE) +  # 允许标签重叠
  
  # 颜色方案 - 与学术标准匹配
  scale_color_manual(
    name = "Metric Type",     # 图例标题
    values = type_colors
  ) +
  
  # 点大小比例尺 - 表示稳定性，使用蓝色
  scale_size_continuous(
    name = "Stability\n(1/Standard Deviation)", 
    range = c(3, 8),  # 点大小范围（最小3，最大8）
    breaks = c(0.5, 1, 2, 4),  # 图例断点
    labels = c("Low", "Medium", "High", "Very High")
  ) +
  # 明确分离图例
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(override.aes = list(color = "steelblue"))
  )+
  # 坐标翻转 - 水平布局
  coord_flip() +
  
  # 标签和标题
  labs(
    title = "Variance Stability Ranking of Evaluation Metrics",
    subtitle = "Lower average rank indicates better variance stability across datasets",
    x = "Evaluation Metric", 
    y = "Average Rank ± Standard Deviation",
    caption = "Error bars represent ±1 standard deviation of rankings across multiple datasets"
  ) +
  
  # 使用黑白主题作为基础（学术标准）
  theme_bw() +
  
  # 全面主题定制，达到出版质量
  theme(
    # 标题样式 - 大号、粗体、居中
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, 
                              margin = margin(b = 10)),
    plot.subtitle = element_text(size = 16, hjust = 0.5, 
                                 margin = margin(b = 15)),
    plot.caption = element_text(size = 12, face = "italic", hjust = 1),
    
    # 坐标轴标题 - 粗体大号
    axis.title.x = element_text(size = 16, face = "bold", 
                                margin = margin(t = 15)),
    axis.title.y = element_text(size = 16, face = "bold", 
                                margin = margin(r = 15)),
    
    # 坐标轴文本 - 增大尺寸提高可读性
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    
    # 坐标轴线 - 加粗
    axis.line = element_line(color = "black", size = 1.2),
    axis.ticks = element_line(color = "black", size = 1),
    
    # 面板边框 - 加粗强调
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    
    # 网格线 - subtle但存在
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    
    # 图例样式 - 显著但不分散注意力
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(1.2, "cm"),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),
    
    # 图形边距 - 足够的留白空间
    plot.margin = margin(20, 20, 20, 20),
    
    # 面板背景
    panel.background = element_rect(fill = "white")
  ) +
  
  # Y轴标签额外定制 - 确保长指标名称完全可见
  scale_x_discrete(expand = expansion(add = 0.5)) +
  
  # 设置Y轴限制，为文本标签留出填充空间
  expand_limits(y = max(final_ranking$Avg_Rank + final_ranking$Std_Rank) * 1.15)

# 显示图形
print(p_ranking_dot)

# 8. 输出推荐结论
cat("\n=== 分析结论 ===\n")
cat("数据集数量:", length(ranking_files), "\n")
cat("总指标数量:", nrow(final_ranking), "\n\n")

# 找出表现最好的指标
top_metrics <- final_ranking %>%
  arrange(Avg_Rank) %>%
  head(10)

cat("指标方差稳定性排名（从最好到最差）：\n")
for(i in 1:nrow(top_metrics)) {
  metric <- top_metrics[i, ]
  cat(i, ". ", metric$Metric, " (", metric$Type, ")\n", sep = "")
  cat("   平均排名: ", round(metric$Avg_Rank, 2), 
      ", 排名标准差: ", round(metric$Std_Rank, 2),
      ", 稳定性得分: ", round(metric$Stability_Score, 3), "\n", sep = "")
}

# 找出排名最稳定的指标（标准差最小）
most_stable <- final_ranking %>%
  arrange(Std_Rank) %>%
  head(5)

cat("\n排名最稳定的指标（标准差最小）：\n")
for(i in 1:nrow(most_stable)) {
  metric <- most_stable[i, ]
  cat(i, ". ", metric$Metric, " (标准差: ", round(metric$Std_Rank, 3), ")\n", sep = "")
}
