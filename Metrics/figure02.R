# 加载必要的包
library(ggplot2)
library(dplyr)
library(scales)

# 创建数据框
data <- data.frame(
  category = c("是", "否"),
  count = c(3, 106)
) %>%
  mutate(
    percentage = count/sum(count) * 100,
    label = paste0(category, "\n", count, " (", round(percentage, 1), "%)")
  )

# 设置Science期刊风格的颜色
science_colors <- c("是" = "#E69F00", "否" = "#0072B2")

# 创建饼图
p <- ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.8) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = science_colors) +
  theme_void(base_size = 24) +  # 增大基础字体
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 20)),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),  # 增大图例文字
    legend.key.size = unit(1.5, "lines"),  # 增加图例键大小
    legend.spacing.y = unit(1, "cm"),   # 增加图例项之间的垂直间距
    legend.margin = margin(t = 10, r = 10, b = 10, l = 10),  # 增加图例外边距
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "文献是否在评估时考虑系统发育关系",
    fill = "类别"
  ) +
  guides(fill = guide_legend(byrow = TRUE))

# 为大的扇形添加白色标签（增大字体）
p <- p + geom_text(
  aes(label = ifelse(count > 5, label, "")),
  position = position_stack(vjust = 0.5),
  size = 8,  # 增大字体
  color = "white",
  fontface = "bold"
)

# 为小的扇形添加带背景的标签和引线（增大字体）
p <- p + 
  # 添加引线
  annotate("segment",
           x = 0.8, xend = 1.2,
           y = 1.5, yend = 2.5,
           color = "black", linewidth = 0.5) +
  # 添加带背景的标签
  annotate("label",
           x = 1.3, y = 2.5,
           label = data$label[1],
           size = 8,  # 增大字体
           color = "black",
           fill = "white",
           label.padding = unit(0.5, "lines"),
           label.size = 0.8,
           fontface = "bold")

# 显示图形
print(p)

# 保存为高分辨率图片
ggsave("phylogenetic_analysis_piechart.png",
       plot = p,
       device = "png",
       width = 10,
       height = 8,
       dpi = 600)