te_outlier_stats <- function(pic_values) {
  q1 <- quantile(pic_values, 0.25, na.rm = TRUE)
  q3 <- quantile(pic_values, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  
  # 识别离群值并计算离群程度
  outliers <- which(pic_values < lower_bound | pic_values > upper_bound)
  outlier_magnitude <- sapply(outliers, function(i) {
    if (pic_values[i] > upper_bound) {
      magnitude <- (pic_values[i] - q3) / iqr
      return(paste0(round(magnitude, 2), "×IQR above Q3"))
    } else {
      magnitude <- (q1 - pic_values[i]) / iqr
      return(paste0(round(magnitude, 2), "×IQR below Q1"))
    }
  })
  
  return(data.frame(
    index = outliers,
    value = pic_values[outliers],
    magnitude = outlier_magnitude
  ))
}

library(caper)
library(ggplot2)
library(gridExtra)
library(ape)

#--------数据1----------
setwd("D:/pine/paper_redo/outlier")
data1 <- read.csv("example1/data/all_mammal_order.csv")[,-1]
tree1 <- read.tree("example1/data/all_mammal_order_tree.tre")
artio_data <- subset(data1, 
                     Order == "Artiodactyla", 
                     select = c("species", "BodyMass_log", "PopDensity_log"))

# 修剪系统发育树
pruned_tree <- drop.tip(tree1, 
                        tip = tree1$tip.label[!tree1$tip.label %in% artio_data$species])

# 确保数据与树完全匹配
matched_data <- artio_data[match(pruned_tree$tip.label, artio_data$species), ]
bm <- artio_data$BodyMass_log
names(bm) <- artio_data$species
pd <- artio_data$PopDensity_log
names(pd) <- artio_data$species

pic1 <- pic(bm,pruned_tree)
pic2 <- pic(pd,pruned_tree)


# 对两个性状分别计算
outliers1 <- calculate_outlier_stats(pic1)
outliers2 <- calculate_outlier_stats(pic2)

if (nrow(outliers1) > 0) {
  outliers1$trait <- "bodymass"
}
if (nrow(outliers2) > 0) {
  outliers2$trait <- "popdensity"  
}
all_outliers <- unique(c(
  which(pic1 < (quantile(pic1, 0.25, na.rm = TRUE) - 1.5 * IQR(pic1, na.rm = TRUE)) | 
          pic1 > (quantile(pic1, 0.75, na.rm = TRUE) + 1.5 * IQR(pic1, na.rm = TRUE))),
  which(pic2 < (quantile(pic2, 0.25, na.rm = TRUE) - 1.5 * IQR(pic2, na.rm = TRUE)) | 
          pic2 > (quantile(pic2, 0.75, na.rm = TRUE) + 1.5 * IQR(pic2, na.rm = TRUE)))
))

outliers <- rbind(outliers1,outliers2)

# 创建PIC散点图数据框
pic_df <- data.frame(
  bodymass_pic = pic1,
  popdensity_pic = pic2,
  is_outlier = FALSE
)

# 标记离群点
pic_df$is_outlier[all_outliers] <- TRUE
# ==================== 图1 ====================

p_enhanced <- ggplot(pic_df, aes(x = bodymass_pic, y = popdensity_pic)) +
  geom_point(aes(color = is_outlier, shape = is_outlier), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                     labels = c("FALSE" = "正常点", "TRUE" = "离群点")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "正常点", "TRUE" = "离群点")) +
  
  # 添加IQR边界线
  geom_vline(xintercept = quantile(pic1, c(0.25, 0.75), na.rm = TRUE)[1] - 1.5 * IQR(pic1, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_vline(xintercept = quantile(pic1, c(0.25, 0.75), na.rm = TRUE)[2] + 1.5 * IQR(pic1, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = quantile(pic2, c(0.25, 0.75), na.rm = TRUE)[1] - 1.5 * IQR(pic2, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = quantile(pic2, c(0.25, 0.75), na.rm = TRUE)[2] + 1.5 * IQR(pic2, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  
  # 为离群点添加index标注
  geom_text(data = subset(pic_df, is_outlier == TRUE),
            aes(label = rownames(subset(pic_df, is_outlier == TRUE))),
            vjust = -0.8, hjust = 0.5, size = 3, color = "red") +
  
  labs(title = paste0("IQR离群值\n",
                      "PICs数目: ", length(pic1)),
       x = "Body Mass PIC",
       y = "Population Density PIC",
       color = "spot", shape = "spot") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )

print(p_enhanced)
write.csv(outliers,file="Artiodactyla.csv")
#-----------数据2------------
data2 <- read.csv("Scales/data/scales2009.csv")
tree2 <- read.tree("Scales/data/tetrapods.tre")
FGSO_data <- subset(data2,
                     select = c("species", "FG", "SO"))
# 修剪系统发育树
pruned_tree <- drop.tip(tree2, 
                        tip = tree2$tip.label[!tree2$tip.label %in% FGSO_data$species])

# 确保数据与树完全匹配
matched_data <- FGSO_data[match(pruned_tree$tip.label, FGSO_data$species), ]
SO <- FGSO_data$SO
names(SO) <- FGSO_data$species
FG <- FGSO_data$FG
names(FG) <- FGSO_data$species

pic1 <- pic(SO,pruned_tree)
pic2 <- pic(FG,pruned_tree)
plot(pic1,pic2)

# 对两个性状分别计算
outliers1 <- calculate_outlier_stats(pic1)
outliers2 <- calculate_outlier_stats(pic2)

if (nrow(outliers1) > 0) {
  outliers1$trait <- "so"
}
if (nrow(outliers2) > 0) {
  outliers2$trait <- "fg"  
}
all_outliers <- unique(c(
  which(pic1 < (quantile(pic1, 0.25, na.rm = TRUE) - 1.5 * IQR(pic1, na.rm = TRUE)) | 
          pic1 > (quantile(pic1, 0.75, na.rm = TRUE) + 1.5 * IQR(pic1, na.rm = TRUE))),
  which(pic2 < (quantile(pic2, 0.25, na.rm = TRUE) - 1.5 * IQR(pic2, na.rm = TRUE)) | 
          pic2 > (quantile(pic2, 0.75, na.rm = TRUE) + 1.5 * IQR(pic2, na.rm = TRUE)))
))

outliers <- rbind(outliers1,outliers2)

# 创建PIC散点图数据框
pic_df <- data.frame(
  so_pic = pic1,
  fg_pic = pic2,
  is_outlier = FALSE
)

# 标记离群点
pic_df$is_outlier[all_outliers] <- TRUE
# ==================== 图2 ====================

p_enhanced <- ggplot(pic_df, aes(x = so_pic, y = fg_pic)) +
  geom_point(aes(color = is_outlier, shape = is_outlier), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                     labels = c("FALSE" = "正常点", "TRUE" = "离群点")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "正常点", "TRUE" = "离群点")) +
  
  # 添加IQR边界线
  geom_vline(xintercept = quantile(pic1, c(0.25, 0.75), na.rm = TRUE)[1] - 1.5 * IQR(pic1, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_vline(xintercept = quantile(pic1, c(0.25, 0.75), na.rm = TRUE)[2] + 1.5 * IQR(pic1, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = quantile(pic2, c(0.25, 0.75), na.rm = TRUE)[1] - 1.5 * IQR(pic2, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = quantile(pic2, c(0.25, 0.75), na.rm = TRUE)[2] + 1.5 * IQR(pic2, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  
  # 为离群点添加index标注
  geom_text(data = subset(pic_df, is_outlier == TRUE),
            aes(label = rownames(subset(pic_df, is_outlier == TRUE))),
            vjust = -0.8, hjust = 0.5, size = 3, color = "red") +
  
  labs(title = paste0("IQR离群值\n",
                      "PICs数目: ", length(pic1)),
       x = "SO PIC",
       y = "FG PIC",
       color = "spot", shape = "spot") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )

print(p_enhanced)
write.csv(outliers,file="FGSO.csv")
#--------数据3---------
data3 <- read.csv("Makino/results/example2.csv")
tree3 <- read.tree("Makino/data/full_species_raxml_tree.ph")
size_data <- subset(data3,
                    select = c("genome", "log_GS", "log_PS"))
# 修剪系统发育树
pruned_tree <- drop.tip(tree3, 
                        tip = tree3$tip.label[!tree3$tip.label %in% size_data$genome])
pruned_tree$node.label <- NULL
# 确保数据与树完全匹配
matched_data <- size_data[match(pruned_tree$tip.label, size_data$genome), ]
pic1 <- pic(matched_data$log_GS, pruned_tree)
pic2 <- pic(matched_data$log_PS, pruned_tree)
plot(pic1,pic2)
# 对两个性状分别计算
outliers1 <- calculate_outlier_stats(pic1)
outliers2 <- calculate_outlier_stats(pic2)

if (nrow(outliers1) > 0) {
  outliers1$trait <- "GS"
}
if (nrow(outliers2) > 0) {
  outliers2$trait <- "PS"  
}
all_outliers <- unique(c(
  which(pic1 < (quantile(pic1, 0.25, na.rm = TRUE) - 1.5 * IQR(pic1, na.rm = TRUE)) | 
          pic1 > (quantile(pic1, 0.75, na.rm = TRUE) + 1.5 * IQR(pic1, na.rm = TRUE))),
  which(pic2 < (quantile(pic2, 0.25, na.rm = TRUE) - 1.5 * IQR(pic2, na.rm = TRUE)) | 
          pic2 > (quantile(pic2, 0.75, na.rm = TRUE) + 1.5 * IQR(pic2, na.rm = TRUE)))
))

outliers <- rbind(outliers1,outliers2)

# 创建PIC散点图数据框
pic_df <- data.frame(
  gs_pic = pic1,
  ps_pic = pic2,
  is_outlier = FALSE
)

# 标记离群点
pic_df$is_outlier[all_outliers] <- TRUE

# ==================== 图3 ====================

p_enhanced <- ggplot(pic_df, aes(x = gs_pic, y = ps_pic)) +
  geom_point(aes(color = is_outlier, shape = is_outlier), size = 3, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                     labels = c("FALSE" = "正常点", "TRUE" = "离群点")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "正常点", "TRUE" = "离群点")) +
  
  # 添加IQR边界线
  geom_vline(xintercept = quantile(pic1, c(0.25, 0.75), na.rm = TRUE)[1] - 1.5 * IQR(pic1, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_vline(xintercept = quantile(pic1, c(0.25, 0.75), na.rm = TRUE)[2] + 1.5 * IQR(pic1, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = quantile(pic2, c(0.25, 0.75), na.rm = TRUE)[1] - 1.5 * IQR(pic2, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = quantile(pic2, c(0.25, 0.75), na.rm = TRUE)[2] + 1.5 * IQR(pic2, na.rm = TRUE), 
             linetype = "dashed", color = "blue", alpha = 0.5) +
  
  # 为离群点添加index标注
  geom_text(data = subset(pic_df, is_outlier == TRUE),
            aes(label = rownames(subset(pic_df, is_outlier == TRUE))),
            vjust = -0.8, hjust = 0.5, size = 3, color = "red") +
  
  labs(title = paste0("IQR离群值\n",
                      "PICs数目: ", length(pic1)),
       x = "log_GenomeSize PIC",
       y = "log_ProgaguleSize PIC",
       color = "spot", shape = "spot") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )

print(p_enhanced)
write.csv(outliers,file="size.csv")

