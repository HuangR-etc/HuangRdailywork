### 看一下当前工作目录
setwd("/home/huangr/projects/2026_04_distance/outputs/tables")
getwd()

### ================================
### 2. 读入四个原始数据文件
### ================================
df_small_balanced <- read.csv("result3_signal_small_balanced.csv", stringsAsFactors = FALSE)
df_small_ladder   <- read.csv("result3_signal_small_ladder.csv", stringsAsFactors = FALSE)
df_large_balanced <- read.csv("result3_signal_large_balanced.csv", stringsAsFactors = FALSE)
df_large_ladder   <- read.csv("result3_signal_large_ladder.csv", stringsAsFactors = FALSE)

### 给每个数据框加一个 Tree 标签，方便后面汇总
df_small_balanced$Tree <- "small_balanced"
df_small_ladder$Tree   <- "small_ladder"
df_large_balanced$Tree <- "large_balanced"
df_large_ladder$Tree   <- "large_ladder"

### 合并成一个总表
signal_df <- rbind(
  df_small_balanced,
  df_small_ladder,
  df_large_balanced,
  df_large_ladder
)

### ================================
### 3. 检查列名
### ================================
colnames(signal_df)

### ================================
### 4. 简单检查 Subset_Type 和 Model
### ================================
unique(signal_df$Subset_Type)
unique(signal_df$Model)

### ================================
### 5. 定义一个经验 p 值计算函数
### 逻辑尽量和你原来的 result3 分析保持一致
### ================================
calculate_empirical_pvalues_manual <- function(signal_df, tree_name, model_name, metric_name) {
  
  dat <- signal_df[signal_df$Tree == tree_name & signal_df$Model == model_name, ]
  
  if (nrow(dat) == 0) {
    return(NULL)
  }
  
  dispersed_vals <- dat[[metric_name]][dat$Subset_Type == "dispersed"]
  clustered_vals <- dat[[metric_name]][dat$Subset_Type == "clustered"]
  random_vals    <- dat[[metric_name]][dat$Subset_Type == "random"]
  
  ### 去掉 NA
  dispersed_vals <- dispersed_vals[!is.na(dispersed_vals)]
  clustered_vals <- clustered_vals[!is.na(clustered_vals)]
  random_vals    <- random_vals[!is.na(random_vals)]
  
  tol <- 1e-15
  
  ### ----------------------------
  ### 情况 1：Moran's I
  ### 按绝对值比较，|I| 越大表示自相关越强
  ### 预期：clustered > random > dispersed（按绝对值）
  ### 所以：
  ### dispersed vs clustered：检验 |dispersed| < |clustered|
  ### dispersed vs random   ：检验 |dispersed| < |random|
  ### random vs clustered   ：检验 |random| < |clustered|
  ### ----------------------------
  if (metric_name == "MoransI") {
    
    dispersed_cmp <- abs(dispersed_vals)
    clustered_cmp <- abs(clustered_vals)
    random_cmp    <- abs(random_vals)
    
    mean_dispersed <- mean(dispersed_vals, na.rm = TRUE)
    mean_clustered <- mean(clustered_vals, na.rm = TRUE)
    mean_random    <- mean(random_vals, na.rm = TRUE)
    
    if (length(dispersed_cmp) > 0 && length(clustered_cmp) > 0) {
      count <- 0
      for (d in dispersed_cmp) {
        count <- count + sum(clustered_cmp <= d | abs(clustered_cmp - d) < tol)
      }
      count <- count / length(dispersed_cmp)
      p_dispersed_vs_clustered <- (count + 1) / (length(clustered_cmp) + 1)
    } else {
      p_dispersed_vs_clustered <- NA
    }
    
    if (length(dispersed_cmp) > 0 && length(random_cmp) > 0) {
      count <- 0
      for (d in dispersed_cmp) {
        count <- count + sum(random_cmp <= d | abs(random_cmp - d) < tol)
      }
      count <- count / length(dispersed_cmp)
      p_dispersed_vs_random <- (count + 1) / (length(random_cmp) + 1)
    } else {
      p_dispersed_vs_random <- NA
    }
    
    if (length(random_cmp) > 0 && length(clustered_cmp) > 0) {
      total_comparisons <- length(random_cmp) * length(clustered_cmp)
      count <- 0
      for (r in random_cmp) {
        count <- count + sum(clustered_cmp <= r | abs(clustered_cmp - r) < tol)
      }
      p_random_vs_clustered <- (count + 1) / (total_comparisons + 1)
    } else {
      p_random_vs_clustered <- NA
    }
  }
  
  ### ----------------------------
  ### 情况 2：ESS
  ### ess_1 / ess_2 越大表示越接近独立
  ### 预期：dispersed > random > clustered
  ### 所以：
  ### dispersed vs clustered：检验 dispersed > clustered
  ### dispersed vs random   ：检验 dispersed > random
  ### random vs clustered   ：检验 random > clustered
  ### ----------------------------
  else if (metric_name %in% c("ess_1", "ess_2")) {
    
    mean_dispersed <- mean(dispersed_vals, na.rm = TRUE)
    mean_clustered <- mean(clustered_vals, na.rm = TRUE)
    mean_random    <- mean(random_vals, na.rm = TRUE)
    
    if (length(dispersed_vals) > 0 && length(clustered_vals) > 0) {
      count <- 0
      for (d in dispersed_vals) {
        count <- count + sum(clustered_vals >= d | abs(clustered_vals - d) < tol)
      }
      count <- count / length(dispersed_vals)
      p_dispersed_vs_clustered <- (count + 1) / (length(clustered_vals) + 1)
    } else {
      p_dispersed_vs_clustered <- NA
    }
    
    if (length(dispersed_vals) > 0 && length(random_vals) > 0) {
      count <- 0
      for (d in dispersed_vals) {
        count <- count + sum(random_vals >= d | abs(random_vals - d) < tol)
      }
      count <- count / length(dispersed_vals)
      p_dispersed_vs_random <- (count + 1) / (length(random_vals) + 1)
    } else {
      p_dispersed_vs_random <- NA
    }
    
    if (length(random_vals) > 0 && length(clustered_vals) > 0) {
      total_comparisons <- length(random_vals) * length(clustered_vals)
      count <- 0
      for (r in random_vals) {
        count <- count + sum(clustered_vals >= r | abs(clustered_vals - r) < tol)
      }
      p_random_vs_clustered <- (count + 1) / (total_comparisons + 1)
    } else {
      p_random_vs_clustered <- NA
    }
  }
  
  ### ----------------------------
  ### 情况 3：其余指标
  ### K / Lambda / mpnns / mean_offdiag_cor / max_offdiag_cor
  ### 数值越大表示信号/依赖越强
  ### 预期：clustered > random > dispersed
  ### 所以：
  ### dispersed vs clustered：检验 dispersed < clustered
  ### dispersed vs random   ：检验 dispersed < random
  ### random vs clustered   ：检验 random < clustered
  ### ----------------------------
  else {
    
    mean_dispersed <- mean(dispersed_vals, na.rm = TRUE)
    mean_clustered <- mean(clustered_vals, na.rm = TRUE)
    mean_random    <- mean(random_vals, na.rm = TRUE)
    
    if (length(dispersed_vals) > 0 && length(clustered_vals) > 0) {
      count <- 0
      for (d in dispersed_vals) {
        count <- count + sum(clustered_vals <= d | abs(clustered_vals - d) < tol)
      }
      count <- count / length(dispersed_vals)
      p_dispersed_vs_clustered <- (count + 1) / (length(clustered_vals) + 1)
    } else {
      p_dispersed_vs_clustered <- NA
    }
    
    if (length(dispersed_vals) > 0 && length(random_vals) > 0) {
      count <- 0
      for (d in dispersed_vals) {
        count <- count + sum(random_vals <= d | abs(random_vals - d) < tol)
      }
      count <- count / length(dispersed_vals)
      p_dispersed_vs_random <- (count + 1) / (length(random_vals) + 1)
    } else {
      p_dispersed_vs_random <- NA
    }
    
    if (length(random_vals) > 0 && length(clustered_vals) > 0) {
      total_comparisons <- length(random_vals) * length(clustered_vals)
      count <- 0
      for (r in random_vals) {
        count <- count + sum(clustered_vals <= r | abs(clustered_vals - r) < tol)
      }
      p_random_vs_clustered <- (count + 1) / (total_comparisons + 1)
    } else {
      p_random_vs_clustered <- NA
    }
  }
  
  out <- data.frame(
    Tree = tree_name,
    Model = model_name,
    Metric = metric_name,
    Mean_Dispersed = mean_dispersed,
    Mean_Clustered = mean_clustered,
    Mean_Random = mean_random,
    P_Dispersed_vs_Clustered = p_dispersed_vs_clustered,
    P_Dispersed_vs_Random = p_dispersed_vs_random,
    P_Random_vs_Clustered = p_random_vs_clustered,
    N_Dispersed = length(dispersed_vals),
    N_Clustered = length(clustered_vals),
    N_Random = length(random_vals),
    stringsAsFactors = FALSE
  )
  
  return(out)
}

### ================================
### 6. 指定要计算的指标
### ================================
metrics_to_test <- c(
  "K",
  "Lambda",
  "MoransI",
  "mpnns",
  "mean_offdiag_cor",
  "max_offdiag_cor",
  "ess_1",
  "ess_2"
)

### ================================
### 7. 提取所有 Tree 和 Model
### ================================
all_trees  <- unique(signal_df$Tree)
all_models <- unique(signal_df$Model)

all_trees
all_models

### ================================
### 8. 循环计算所有 Tree × Model × Metric 的经验 p 值
### ================================
empirical_pvalues_all <- data.frame()

for (tree_name in all_trees) {
  for (model_name in all_models) {
    for (metric_name in metrics_to_test) {
      
      res <- calculate_empirical_pvalues_manual(
        signal_df   = signal_df,
        tree_name   = tree_name,
        model_name  = model_name,
        metric_name = metric_name
      )
      
      if (!is.null(res)) {
        empirical_pvalues_all <- rbind(empirical_pvalues_all, res)
      }
    }
  }
}

### 看结果前几行
head(empirical_pvalues_all)

### ================================
### 9. 查看完整结果
### ================================
empirical_pvalues_all

### ================================
### 10. 按 Tree 和 Model 排序，便于阅读
### ================================
empirical_pvalues_all <- empirical_pvalues_all[
  order(empirical_pvalues_all$Tree,
        empirical_pvalues_all$Model,
        empirical_pvalues_all$Metric),
]

empirical_pvalues_all

### ================================
### 11. 写出结果到 csv
### ================================
write.csv(
  empirical_pvalues_all,
  file = "manual_empirical_pvalues_all.csv",
  row.names = FALSE
)

### ================================
### 12. 如果你只想看某一棵树，例如 large_balanced
### ================================
subset(empirical_pvalues_all, Tree == "large_balanced")

### ================================
### 13. 如果你只想看某个指标，例如 ess_1
### ================================
subset(empirical_pvalues_all, Metric == "ess_1")

### ================================
### 14. 如果你想再检查每个组的均值排序
### 例如看某个 Tree × Model
### ================================
subset(
  empirical_pvalues_all,
  Tree == "small_balanced" & Model == "Lambda1"
)

### ================================
### 15. 可选：加一个“是否符合预期排序”的检查
### 对不同指标方向分别判断
### ================================
check_pattern <- function(mean_d, mean_r, mean_c, metric) {
  
  tol <- 1e-10
  
  if (metric == "MoransI") {
    return( abs(mean_c) > abs(mean_r) - tol && abs(mean_r) > abs(mean_d) - tol )
  }
  
  if (metric %in% c("ess_1", "ess_2")) {
    return( mean_d > mean_r - tol && mean_r > mean_c - tol )
  }
  
  return( mean_c > mean_r - tol && mean_r > mean_d - tol )
}

empirical_pvalues_all$Pattern_Holds <- mapply(
  check_pattern,
  empirical_pvalues_all$Mean_Dispersed,
  empirical_pvalues_all$Mean_Random,
  empirical_pvalues_all$Mean_Clustered,
  empirical_pvalues_all$Metric
)

empirical_pvalues_all

### ================================
### 16. 把带 pattern 的结果也保存
### ================================
write.csv(
  empirical_pvalues_all,
  file = "empirical_pvalues_with_pattern.csv",
  row.names = FALSE
)

library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(ggrepel)

### 如果 empirical_pvalues_all 已经在环境里，就直接用
df <- empirical_pvalues_all

### 如果没有，可以手动读入
# df <- read.csv("manual_empirical_pvalues_with_pattern.csv", stringsAsFactors = FALSE)

### 只保留你现在关心的列
df_plot <- df %>%
  select(Tree, Model, Metric, P_Dispersed_vs_Random, Pattern_Holds)

### family 列
df_plot <- df_plot %>%
  mutate(
    Family = case_when(
      str_detect(Model, "^Lambda") ~ "Lambda",
      str_detect(Model, "^OU_") ~ "OU",
      TRUE ~ "Other"
    )
  )

### 模型短标签（图上更紧凑）
df_plot <- df_plot %>%
  mutate(
    Model_short = case_when(
      Model == "Lambda1" ~ "L1",
      Model == "Lambda0.7" ~ "L0.7",
      Model == "OU_alpha0.2" ~ "OU 0.2",
      Model == "OU_alpha1" ~ "OU 1",
      Model == "OU_alpha5" ~ "OU 5",
      TRUE ~ Model
    )
  )

### metric 顺序：建议先 dependence diagnostics，再 signal metrics
metric_order <- c(
  "mean_offdiag_cor",
  "max_offdiag_cor",
  "ess_1",
  "ess_2",
  "K",
  "Lambda",
  "MoransI",
  "mpnns"
)

tree_order <- c(
  "large_balanced",
  "large_ladder",
  "small_balanced",
  "small_ladder"
)

family_order <- c("Lambda", "OU")

df_plot <- df_plot %>%
  mutate(
    Tree = factor(Tree, levels = tree_order),
    Family = factor(Family, levels = family_order),
    Metric = factor(Metric, levels = metric_order)
  )

### 给 metric 换成更适合图上的标签
metric_labels <- c(
  mean_offdiag_cor = "Mean offdiag cor",
  max_offdiag_cor  = "Max offdiag cor",
  ess_1            = "ESS 1",
  ess_2            = "ESS 2",
  K                = "K",
  Lambda           = "Lambda",
  MoransI          = "Moran's I",
  mpnns            = "mpnns"
)

### p 值标签函数
fmt_p <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(
      x < 0.001,
      format(x, scientific = TRUE, digits = 2),
      sprintf("%.3f", x)
    )
  )
}

df_plot <- df_plot %>%
  mutate(
    p_label = fmt_p(P_Dispersed_vs_Random),
    neglog10_p = -log10(pmax(P_Dispersed_vs_Random, 1e-300))
  )
p1 <- ggplot(df_plot, aes(x = Metric, y = Model_short, fill = neglog10_p)) +
  geom_tile(color = "grey85", linewidth = 0.4) +
  geom_text(aes(label = p_label), size = 3.1) +
  facet_grid(Tree ~ Family, scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = metric_labels) +
  scale_fill_viridis_c(
    option = "C",
    name = expression(-log[10](p)),
    direction = 1
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "P_Dispersed_vs_Random across trees, metrices, and model families",
    subtitle = "Each cell is one actual p-value; no averaging"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "right"
  )


