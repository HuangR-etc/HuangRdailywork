# 保存所有关键结果的列表
key_results <- list(
  # 最终比较表格
  final_comparison = final_comparison,
  
  # 模拟次数
  n_simulations = n_simulations,
  
  # 方差分析结果
  variance_results = variance_results,
  
  # 分割分析结果
  all_random_split_results = all_random_split_results,
  all_block_split_results = all_block_split_results,
  
  # 场景颜色定义
  scenario_colors = scenario_colors,
  
  # 指标分类（便于后续绘图）
  higher_better_metrics <- c(
    "ols_r2", "adjusted_r2", "pearson_corr", "spearman_corr", "cn_smape",
    "pic_r2", "pic_pearson", "pic_spearman", "whitened_r2", "resid_r2", "lik_r2",
    "R_whole","R_ari","R_geo","R_median"
  ),
  
  lower_better_metrics <- c(
    "mse", "rmse","weighted_rmse", "mae", "median_ae", "mape", "smape",
    "pic_rmse", "pic_mae", "whitened_rmse", "whitened_mae",
    "see","percent_see"
  ),
  
  # 元信息
  save_timestamp = Sys.time(),
  n_simulations_done = n_simulations
)

# 保存到RData文件
# save(key_results, file = "time_series_analysis_results.RData")
save(key_results, file = "space_series_analysis_results.RData")
# save(key_results, file = "phy_series_analysis_results.RData")
#--------加载保存的变量---------
# 加载保存的变量
load("time_series_analysis_results.RData")
load("space_series_analysis_results.RData")
load("phy_series_analysis_results.RData")
# 重新提取必要的变量
final_comparison <- key_results$final_comparison
n_simulations <- key_results$n_simulations
variance_results <- key_results$variance_results
all_random_split_results <-  key_results$all_random_split_results
all_block_split_results <-  key_results$all_block_split_results
scenario_colors <- key_results$scenario_colors
higher_better_metrics <- key_results$higher_better_metrics
lower_better_metrics <- key_results$lower_better_metrics

# 重新创建比较均值和标准误数据框
# 创建简化的比较数据框（只包含均值）
comparison_means <- data.frame(Scenario = final_comparison$Scenario)
comparison_se <- data.frame(Scenario = final_comparison$Scenario)

metric_names <- gsub("_mean", "", grep("_mean$", names(final_comparison), value = TRUE))

for (metric in metric_names) {
  mean_col <- paste0(metric, "_mean")
  se_col <- paste0(metric, "_se")
  
  if (mean_col %in% names(final_comparison)) {
    comparison_means[[metric]] <- final_comparison[[mean_col]]
  }
  if (se_col %in% names(final_comparison)) {
    comparison_se[[metric]] <- final_comparison[[se_col]]
  }
}

# 从9.改进的可视化方案 部分开始运行所有绘图代码