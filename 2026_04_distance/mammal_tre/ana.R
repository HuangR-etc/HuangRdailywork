# 加载ape包
library(ape)
setwd("D:/mammal_tre")
# 读取树文件
tree <- read.nexus("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre")
tax <- read.csv("taxonomy_mamPhy_5911species_toPublish.csv")
sampled <- tax$tiplabel[tax$samp == "sampled"]

setdiff(tree$tip.label, sampled)
# "_Anolis_carolinensis"
# "Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA"

drop.tips <- c(
  "_Anolis_carolinensis",
  "Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA"
)

tr_4098 <- drop.tip(tree, drop.tips)

Ntip(tr_4098)

tax_carn <- tax[
  toupper(tax[["ord"]]) == "CARNIVORA" &
    tax[["extinct."]] == 0,
]

# 进一步保留 sampled 物种
if ("samp" %in% names(tax)) {
  tax_carn <- tax_carn[tax_carn$samp == "sampled", ]
}

# 6. 与树的 tip 做交集
candidate_tips <- intersect(tax_carn[["tiplabel"]], tree$tip.label)
carn_tree <- drop.tip(tree, setdiff(tree$tip.label, candidate_tips))

# 8. 候选池大小，即文中的 [N]
N <- Ntip(carn_tree)
#Carnivora食肉目：261物种
N


# 提取鼠科 Muridae
tax_mur <- tax[
  toupper(tax[["fam"]]) == "MURIDAE" &
    tax[["extinct."]] == 0,
]

# 如果有 samp 字段，只保留 sampled 物种
if ("samp" %in% names(tax)) {
  tax_mur <- tax_mur[tax_mur$samp == "sampled", ]
}

# 与清理后的 4098-tip 树取交集
mur_tips <- intersect(tax_mur[["tiplabel"]], tr_4098$tip.label)

# 提取鼠科树
mur_tree <- drop.tip(tr_4098, setdiff(tr_4098$tip.label, mur_tips))

# 查看鼠科物种数
Ntip(mur_tree)
source("D:/R/distance_metrics.R")
source("D:/R/objective_compare.R")
source("D:/R/subset_greedy.R")
source("D:/R/subset_exchange.R")
source("D:/R/subset_random.R")
# ============================================================
# NEW FUNCTIONS:
# 新版 empirical analysis 所需函数
# ============================================================

# ---------- 1. MaxPD ----------
calc_maxpd <- function(dist_mat, subset) {
  if (length(subset) < 2) return(0)
  
  sub_dist <- dist_mat[subset, subset]
  upper_values <- sub_dist[upper.tri(sub_dist)]
  
  max(upper_values, na.rm = TRUE)
}

# ---------- 2. 扩展指标：MinPD, MeanPD, MeanNND, MaxPD ----------
calc_subset_metrics_extended <- function(dist_mat, subset) {
  base_metrics <- calc_subset_metrics(dist_mat, subset)
  
  list(
    MinPD = base_metrics$MinPD,
    MeanPD = base_metrics$MeanPD,
    MeanNND = base_metrics$MeanNND,
    MaxPD = calc_maxpd(dist_mat, subset)
  )
}

# ---------- 3. 多个 subset 的扩展指标 ----------
calc_multiple_subsets_metrics_extended <- function(dist_mat, subsets) {
  
  out <- data.frame(
    SubsetID = seq_along(subsets),
    MinPD = NA_real_,
    MeanPD = NA_real_,
    MeanNND = NA_real_,
    MaxPD = NA_real_
  )
  
  for (i in seq_along(subsets)) {
    m <- calc_subset_metrics_extended(dist_mat, subsets[[i]])
    
    out$MinPD[i] <- m$MinPD
    out$MeanPD[i] <- m$MeanPD
    out$MeanNND[i] <- m$MeanNND
    out$MaxPD[i] <- m$MaxPD
  }
  
  out
}

# ---------- 4. 分散子集起始物种：平均距离最大 ----------
find_max_mean_distance_species <- function(dist_obj) {
  
  d <- dist_obj$dist_mat
  
  diag(d) <- NA
  
  mean_dist <- rowMeans(d, na.rm = TRUE)
  
  start_idx <- which.max(mean_dist)
  
  list(
    start_idx = start_idx,
    start_name = dist_obj$tip_labels[start_idx],
    mean_distance = mean_dist[start_idx],
    all_mean_distances = mean_dist
  )
}

# ---------- 5. 从固定起始物种开始 greedy 构造分散子集 ----------
build_dispersed_greedy_fixed_start <- function(dist_obj, subset_size) {
  
  all_tips <- seq_along(dist_obj$tip_labels)
  
  start_info <- find_max_mean_distance_species(dist_obj)
  
  current_subset <- start_info$start_idx
  available <- setdiff(all_tips, current_subset)
  
  while (length(current_subset) < subset_size) {
    
    best_candidate <- NULL
    best_metrics <- NULL
    
    for (candidate in available) {
      
      temp_subset <- c(current_subset, candidate)
      temp_metrics <- calc_subset_metrics(dist_obj$dist_mat, temp_subset)
      
      if (is.null(best_candidate)) {
        best_candidate <- candidate
        best_metrics <- temp_metrics
      } else {
        if (is_better_lexico_max(temp_metrics, best_metrics)) {
          best_candidate <- candidate
          best_metrics <- temp_metrics
        }
      }
    }
    
    current_subset <- c(current_subset, best_candidate)
    available <- setdiff(all_tips, current_subset)
  }
  
  final_metrics <- calc_subset_metrics(dist_obj$dist_mat, current_subset)
  
  list(
    subset = current_subset,
    subset_names = dist_obj$tip_labels[current_subset],
    metrics = final_metrics,
    start_info = start_info,
    algorithm = "dispersed_fixed_max_mean_distance_greedy"
  )
}

# ---------- 6. 新版分散子集完整算法：固定起点 + greedy + exchange ----------
run_dispersed_fixed_start_algorithm <- function(dist_obj, subset_size) {
  
  cat("Running dispersed fixed-start algorithm...\n")
  
  cat("  Phase 0: choosing fixed start species\n")
  start_info <- find_max_mean_distance_species(dist_obj)
  cat("    Start species:", start_info$start_name, "\n")
  cat("    Mean distance:", start_info$mean_distance, "\n")
  
  cat("  Phase 1: greedy forward selection\n")
  greedy_result <- build_dispersed_greedy_fixed_start(
    dist_obj = dist_obj,
    subset_size = subset_size
  )
  
  cat("    Greedy result: MinPD =", greedy_result$metrics$MinPD,
      "MeanPD =", greedy_result$metrics$MeanPD,
      "MeanNND =", greedy_result$metrics$MeanNND, "\n")
  
  cat("  Phase 2: exchange refinement\n")
  exchange_result <- refine_subset_exchange(
    dist_obj = dist_obj,
    current_subset = greedy_result$subset,
    maximize = TRUE,
    single_objective = NULL
  )
  
  cat("    Exchange result: MinPD =", exchange_result$metrics$MinPD,
      "MeanPD =", exchange_result$metrics$MeanPD,
      "MeanNND =", exchange_result$metrics$MeanNND, "\n")
  
  list(
    start_info = start_info,
    greedy_result = greedy_result,
    exchange_result = exchange_result,
    final_subset = exchange_result$subset,
    final_subset_names = exchange_result$subset_names,
    final_metrics = exchange_result$metrics,
    algorithm = "dispersed_fixed_start_greedy_exchange"
  )
}

# ---------- 7. 聚集子集候选生成：每个物种作为种子，找最近 s - 1 个邻居 ----------
generate_clustered_seed_neighborhoods <- function(dist_obj, subset_size) {
  
  d <- dist_obj$dist_mat
  all_tips <- seq_along(dist_obj$tip_labels)
  
  candidates <- vector("list", length(all_tips))
  
  for (seed in all_tips) {
    
    distances_from_seed <- d[seed, ]
    distances_from_seed[seed] <- Inf
    
    nearest_neighbors <- order(distances_from_seed)[1:(subset_size - 1)]
    
    subset_now <- c(seed, nearest_neighbors)
    
    candidates[[seed]] <- sort(unique(subset_now))
  }
  
  # 去重：不同 seed 可能得到同一个 subset
  keys <- sapply(candidates, function(x) paste(sort(x), collapse = "|"))
  
  unique_keys <- unique(keys)
  
  unique_candidates <- candidates[match(unique_keys, keys)]
  
  seed_id_for_unique <- match(unique_keys, keys)
  
  list(
    subsets = unique_candidates,
    seed_idx = seed_id_for_unique,
    seed_names = dist_obj$tip_labels[seed_id_for_unique],
    n_raw = length(candidates),
    n_unique = length(unique_candidates)
  )
}

# ---------- 8. 聚集子集排序：MeanPD min -> MeanNND min -> MaxPD min ----------
select_best_clustered_neighborhood <- function(dist_obj, subset_size) {
  
  cand <- generate_clustered_seed_neighborhoods(
    dist_obj = dist_obj,
    subset_size = subset_size
  )
  
  metrics_df <- calc_multiple_subsets_metrics_extended(
    dist_mat = dist_obj$dist_mat,
    subsets = cand$subsets
  )
  
  metrics_df$Seed_Idx <- cand$seed_idx
  metrics_df$Seed_Name <- cand$seed_names
  
  # 聚集子集新版排序：
  # 1. MeanPD 最小
  # 2. MeanNND 最小
  # 3. MaxPD 最小
  ord <- order(
    metrics_df$MeanPD,
    metrics_df$MeanNND,
    metrics_df$MaxPD
  )
  
  best_row <- metrics_df[ord[1], ]
  best_subset <- cand$subsets[[ord[1]]]
  
  list(
    final_subset = best_subset,
    final_subset_names = dist_obj$tip_labels[best_subset],
    final_metrics = calc_subset_metrics_extended(dist_obj$dist_mat, best_subset),
    candidate_metrics = metrics_df,
    ordered_candidate_metrics = metrics_df[ord, ],
    n_raw_candidates = cand$n_raw,
    n_unique_candidates = cand$n_unique,
    best_seed_idx = best_row$Seed_Idx,
    best_seed_name = best_row$Seed_Name,
    algorithm = "clustered_seed_nearest_neighbors_meanpd_meannnd_maxpd"
  )
}

# ---------- 9. p value 和 SES ----------
calc_p_high <- function(obs, null_values) {
  (1 + sum(null_values >= obs)) / (1 + length(null_values))
}

calc_p_low <- function(obs, null_values) {
  (1 + sum(null_values <= obs)) / (1 + length(null_values))
}

calc_ses <- function(obs, null_values) {
  (obs - mean(null_values)) / sd(null_values)
}
# ============================================================
# 5. 基本检查
# ============================================================

is.rooted(carn_tree)
is.ultrametric(carn_tree, tol = 1e-7)
#在允许1e-6量级的误差时，是超度量树
is.ultrametric(carn_tree, tol = 1e-6)

# 如果 ultrametric 有数值误差，可先不处理；
# 这里只是检查，不建议随意强制修改 branch length。
summary(carn_tree$edge.length)

# 保存候选池
write.csv(
  data.frame(tiplabel = carn_tree$tip.label),
  "Carnivora_candidate_pool_261_DNAonly_extant_sampled.csv",
  row.names = FALSE
)

# ============================================================
# 6. 建立距离对象
# ============================================================

dist_obj <- create_distance_object(carn_tree)

length(dist_obj$tip_labels)
dim(dist_obj$dist_mat)

all(rownames(dist_obj$dist_mat) == dist_obj$tip_labels)
all(colnames(dist_obj$dist_mat) == dist_obj$tip_labels)

# ============================================================
# 7. 设置经验分析参数
# ============================================================

subset_size <- 20

# 调试时用 100 或 1000；正式结果建议 10000
n_null_reps <- 1000

# greedy 随机起点次数；调试可设 5，正式可设 50 或 100
n_greedy_starts <- 5

subset_size < N

# ============================================================
# 8. 构造 phylogenetically dispersed subset
# NEW: fixed start + greedy + exchange
# ============================================================

disp_result <- run_dispersed_fixed_start_algorithm(
  dist_obj = dist_obj,
  subset_size = subset_size
)

disp_subset_idx <- disp_result$final_subset
disp_subset_names <- disp_result$final_subset_names
disp_metrics <- disp_result$final_metrics

length(disp_subset_idx)
disp_metrics
disp_subset_names

# 检查起始物种
disp_result$start_info$start_name
disp_result$start_info$mean_distance
# ============================================================
# 9. 构造 phylogenetically clustered subset
# NEW: seed species + nearest s - 1 neighbors
# ============================================================

clust_result <- select_best_clustered_neighborhood(
  dist_obj = dist_obj,
  subset_size = subset_size
)

clust_subset_idx <- clust_result$final_subset
clust_subset_names <- clust_result$final_subset_names
clust_metrics <- clust_result$final_metrics

length(clust_subset_idx)
clust_metrics
clust_subset_names

# 查看候选聚集子集数量
clust_result$n_raw_candidates
clust_result$n_unique_candidates

# 查看最优聚集子集的 seed
clust_result$best_seed_name

# 查看前 10 个最聚集候选
head(clust_result$ordered_candidate_metrics, 10)

write.csv(
  clust_result$candidate_metrics,
  paste0("Carnivora_clustered_seed_candidate_metrics_s", subset_size, ".csv"),
  row.names = FALSE
)
# ============================================================
# 10. 生成同一个 random baseline
# ============================================================

set.seed(20260428)

random_subsets_idx <- sample_random_subsets(
  dist_obj = dist_obj,
  subset_size = subset_size,
  n_reps = n_null_reps,
  replace = FALSE
)

length(random_subsets_idx)
length(random_subsets_idx[[1]])


# ============================================================
# 11. 计算 random baseline 的距离指标
# NEW: include MaxPD
# ============================================================

random_dist_metrics <- calc_multiple_subsets_metrics_extended(
  dist_mat = dist_obj$dist_mat,
  subsets = random_subsets_idx
)

head(random_dist_metrics)
summary(random_dist_metrics)

write.csv(
  random_dist_metrics,
  paste0("Carnivora_random_distance_metrics_s", subset_size, ".csv"),
  row.names = FALSE
)

# ============================================================
# 12. 检查 observed distance metrics
# ============================================================

disp_metrics <- calc_subset_metrics_extended(dist_obj$dist_mat, disp_subset_idx)
clust_metrics <- calc_subset_metrics_extended(dist_obj$dist_mat, clust_subset_idx)

disp_metrics
clust_metrics

# ============================================================
# 13A. Dispersed subset vs random baseline
# 指标越大越分散：high-tail p-value
# ============================================================

disp_distance_summary <- data.frame(
  Subset_Type = "dispersed",
  Metric = c("MinPD", "MeanPD", "MeanNND"),
  
  Observed = c(
    disp_metrics$MinPD,
    disp_metrics$MeanPD,
    disp_metrics$MeanNND
  ),
  
  Random_Mean = c(
    mean(random_dist_metrics$MinPD),
    mean(random_dist_metrics$MeanPD),
    mean(random_dist_metrics$MeanNND)
  ),
  
  Random_SD = c(
    sd(random_dist_metrics$MinPD),
    sd(random_dist_metrics$MeanPD),
    sd(random_dist_metrics$MeanNND)
  ),
  
  Random_Q025 = c(
    quantile(random_dist_metrics$MinPD, 0.025),
    quantile(random_dist_metrics$MeanPD, 0.025),
    quantile(random_dist_metrics$MeanNND, 0.025)
  ),
  
  Random_Q975 = c(
    quantile(random_dist_metrics$MinPD, 0.975),
    quantile(random_dist_metrics$MeanPD, 0.975),
    quantile(random_dist_metrics$MeanNND, 0.975)
  ),
  
  SES = c(
    calc_ses(disp_metrics$MinPD, random_dist_metrics$MinPD),
    calc_ses(disp_metrics$MeanPD, random_dist_metrics$MeanPD),
    calc_ses(disp_metrics$MeanNND, random_dist_metrics$MeanNND)
  ),
  
  P_value = c(
    calc_p_high(disp_metrics$MinPD, random_dist_metrics$MinPD),
    calc_p_high(disp_metrics$MeanPD, random_dist_metrics$MeanPD),
    calc_p_high(disp_metrics$MeanNND, random_dist_metrics$MeanNND)
  )
)

disp_distance_summary

# ============================================================
# 13B. Clustered subset vs random baseline
# 指标越小越聚集：low-tail p-value
# NEW metrics: MeanPD, MeanNND, MaxPD
# ============================================================

clust_distance_summary <- data.frame(
  Subset_Type = "clustered",
  Metric = c("MeanPD", "MeanNND", "MaxPD"),
  
  Observed = c(
    clust_metrics$MeanPD,
    clust_metrics$MeanNND,
    clust_metrics$MaxPD
  ),
  
  Random_Mean = c(
    mean(random_dist_metrics$MeanPD),
    mean(random_dist_metrics$MeanNND),
    mean(random_dist_metrics$MaxPD)
  ),
  
  Random_SD = c(
    sd(random_dist_metrics$MeanPD),
    sd(random_dist_metrics$MeanNND),
    sd(random_dist_metrics$MaxPD)
  ),
  
  Random_Q025 = c(
    quantile(random_dist_metrics$MeanPD, 0.025),
    quantile(random_dist_metrics$MeanNND, 0.025),
    quantile(random_dist_metrics$MaxPD, 0.025)
  ),
  
  Random_Q975 = c(
    quantile(random_dist_metrics$MeanPD, 0.975),
    quantile(random_dist_metrics$MeanNND, 0.975),
    quantile(random_dist_metrics$MaxPD, 0.975)
  ),
  
  SES = c(
    calc_ses(clust_metrics$MeanPD, random_dist_metrics$MeanPD),
    calc_ses(clust_metrics$MeanNND, random_dist_metrics$MeanNND),
    calc_ses(clust_metrics$MaxPD, random_dist_metrics$MaxPD)
  ),
  
  P_value = c(
    calc_p_low(clust_metrics$MeanPD, random_dist_metrics$MeanPD),
    calc_p_low(clust_metrics$MeanNND, random_dist_metrics$MeanNND),
    calc_p_low(clust_metrics$MaxPD, random_dist_metrics$MaxPD)
  )
)

clust_distance_summary

# ============================================================
# 13C. 合并距离结果
# ============================================================

distance_summary_new <- rbind(
  disp_distance_summary,
  clust_distance_summary
)

distance_summary_new

write.csv(
  distance_summary_new,
  paste0("Carnivora_distance_summary_new_s", subset_size, ".csv"),
  row.names = FALSE
)
# ============================================================
# 14. 计算 Carnivora full candidate tree 的 BM covariance matrix
# ============================================================

V_carn <- vcv.phylo(carn_tree, corr = FALSE)

dim(V_carn)
all(rownames(V_carn) == carn_tree$tip.label)
all(colnames(V_carn) == carn_tree$tip.label)

# ============================================================
# 15. 定义 RS-based dependence metrics
# ============================================================

calc_dependence_from_V <- function(V_full, subset_names) {
  
  V_sub <- V_full[subset_names, subset_names, drop = FALSE]
  
  R_sub <- cov2cor(V_sub)
  
  off_values <- R_sub[upper.tri(R_sub)]
  
  n <- length(subset_names)
  one_vec <- rep(1, n)
  
  neff_mean <- as.numeric(
    t(one_vec) %*% solve(R_sub, one_vec)
  )
  
  data.frame(
    off_mean = mean(off_values),
    rmax = max(off_values),
    neff_mean = neff_mean
  )
}

# ============================================================
# 16. 计算 dispersed 和 clustered 的 dependence diagnostics
# ============================================================

disp_dep <- calc_dependence_from_V(V_carn, disp_subset_names)
clust_dep <- calc_dependence_from_V(V_carn, clust_subset_names)

disp_dep
clust_dep

# ============================================================
# 17. 计算 random subsets 的 dependence diagnostics
# ============================================================

random_subsets_names <- lapply(random_subsets_idx, function(x) {
  dist_obj$tip_labels[x]
})

random_dep_list <- vector("list", length(random_subsets_names))

for (i in seq_along(random_subsets_names)) {
  random_dep_list[[i]] <- calc_dependence_from_V(
    V_full = V_carn,
    subset_names = random_subsets_names[[i]]
  )
  
  if (i %% 1000 == 0) {
    cat("Finished", i, "random dependence diagnostics\n")
  }
}

random_dep_metrics <- do.call(rbind, random_dep_list)

random_dep_metrics$SubsetID <- seq_len(nrow(random_dep_metrics))

head(random_dep_metrics)
summary(random_dep_metrics)

write.csv(random_dep_metrics,
          "Carnivora_random_dependence_metrics_s20.csv",
          row.names = FALSE)

# ============================================================
# 18. Dependence diagnostics summary
# ============================================================

dependence_summary_new <- data.frame(
  Subset_Type = rep(c("dispersed", "clustered"), each = 3),
  
  Metric = rep(c("off_mean", "rmax", "neff_mean"), times = 2),
  
  Observed = c(
    disp_dep$off_mean,
    disp_dep$rmax,
    disp_dep$neff_mean,
    clust_dep$off_mean,
    clust_dep$rmax,
    clust_dep$neff_mean
  ),
  
  Random_Mean = c(
    mean(random_dep_metrics$off_mean),
    mean(random_dep_metrics$rmax),
    mean(random_dep_metrics$neff_mean),
    mean(random_dep_metrics$off_mean),
    mean(random_dep_metrics$rmax),
    mean(random_dep_metrics$neff_mean)
  ),
  
  Random_SD = c(
    sd(random_dep_metrics$off_mean),
    sd(random_dep_metrics$rmax),
    sd(random_dep_metrics$neff_mean),
    sd(random_dep_metrics$off_mean),
    sd(random_dep_metrics$rmax),
    sd(random_dep_metrics$neff_mean)
  ),
  
  Random_Q025 = c(
    quantile(random_dep_metrics$off_mean, 0.025),
    quantile(random_dep_metrics$rmax, 0.025),
    quantile(random_dep_metrics$neff_mean, 0.025),
    quantile(random_dep_metrics$off_mean, 0.025),
    quantile(random_dep_metrics$rmax, 0.025),
    quantile(random_dep_metrics$neff_mean, 0.025)
  ),
  
  Random_Q975 = c(
    quantile(random_dep_metrics$off_mean, 0.975),
    quantile(random_dep_metrics$rmax, 0.975),
    quantile(random_dep_metrics$neff_mean, 0.975),
    quantile(random_dep_metrics$off_mean, 0.975),
    quantile(random_dep_metrics$rmax, 0.975),
    quantile(random_dep_metrics$neff_mean, 0.975)
  ),
  
  SES = c(
    calc_ses(disp_dep$off_mean, random_dep_metrics$off_mean),
    calc_ses(disp_dep$rmax, random_dep_metrics$rmax),
    calc_ses(disp_dep$neff_mean, random_dep_metrics$neff_mean),
    calc_ses(clust_dep$off_mean, random_dep_metrics$off_mean),
    calc_ses(clust_dep$rmax, random_dep_metrics$rmax),
    calc_ses(clust_dep$neff_mean, random_dep_metrics$neff_mean)
  ),
  
  P_value = c(
    calc_p_low(disp_dep$off_mean, random_dep_metrics$off_mean),
    calc_p_low(disp_dep$rmax, random_dep_metrics$rmax),
    calc_p_high(disp_dep$neff_mean, random_dep_metrics$neff_mean),
    calc_p_high(clust_dep$off_mean, random_dep_metrics$off_mean),
    calc_p_high(clust_dep$rmax, random_dep_metrics$rmax),
    calc_p_low(clust_dep$neff_mean, random_dep_metrics$neff_mean)
  )
)

dependence_summary_new

write.csv(
  dependence_summary_new,
  paste0("Carnivora_dependence_summary_new_s", subset_size, ".csv"),
  row.names = FALSE
)
# ============================================================
# 19. 提取 selected species 的 taxonomy 信息
# ============================================================

tax_cols <- intersect(
  c("tiplabel", "Species_Name", "sciName", "family", "fam", "ord", "samp", "extinct."),
  names(tax)
)

tax_cols

disp_tax <- tax[match(disp_subset_names, tax$tiplabel), tax_cols]
clust_tax <- tax[match(clust_subset_names, tax$tiplabel), tax_cols]

disp_tax$Subset_Type <- "dispersed"
clust_tax$Subset_Type <- "clustered"

selected_tax <- rbind(disp_tax, clust_tax)

selected_tax

write.csv(selected_tax,
          "Carnivora_selected_species_taxonomy_s20.csv",
          row.names = FALSE)
# ============================================================
# 20. 画树函数（改进版：将点沿径向向外移动）
# ============================================================

plot_subset_on_carn_tree <- function(tree, subset_names, file_name, main_title,
                                     point_cex = 1.0,
                                     outward_offset_frac = 0.03) {
  
  pdf(file_name, width = 9, height = 9)
  
  par(mar = c(1, 1, 3, 1))
  
  # 先画树
  plot.phylo(
    tree,
    type = "fan",
    show.tip.label = FALSE,
    no.margin = FALSE,
    cex = 0.3
  )
  
  # 读取 ape 存储的绘图坐标
  pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  n_tip <- Ntip(tree)
  tip_idx <- which(tree$tip.label %in% subset_names)
  
  # tip 坐标
  tip_x <- pp$xx[tip_idx]
  tip_y <- pp$yy[tip_idx]
  
  # 根节点坐标（通常是第一个内部节点，即 n_tip + 1）
  root_x <- pp$xx[n_tip + 1]
  root_y <- pp$yy[n_tip + 1]
  
  # 从根指向各 tip 的方向向量
  dx <- tip_x - root_x
  dy <- tip_y - root_y
  
  r <- sqrt(dx^2 + dy^2)
  
  # 防止除零
  r[r == 0] <- 1
  
  # outward offset：按整棵树尺度的一定比例外移
  # 这里用 tip 到根的最大半径来确定平移量
  all_tip_x <- pp$xx[1:n_tip]
  all_tip_y <- pp$yy[1:n_tip]
  all_r <- sqrt((all_tip_x - root_x)^2 + (all_tip_y - root_y)^2)
  offset <- max(all_r) * outward_offset_frac
  
  # 平移后的点坐标
  tip_x_out <- (tip_x + offset * dx / r )*0.985
  tip_y_out <- (tip_y + offset * dy / r )*0.985
  
  # 画点
  points(
    tip_x_out,
    tip_y_out,
    pch = 21,
    bg = "black",
    cex = 1.1,
    xpd = NA
  )
  
  title(main_title, line = 1)
  
  dev.off()
}
# ============================================================
# 21. 树上标出 dispersed subset
# ============================================================

plot_subset_on_carn_tree(
  tree = carn_tree,
  subset_names = disp_subset_names,
  file_name = paste0("Carnivora_tree_dispersed_subset_new_s", subset_size, ".pdf"),
  main_title = paste0("Phylogenetically dispersed subset, s = ", subset_size)
)
# ============================================================
# 22. 树上标出 clustered subset
# ============================================================

plot_subset_on_carn_tree(
  tree = carn_tree,
  subset_names = clust_subset_names,
  file_name = paste0("Carnivora_tree_clustered_subset_new_s", subset_size, ".pdf"),
  main_title = paste0("Phylogenetically clustered subset, s = ", subset_size)
)
# ============================================================
# 23. Distance metric random distributions
# ============================================================

pdf(
  paste0("Carnivora_random_distance_distributions_new_s", subset_size, ".pdf"),
  width = 10,
  height = 8
)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

hist(random_dist_metrics$MinPD,
     breaks = 40,
     main = "Random baseline: MinPD",
     xlab = "MinPD")
abline(v = disp_metrics$MinPD, lwd = 2, lty = 2)

hist(random_dist_metrics$MeanPD,
     breaks = 40,
     main = "Random baseline: MeanPD",
     xlab = "MeanPD")
abline(v = disp_metrics$MeanPD, lwd = 2, lty = 2)
abline(v = clust_metrics$MeanPD, lwd = 2, lty = 3)

hist(random_dist_metrics$MeanNND,
     breaks = 40,
     main = "Random baseline: MeanNND",
     xlab = "MeanNND")
abline(v = disp_metrics$MeanNND, lwd = 2, lty = 2)
abline(v = clust_metrics$MeanNND, lwd = 2, lty = 3)

hist(random_dist_metrics$MaxPD,
     breaks = 40,
     main = "Random baseline: MaxPD",
     xlab = "MaxPD")
abline(v = clust_metrics$MaxPD, lwd = 2, lty = 3)

plot.new()
legend(
  "center",
  legend = c("Dispersed", "Clustered"),
  lty = c(2, 3),
  lwd = 2,
  bty = "n"
)

dev.off()
# ============================================================
# 24. Dependence metric random distributions
# ============================================================

pdf(
  paste0("Carnivora_random_dependence_distributions_new_s", subset_size, ".pdf"),
  width = 10,
  height = 8
)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

hist(random_dep_metrics$off_mean,
     breaks = 40,
     main = "Random baseline: off_mean",
     xlab = "Mean off-diagonal correlation")
abline(v = disp_dep$off_mean, lwd = 2, lty = 2)
abline(v = clust_dep$off_mean, lwd = 2, lty = 3)

hist(random_dep_metrics$rmax,
     breaks = 40,
     main = "Random baseline: rmax",
     xlab = "Maximum off-diagonal correlation")
abline(v = disp_dep$rmax, lwd = 2, lty = 2)
abline(v = clust_dep$rmax, lwd = 2, lty = 3)

hist(random_dep_metrics$neff_mean,
     breaks = 40,
     main = "Random baseline: neff_mean",
     xlab = "Mean-based effective sample size")
abline(v = disp_dep$neff_mean, lwd = 2, lty = 2)
abline(v = clust_dep$neff_mean, lwd = 2, lty = 3)

plot.new()
legend(
  "center",
  legend = c("Dispersed", "Clustered"),
  lty = c(2, 3),
  lwd = 2,
  bty = "n"
)

dev.off()

# ============================================================
# 25. Save random baseline raw distributions and summaries
# 保存随机基线的原始分布和统计汇总
# ============================================================

# ---------- 1. 确认对象存在 ----------
stopifnot(exists("random_dist_metrics"))
stopifnot(exists("random_dep_metrics"))
stopifnot(exists("disp_metrics"))
stopifnot(exists("clust_metrics"))
stopifnot(exists("disp_dep"))
stopifnot(exists("clust_dep"))

# ---------- 2. 确认 p-value 和 SES 函数存在 ----------
calc_p_high <- function(obs, null_values) {
  (1 + sum(null_values >= obs)) / (1 + length(null_values))
}

calc_p_low <- function(obs, null_values) {
  (1 + sum(null_values <= obs)) / (1 + length(null_values))
}

calc_ses <- function(obs, null_values) {
  (obs - mean(null_values)) / sd(null_values)
}

# ============================================================
# 25.1 保存 random baseline 的原始分布
# ============================================================

write.csv(
  random_dist_metrics,
  paste0("Carnivora_random_distance_baseline_raw_s", subset_size, ".csv"),
  row.names = FALSE
)

write.csv(
  random_dep_metrics,
  paste0("Carnivora_random_dependence_baseline_raw_s", subset_size, ".csv"),
  row.names = FALSE
)
# ============================================================
# 25.2 距离指标 random baseline summary
# ============================================================

summarize_null_metric <- function(x) {
  data.frame(
    N_random = length(x),
    Mean = mean(x),
    SD = sd(x),
    Median = median(x),
    Min = min(x),
    Q025 = as.numeric(quantile(x, 0.025)),
    Q25 = as.numeric(quantile(x, 0.25)),
    Q75 = as.numeric(quantile(x, 0.75)),
    Q975 = as.numeric(quantile(x, 0.975)),
    Max = max(x)
  )
}

distance_baseline_summary <- rbind(
  data.frame(
    Metric = "MinPD",
    summarize_null_metric(random_dist_metrics$MinPD)
  ),
  data.frame(
    Metric = "MeanPD",
    summarize_null_metric(random_dist_metrics$MeanPD)
  ),
  data.frame(
    Metric = "MeanNND",
    summarize_null_metric(random_dist_metrics$MeanNND)
  ),
  data.frame(
    Metric = "MaxPD",
    summarize_null_metric(random_dist_metrics$MaxPD)
  )
)

distance_baseline_summary

write.csv(
  distance_baseline_summary,
  paste0("Carnivora_random_distance_baseline_summary_s", subset_size, ".csv"),
  row.names = FALSE
)
# ============================================================
# 25.3 依赖诊断 random baseline summary
# ============================================================

dependence_baseline_summary <- rbind(
  data.frame(
    Metric = "off_mean",
    summarize_null_metric(random_dep_metrics$off_mean)
  ),
  data.frame(
    Metric = "rmax",
    summarize_null_metric(random_dep_metrics$rmax)
  ),
  data.frame(
    Metric = "neff_mean",
    summarize_null_metric(random_dep_metrics$neff_mean)
  )
)

dependence_baseline_summary

write.csv(
  dependence_baseline_summary,
  paste0("Carnivora_random_dependence_baseline_summary_s", subset_size, ".csv"),
  row.names = FALSE
)
# ============================================================
# 25.4 距离指标：observed subsets vs random baseline
# ============================================================

distance_observed_vs_random <- rbind(
  
  # Dispersed subset: larger values are more extreme
  data.frame(
    Subset_Type = "dispersed",
    Metric = "MinPD",
    Observed = disp_metrics$MinPD,
    Baseline_Mean = mean(random_dist_metrics$MinPD),
    Baseline_SD = sd(random_dist_metrics$MinPD),
    Baseline_Q025 = as.numeric(quantile(random_dist_metrics$MinPD, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dist_metrics$MinPD, 0.975)),
    SES = calc_ses(disp_metrics$MinPD, random_dist_metrics$MinPD),
    P_value = calc_p_high(disp_metrics$MinPD, random_dist_metrics$MinPD),
    Tail = "high"
  ),
  
  data.frame(
    Subset_Type = "dispersed",
    Metric = "MeanPD",
    Observed = disp_metrics$MeanPD,
    Baseline_Mean = mean(random_dist_metrics$MeanPD),
    Baseline_SD = sd(random_dist_metrics$MeanPD),
    Baseline_Q025 = as.numeric(quantile(random_dist_metrics$MeanPD, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dist_metrics$MeanPD, 0.975)),
    SES = calc_ses(disp_metrics$MeanPD, random_dist_metrics$MeanPD),
    P_value = calc_p_high(disp_metrics$MeanPD, random_dist_metrics$MeanPD),
    Tail = "high"
  ),
  
  data.frame(
    Subset_Type = "dispersed",
    Metric = "MeanNND",
    Observed = disp_metrics$MeanNND,
    Baseline_Mean = mean(random_dist_metrics$MeanNND),
    Baseline_SD = sd(random_dist_metrics$MeanNND),
    Baseline_Q025 = as.numeric(quantile(random_dist_metrics$MeanNND, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dist_metrics$MeanNND, 0.975)),
    SES = calc_ses(disp_metrics$MeanNND, random_dist_metrics$MeanNND),
    P_value = calc_p_high(disp_metrics$MeanNND, random_dist_metrics$MeanNND),
    Tail = "high"
  ),
  
  # Clustered subset: smaller values are more extreme
  data.frame(
    Subset_Type = "clustered",
    Metric = "MeanPD",
    Observed = clust_metrics$MeanPD,
    Baseline_Mean = mean(random_dist_metrics$MeanPD),
    Baseline_SD = sd(random_dist_metrics$MeanPD),
    Baseline_Q025 = as.numeric(quantile(random_dist_metrics$MeanPD, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dist_metrics$MeanPD, 0.975)),
    SES = calc_ses(clust_metrics$MeanPD, random_dist_metrics$MeanPD),
    P_value = calc_p_low(clust_metrics$MeanPD, random_dist_metrics$MeanPD),
    Tail = "low"
  ),
  
  data.frame(
    Subset_Type = "clustered",
    Metric = "MeanNND",
    Observed = clust_metrics$MeanNND,
    Baseline_Mean = mean(random_dist_metrics$MeanNND),
    Baseline_SD = sd(random_dist_metrics$MeanNND),
    Baseline_Q025 = as.numeric(quantile(random_dist_metrics$MeanNND, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dist_metrics$MeanNND, 0.975)),
    SES = calc_ses(clust_metrics$MeanNND, random_dist_metrics$MeanNND),
    P_value = calc_p_low(clust_metrics$MeanNND, random_dist_metrics$MeanNND),
    Tail = "low"
  ),
  
  data.frame(
    Subset_Type = "clustered",
    Metric = "MaxPD",
    Observed = clust_metrics$MaxPD,
    Baseline_Mean = mean(random_dist_metrics$MaxPD),
    Baseline_SD = sd(random_dist_metrics$MaxPD),
    Baseline_Q025 = as.numeric(quantile(random_dist_metrics$MaxPD, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dist_metrics$MaxPD, 0.975)),
    SES = calc_ses(clust_metrics$MaxPD, random_dist_metrics$MaxPD),
    P_value = calc_p_low(clust_metrics$MaxPD, random_dist_metrics$MaxPD),
    Tail = "low"
  )
)

distance_observed_vs_random

write.csv(
  distance_observed_vs_random,
  paste0("Carnivora_distance_observed_vs_random_s", subset_size, ".csv"),
  row.names = FALSE
)

# ============================================================
# 25.5 依赖诊断：observed subsets vs random baseline
# ============================================================

dependence_observed_vs_random <- rbind(
  
  # Dispersed subset: off_mean and rmax lower are favorable; neff higher is favorable
  data.frame(
    Subset_Type = "dispersed",
    Metric = "off_mean",
    Observed = disp_dep$off_mean,
    Baseline_Mean = mean(random_dep_metrics$off_mean),
    Baseline_SD = sd(random_dep_metrics$off_mean),
    Baseline_Q025 = as.numeric(quantile(random_dep_metrics$off_mean, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dep_metrics$off_mean, 0.975)),
    SES = calc_ses(disp_dep$off_mean, random_dep_metrics$off_mean),
    P_value = calc_p_low(disp_dep$off_mean, random_dep_metrics$off_mean),
    Tail = "low"
  ),
  
  data.frame(
    Subset_Type = "dispersed",
    Metric = "rmax",
    Observed = disp_dep$rmax,
    Baseline_Mean = mean(random_dep_metrics$rmax),
    Baseline_SD = sd(random_dep_metrics$rmax),
    Baseline_Q025 = as.numeric(quantile(random_dep_metrics$rmax, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dep_metrics$rmax, 0.975)),
    SES = calc_ses(disp_dep$rmax, random_dep_metrics$rmax),
    P_value = calc_p_low(disp_dep$rmax, random_dep_metrics$rmax),
    Tail = "low"
  ),
  
  data.frame(
    Subset_Type = "dispersed",
    Metric = "neff_mean",
    Observed = disp_dep$neff_mean,
    Baseline_Mean = mean(random_dep_metrics$neff_mean),
    Baseline_SD = sd(random_dep_metrics$neff_mean),
    Baseline_Q025 = as.numeric(quantile(random_dep_metrics$neff_mean, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dep_metrics$neff_mean, 0.975)),
    SES = calc_ses(disp_dep$neff_mean, random_dep_metrics$neff_mean),
    P_value = calc_p_high(disp_dep$neff_mean, random_dep_metrics$neff_mean),
    Tail = "high"
  ),
  
  # Clustered subset: off_mean and rmax higher are favorable; neff lower is favorable
  data.frame(
    Subset_Type = "clustered",
    Metric = "off_mean",
    Observed = clust_dep$off_mean,
    Baseline_Mean = mean(random_dep_metrics$off_mean),
    Baseline_SD = sd(random_dep_metrics$off_mean),
    Baseline_Q025 = as.numeric(quantile(random_dep_metrics$off_mean, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dep_metrics$off_mean, 0.975)),
    SES = calc_ses(clust_dep$off_mean, random_dep_metrics$off_mean),
    P_value = calc_p_high(clust_dep$off_mean, random_dep_metrics$off_mean),
    Tail = "high"
  ),
  
  data.frame(
    Subset_Type = "clustered",
    Metric = "rmax",
    Observed = clust_dep$rmax,
    Baseline_Mean = mean(random_dep_metrics$rmax),
    Baseline_SD = sd(random_dep_metrics$rmax),
    Baseline_Q025 = as.numeric(quantile(random_dep_metrics$rmax, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dep_metrics$rmax, 0.975)),
    SES = calc_ses(clust_dep$rmax, random_dep_metrics$rmax),
    P_value = calc_p_high(clust_dep$rmax, random_dep_metrics$rmax),
    Tail = "high"
  ),
  
  data.frame(
    Subset_Type = "clustered",
    Metric = "neff_mean",
    Observed = clust_dep$neff_mean,
    Baseline_Mean = mean(random_dep_metrics$neff_mean),
    Baseline_SD = sd(random_dep_metrics$neff_mean),
    Baseline_Q025 = as.numeric(quantile(random_dep_metrics$neff_mean, 0.025)),
    Baseline_Q975 = as.numeric(quantile(random_dep_metrics$neff_mean, 0.975)),
    SES = calc_ses(clust_dep$neff_mean, random_dep_metrics$neff_mean),
    P_value = calc_p_low(clust_dep$neff_mean, random_dep_metrics$neff_mean),
    Tail = "low"
  )
)

dependence_observed_vs_random

write.csv(
  dependence_observed_vs_random,
  paste0("Carnivora_dependence_observed_vs_random_s", subset_size, ".csv"),
  row.names = FALSE
)

# ============================================================
# 25.6 合并全部 observed-vs-random 结果
# ============================================================

all_observed_vs_random <- rbind(
  data.frame(
    Analysis = "distance",
    distance_observed_vs_random,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "dependence",
    dependence_observed_vs_random,
    stringsAsFactors = FALSE
  )
)

all_observed_vs_random

write.csv(
  all_observed_vs_random,
  paste0("Carnivora_all_observed_vs_random_summary_s", subset_size, ".csv"),
  row.names = FALSE
)
# ============================================================
# 26. Idealized tree sanity checks
# 理想化系统发育树上的算法合理性检查
# ============================================================

# 本部分假设你前面已经 source 了：
# distance_metrics.R
# objective_compare.R
# subset_exchange.R
# 并且已经定义了新版函数：
# run_dispersed_fixed_start_algorithm()
# select_best_clustered_neighborhood()

# 检查必要函数是否存在
stopifnot(exists("create_distance_object"))
stopifnot(exists("calc_subset_metrics"))
stopifnot(exists("refine_subset_exchange"))
stopifnot(exists("run_dispersed_fixed_start_algorithm"))
stopifnot(exists("select_best_clustered_neighborhood"))

# ============================================================
# 26.1 Helper functions for idealized-tree checks
# ============================================================

# ---------- 1. MaxPD ----------
if (!exists("calc_maxpd")) {
  calc_maxpd <- function(dist_mat, subset) {
    if (length(subset) < 2) return(0)
    sub_dist <- dist_mat[subset, subset]
    max(sub_dist[upper.tri(sub_dist)], na.rm = TRUE)
  }
}

# ---------- 2. 扩展指标 ----------
calc_subset_metrics_extended_local <- function(dist_mat, subset) {
  m <- calc_subset_metrics(dist_mat, subset)
  data.frame(
    MinPD = m$MinPD,
    MeanPD = m$MeanPD,
    MeanNND = m$MeanNND,
    MaxPD = calc_maxpd(dist_mat, subset)
  )
}

# ---------- 3. 获取某个 node 的所有 descendant tips ----------
get_descendant_tips <- function(tree, node) {
  
  n_tip <- Ntip(tree)
  
  children <- tree$edge[tree$edge[, 1] == node, 2]
  tips <- children[children <= n_tip]
  internal_nodes <- children[children > n_tip]
  
  while (length(internal_nodes) > 0) {
    new_children <- unlist(lapply(internal_nodes, function(nd) {
      tree$edge[tree$edge[, 1] == nd, 2]
    }))
    
    tips <- c(tips, new_children[new_children <= n_tip])
    internal_nodes <- new_children[new_children > n_tip]
  }
  
  tree$tip.label[sort(tips)]
}

# ---------- 4. 从 root 向下走指定 edge-depth，得到这一层的节点 ----------
get_nodes_at_depth_from_root <- function(tree, depth) {
  
  current_nodes <- Ntip(tree) + 1
  
  if (depth == 0) return(current_nodes)
  
  for (i in seq_len(depth)) {
    current_nodes <- unlist(lapply(current_nodes, function(nd) {
      tree$edge[tree$edge[, 1] == nd, 2]
    }))
  }
  
  current_nodes
}

# ---------- 5. 给一组 clade nodes 建立 tip -> clade 的对应表 ----------
make_clade_membership <- function(tree, nodes, group_prefix = "clade") {
  
  out <- vector("list", length(nodes))
  
  for (i in seq_along(nodes)) {
    tips_i <- get_descendant_tips(tree, nodes[i])
    
    out[[i]] <- data.frame(
      tiplabel = tips_i,
      group = paste0(group_prefix, "_", sprintf("%02d", i)),
      node = nodes[i],
      clade_size = length(tips_i),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

# ---------- 6. 找到所有刚好包含 k 个 tips 的 internal nodes ----------
find_nodes_with_k_tips <- function(tree, k) {
  
  internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)
  
  sizes <- sapply(internal_nodes, function(nd) {
    length(get_descendant_tips(tree, nd))
  })
  
  internal_nodes[sizes == k]
}

# ---------- 7. 画图：所有树枝都用细黑线；被选中的 tips 用加粗黑点 ----------
plot_selected_tips_on_tree <- function(tree, selected_tips, main_title) {
  
  selected_idx <- match(selected_tips, tree$tip.label)
  
  plot.phylo(
    tree,
    type = "phylogram",
    show.tip.label = FALSE,
    edge.color = "black",
    edge.width = 0.8,
    no.margin = FALSE,
    cex = 0.4
  )
  
  tiplabels(
    pch = 16,
    col = "black",
    cex = 1.2,
    tip = selected_idx
  )
  
  title(main_title, line = 1)
}

# ---------- 8. 构造严格超度量的 ladder-like / pectinate tree ----------
# 通过手工构造 hclust 对象，再用 as.phylo 转为 rooted ultrametric tree
make_ultrametric_ladder_tree <- function(n_tips) {
  
  if (n_tips < 2) stop("n_tips must be >= 2")
  
  merge_mat <- matrix(NA_integer_, nrow = n_tips - 1, ncol = 2)
  heights <- numeric(n_tips - 1)
  
  # 第一次合并 tip1 与 tip2
  merge_mat[1, ] <- c(-1, -2)
  heights[1] <- 1
  
  # 之后每次把前一个 cluster 与下一个新 tip 合并，构成 pectinate / ladder tree
  if (n_tips >= 3) {
    for (i in 2:(n_tips - 1)) {
      merge_mat[i, ] <- c(i - 1, -(i + 1))
      heights[i] <- i
    }
  }
  
  hc <- list(
    merge = merge_mat,
    height = heights,
    order = 1:n_tips,
    labels = paste0("L", sprintf("%03d", 1:n_tips)),
    method = "complete",
    call = match.call()
  )
  
  class(hc) <- "hclust"
  
  tr <- as.phylo(hc)
  
  # 缩放到 root-to-tip = 1，仍保持超度量
  tr$edge.length <- tr$edge.length / max(node.depth.edgelength(tr)[1:Ntip(tr)])
  
  tr
}

# ============================================================
# 26.2 Generate idealized ultrametric phylogenetic trees
# ============================================================

sanity_n_tips <- 32
sanity_s <- 8

# ---------- 严格平衡超度量树 ----------
balanced_tree <- stree(sanity_n_tips, type = "balanced")
balanced_tree$tip.label <- paste0("B", sprintf("%03d", seq_len(sanity_n_tips)))

# 为 balanced tree 赋 branch lengths，使之成为严格 ultrametric
# 这里所有 edge 长度先设为 1，再整体缩放；对 balanced stree 这样做仍是 ultrametric
balanced_tree$edge.length <- rep(1, nrow(balanced_tree$edge))
balanced_tree$edge.length <- balanced_tree$edge.length /
  max(node.depth.edgelength(balanced_tree)[1:Ntip(balanced_tree)])

# ---------- 阶梯状（ladder-like / pectinate）严格超度量树 ----------
ladder_tree <- make_ultrametric_ladder_tree(sanity_n_tips)

# ---------- 基本检查 ----------
Ntip(balanced_tree)
Ntip(ladder_tree)

is.rooted(balanced_tree)
is.rooted(ladder_tree)

is.ultrametric(balanced_tree, tol = 1e-8)
is.ultrametric(ladder_tree, tol = 1e-8)

# 建立 distance object
balanced_dist_obj <- create_distance_object(balanced_tree)
ladder_dist_obj <- create_distance_object(ladder_tree)

dim(balanced_dist_obj$dist_mat)
dim(ladder_dist_obj$dist_mat)

# ============================================================
# 26.3 Experiment 1:
# Dispersed subset on 128-tip balanced tree
# Expected: one selected tip from each of eight 16-tip clades
# ============================================================

balanced_disp_result <- run_dispersed_fixed_start_algorithm(
  dist_obj = balanced_dist_obj,
  subset_size = sanity_s
)

balanced_disp_idx <- balanced_disp_result$final_subset
balanced_disp_names <- balanced_disp_result$final_subset_names
balanced_disp_metrics <- calc_subset_metrics_extended_local(
  balanced_dist_obj$dist_mat,
  balanced_disp_idx
)

balanced_disp_names
balanced_disp_metrics

# 检查是否每个 16-tip clade 选 1 个
# root 下第 3 层对应 8 个 major 16-tip clades
major16_nodes <- get_nodes_at_depth_from_root(balanced_tree, depth = 3)

length(major16_nodes)

major16_membership <- make_clade_membership(
  tree = balanced_tree,
  nodes = major16_nodes,
  group_prefix = "major16"
)

table(major16_membership$group)

balanced_disp_major16_groups <- major16_membership$group[
  match(balanced_disp_names, major16_membership$tiplabel)
]

balanced_disp_major16_counts <- table(balanced_disp_major16_groups)

balanced_disp_major16_counts

balanced_disp_pass <- (
  length(balanced_disp_major16_counts) == 8 &&
    all(as.integer(balanced_disp_major16_counts) == 1)
)

balanced_disp_pass

# ============================================================
# 26.4 Experiment 2:
# Clustered subset on 128-tip balanced tree
# Expected: all eight selected tips from one 8-tip clade
# ============================================================

balanced_clust_result <- select_best_clustered_neighborhood(
  dist_obj = balanced_dist_obj,
  subset_size = sanity_s
)

balanced_clust_idx <- balanced_clust_result$final_subset
balanced_clust_names <- balanced_clust_result$final_subset_names
balanced_clust_metrics <- calc_subset_metrics_extended_local(
  balanced_dist_obj$dist_mat,
  balanced_clust_idx
)

balanced_clust_names
balanced_clust_metrics
balanced_clust_result$best_seed_name
balanced_clust_result$n_raw_candidates
balanced_clust_result$n_unique_candidates

# 检查是否属于同一个 8-tip clade
# root 下第 4 层对应 16 个 8-tip clades
local8_nodes_balanced <- get_nodes_at_depth_from_root(balanced_tree, depth = 4)

length(local8_nodes_balanced)

local8_membership_balanced <- make_clade_membership(
  tree = balanced_tree,
  nodes = local8_nodes_balanced,
  group_prefix = "local8"
)

table(local8_membership_balanced$group)

balanced_clust_local8_groups <- local8_membership_balanced$group[
  match(balanced_clust_names, local8_membership_balanced$tiplabel)
]

balanced_clust_local8_counts <- table(balanced_clust_local8_groups)

balanced_clust_local8_counts

balanced_clust_pass <- (
  length(balanced_clust_local8_counts) == 1 &&
    as.integer(balanced_clust_local8_counts[1]) == sanity_s
)

balanced_clust_pass

# ============================================================
# 26.5 Experiment 3:
# Clustered subset on 128-tip ladder-like ultrametric tree
# Expected: eight most recently diverged terminal tips
# ============================================================

ladder_clust_result <- select_best_clustered_neighborhood(
  dist_obj = ladder_dist_obj,
  subset_size = sanity_s
)

ladder_clust_idx <- ladder_clust_result$final_subset
ladder_clust_names <- ladder_clust_result$final_subset_names
ladder_clust_metrics <- calc_subset_metrics_extended_local(
  ladder_dist_obj$dist_mat,
  ladder_clust_idx
)

ladder_clust_names
ladder_clust_metrics
ladder_clust_result$best_seed_name
ladder_clust_result$n_raw_candidates
ladder_clust_result$n_unique_candidates

# 找出阶梯状树中“最近分化的 8-tip clade”
# 定义：在所有包含 8 个 descendants 的 internal nodes 中，
# 距离 root 最远的那个 node，即最晚形成的 8-tip group
ladder_8tip_nodes <- find_nodes_with_k_tips(ladder_tree, k = sanity_s)

ladder_8tip_nodes

ladder_node_depths <- node.depth.edgelength(ladder_tree)

most_recent_8tip_node <- ladder_8tip_nodes[
  which.max(ladder_node_depths[ladder_8tip_nodes])
]

most_recent_8tip_node

expected_ladder_recent8_names <- get_descendant_tips(
  ladder_tree,
  most_recent_8tip_node
)

expected_ladder_recent8_names

# 检查算法结果是否等于最近分化的 8 个 tips
ladder_clust_pass <- setequal(
  ladder_clust_names,
  expected_ladder_recent8_names
)

ladder_clust_pass

# 查看差异；如果 character(0)，说明完全一致
setdiff(ladder_clust_names, expected_ladder_recent8_names)
setdiff(expected_ladder_recent8_names, ladder_clust_names)

# ============================================================
# 26.6 Summary tables
# ============================================================

sanity_check_summary <- data.frame(
  Experiment = c(
    "Balanced tree: dispersed subset",
    "Balanced tree: clustered subset",
    "Ladder tree: clustered subset"
  ),
  
  Tree_Type = c(
    "balanced_128",
    "balanced_128",
    "ladder_128_ultrametric"
  ),
  
  Subset_Type = c(
    "dispersed",
    "clustered",
    "clustered"
  ),
  
  Expected_Pattern = c(
    "one selected species from each of eight 16-tip clades",
    "all selected species from one 8-tip clade",
    "selected species match the most recently diverged 8-tip clade"
  ),
  
  Pass = c(
    balanced_disp_pass,
    balanced_clust_pass,
    ladder_clust_pass
  ),
  
  stringsAsFactors = FALSE
)

sanity_check_summary

write.csv(
  sanity_check_summary,
  "Sanity_check_summary_idealized_trees_s8.csv",
  row.names = FALSE
)

# 距离指标汇总
sanity_metrics_summary <- rbind(
  data.frame(
    Experiment = "Balanced tree: dispersed subset",
    Subset_Type = "dispersed",
    balanced_disp_metrics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Experiment = "Balanced tree: clustered subset",
    Subset_Type = "clustered",
    balanced_clust_metrics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Experiment = "Ladder tree: clustered subset",
    Subset_Type = "clustered",
    ladder_clust_metrics,
    stringsAsFactors = FALSE
  )
)

sanity_metrics_summary

write.csv(
  sanity_metrics_summary,
  "Sanity_check_distance_metrics_idealized_trees_s8.csv",
  row.names = FALSE
)

# 被选 tip 列表
sanity_selected_tips <- rbind(
  data.frame(
    Experiment = "Balanced tree: dispersed subset",
    Subset_Type = "dispersed",
    Rank = seq_along(balanced_disp_names),
    tiplabel = balanced_disp_names,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Experiment = "Balanced tree: clustered subset",
    Subset_Type = "clustered",
    Rank = seq_along(balanced_clust_names),
    tiplabel = balanced_clust_names,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Experiment = "Ladder tree: clustered subset",
    Subset_Type = "clustered",
    Rank = seq_along(ladder_clust_names),
    tiplabel = ladder_clust_names,
    stringsAsFactors = FALSE
  )
)

sanity_selected_tips

write.csv(
  sanity_selected_tips,
  "Sanity_check_selected_tips_idealized_trees_s8.csv",
  row.names = FALSE
)

# ============================================================
# 26.7 Figure 1 panels as separate pdf files
# 所有树枝均为细黑线；选中 tips 为较大的黑点
# ============================================================

pdf("Fig1a_balanced_tree_dispersed_subset_s8.pdf",
    width = 7, height = 7)

par(mar = c(1, 1, 3, 1))

plot_selected_tips_on_tree(
  tree = balanced_tree,
  selected_tips = balanced_disp_names,
  main_title = "a) Balanced tree: dispersed subset"
)

dev.off()

pdf("Fig1b_balanced_tree_clustered_subset_s8.pdf",
    width = 7, height = 7)

par(mar = c(1, 1, 3, 1))

plot_selected_tips_on_tree(
  tree = balanced_tree,
  selected_tips = balanced_clust_names,
  main_title = "b) Balanced tree: clustered subset"
)

dev.off()

pdf("Fig1c_ladder_tree_clustered_subset_s8.pdf",
    width = 7, height = 7)

par(mar = c(1, 1, 3, 1))

plot_selected_tips_on_tree(
  tree = ladder_tree,
  selected_tips = ladder_clust_names,
  main_title = "c) Ladder tree: clustered subset"
)

dev.off()

# ============================================================
# 26.8 Figure 1 combined PDF
# ============================================================

pdf("Figure1_idealized_tree_sanity_checks_s8.pdf",
    width = 15, height = 5)

par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))

plot_selected_tips_on_tree(
  tree = balanced_tree,
  selected_tips = balanced_disp_names,
  main_title = "a) Balanced: dispersed"
)

plot_selected_tips_on_tree(
  tree = balanced_tree,
  selected_tips = balanced_clust_names,
  main_title = "b) Balanced: clustered"
)

plot_selected_tips_on_tree(
  tree = ladder_tree,
  selected_tips = ladder_clust_names,
  main_title = "c) Ladder: clustered"
)

dev.off()

# ============================================================
# 26.9 Console-readable checks
# ============================================================

cat("\n============================================================\n")
cat("Idealized tree sanity checks\n")
cat("============================================================\n\n")

cat("Ultrametric checks:\n")
cat("Balanced tree ultrametric:", is.ultrametric(balanced_tree, tol = 1e-8), "\n")
cat("Ladder tree ultrametric:", is.ultrametric(ladder_tree, tol = 1e-8), "\n\n")

cat("Experiment 1: Balanced tree dispersed subset\n")
cat("Expected: one selected tip from each of eight 16-tip clades\n")
print(balanced_disp_major16_counts)
cat("Pass:", balanced_disp_pass, "\n\n")

cat("Experiment 2: Balanced tree clustered subset\n")
cat("Expected: all selected tips from one 8-tip clade\n")
print(balanced_clust_local8_counts)
cat("Pass:", balanced_clust_pass, "\n\n")

cat("Experiment 3: Ladder tree clustered subset\n")
cat("Expected: selected tips match most recently diverged 8-tip clade\n")
cat("Expected tips:\n")
print(expected_ladder_recent8_names)
cat("Selected tips:\n")
print(ladder_clust_names)
cat("Pass:", ladder_clust_pass, "\n\n")

cat("Distance metrics:\n")
print(sanity_metrics_summary)


# ============================================================
# 27. Lambda-transformed BM sensitivity analysis
# 固定 Carnivora dispersed 20-species subset，只改变 lambda
# ============================================================

stopifnot(exists("V_carn"))
stopifnot(exists("disp_subset_names"))
stopifnot(exists("subset_size"))

# 需要测试的 lambda 值
lambda_values <- c(1.00, 0.75, 0.50, 0.25, 0.10, 0)

# ------------------------------------------------------------
# 27.1 Pagel's lambda covariance transformation
# 对 BM covariance matrix 的非对角元素乘以 lambda
# 对角线保持不变
# ------------------------------------------------------------

lambda_transform_cov <- function(V_bm, lambda) {
  
  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("lambda must be a single numeric value.")
  }
  
  if (lambda < 0 || lambda > 1) {
    stop("lambda should be between 0 and 1 for this sensitivity analysis.")
  }
  
  V_lam <- V_bm
  
  offdiag <- row(V_lam) != col(V_lam)
  V_lam[offdiag] <- lambda * V_lam[offdiag]
  
  V_lam
}

# ------------------------------------------------------------
# 27.2 从 correlation matrix 计算 dependence diagnostics
# 注意：
# neff_mean_current_code = 你当前代码使用的 n^2 / (1' R 1)
# neff_mean_gls_inverse = 论文中 1' R^{-1} 1 对应的形式
# ------------------------------------------------------------

calc_dependence_from_R <- function(R_sub) {
  
  if (!isSymmetric(R_sub, tol = 1e-8)) {
    warning("R_sub is not exactly symmetric; check numerical precision.")
  }
  
  off_values <- R_sub[upper.tri(R_sub)]
  
  n <- nrow(R_sub)
  one_vec <- rep(1, n)
  
  neff_mean <- tryCatch(
    {
      as.numeric(t(one_vec) %*% solve(R_sub, one_vec))
    },
    error = function(e) {
      warning("R_sub could not be inverted.")
      NA_real_
    }
  )
  
  data.frame(
    off_mean = mean(off_values),
    rmax = max(off_values),
    neff_mean = neff_mean
  )
}

# ------------------------------------------------------------
# 27.3 对固定 dispersed subset 在不同 lambda 下计算 diagnostics
# ------------------------------------------------------------

calc_lambda_diagnostics_for_subset <- function(V_full, subset_names, lambda) {
  
  V_sub_bm <- V_full[subset_names, subset_names, drop = FALSE]
  
  V_sub_lam <- lambda_transform_cov(
    V_bm = V_sub_bm,
    lambda = lambda
  )
  
  R_sub_lam <- cov2cor(V_sub_lam)
  
  diag_out <- calc_dependence_from_R(R_sub_lam)
  
  data.frame(
    lambda = lambda,
    n_species = length(subset_names),
    diag_out
  )
}

disp_lambda_results <- do.call(
  rbind,
  lapply(lambda_values, function(lam) {
    calc_lambda_diagnostics_for_subset(
      V_full = V_carn,
      subset_names = disp_subset_names,
      lambda = lam
    )
  })
)

disp_lambda_results$Subset_Type <- "dispersed"
disp_lambda_results$Subset_Size <- subset_size

# 调整列顺序
disp_lambda_results <- disp_lambda_results[
  ,
  c(
    "Subset_Type",
    "Subset_Size",
    "lambda",
    "n_species",
    "off_mean",
    "rmax",
    "neff_mean"
  )
]

disp_lambda_results

# ------------------------------------------------------------
# 27.5 保存结果
# ------------------------------------------------------------

write.csv(
  disp_lambda_results,
  paste0("Carnivora_dispersed_lambda_dependence_s", subset_size, ".csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 27.6 画 lambda sensitivity 图
# ------------------------------------------------------------

plot_df <- disp_lambda_results[
  order(disp_lambda_results$lambda),
]

pdf(
  paste0("Carnivora_dispersed_lambda_dependence_s", subset_size, ".pdf"),
  width = 10,
  height = 4
)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot(
  plot_df$lambda,
  plot_df$off_mean,
  type = "b",
  xlab = expression(lambda),
  ylab = "off_mean",
  main = "Mean off-diagonal correlation",
  xaxt = "n"
)
axis(1, at = sort(lambda_values))

plot(
  plot_df$lambda,
  plot_df$rmax,
  type = "b",
  xlab = expression(lambda),
  ylab = "rmax",
  main = "Maximum off-diagonal correlation",
  xaxt = "n"
)
axis(1, at = sort(lambda_values))

plot(
  plot_df$lambda,
  plot_df$neff_mean,
  type = "b",
  xlab = expression(lambda),
  ylab = "neff_mean",
  main = "Effective sample size",
  xaxt = "n"
)
axis(1, at = sort(lambda_values))

dev.off()

# ------------------------------------------------------------
# 27.7 自动生成可粘贴进正文的句子
# 使用和你当前代码一致的 neff_mean_current_code
# ------------------------------------------------------------

get_lambda_result <- function(df, lam, colname) {
  df[abs(df$lambda - lam) < 1e-10, colname]
}

cat("\nSentence for manuscript using neff_mean_current_code:\n\n")

cat(sprintf(
  paste0(
    "The mean-based effective sample size increased progressively as ",
    "lambda decreased, from %.2f at lambda = 1.00 to %.2f at lambda = 0.75, ",
    "%.2f at lambda = 0.50, %.2f at lambda = 0.25, and %.2f at lambda = 0.10. ",
    "At lambda = 0, neff_mean reached %.2f.\n"
  ),
  get_lambda_result(disp_lambda_results, 1.00, "neff_mean_current_code"),
  get_lambda_result(disp_lambda_results, 0.75, "neff_mean_current_code"),
  get_lambda_result(disp_lambda_results, 0.50, "neff_mean_current_code"),
  get_lambda_result(disp_lambda_results, 0.25, "neff_mean_current_code"),
  get_lambda_result(disp_lambda_results, 0.10, "neff_mean_current_code"),
  get_lambda_result(disp_lambda_results, 0.00, "neff_mean_current_code")
))

# ============================================================
# TEMP Fig. 1d:
# Dispersed subset on 128-tip ladder-like ultrametric tree
# 仅用于临时查看，不影响前面 Fig. 1a-c
# ============================================================

# ------------------------------------------------------------
# 1. 在 ladder-like ultrametric tree 上构造 dispersed subset
# ------------------------------------------------------------

ladder_disp_result <- run_dispersed_fixed_start_algorithm(
  dist_obj = ladder_dist_obj,
  subset_size = sanity_s
)

ladder_disp_idx <- ladder_disp_result$final_subset
ladder_disp_names <- ladder_disp_result$final_subset_names
ladder_disp_metrics <- calc_subset_metrics_extended_local(
  ladder_dist_obj$dist_mat,
  ladder_disp_idx
)

ladder_disp_names
ladder_disp_metrics

# 检查起始物种
ladder_disp_result$start_info$start_name
ladder_disp_result$start_info$mean_distance

# ------------------------------------------------------------
# 2. 保存 ladder dispersed subset 的距离指标和物种列表
# ------------------------------------------------------------

ladder_disp_metrics_table <- data.frame(
  Experiment = "Ladder tree: dispersed subset",
  Subset_Type = "dispersed",
  ladder_disp_metrics,
  stringsAsFactors = FALSE
)

ladder_disp_metrics_table

write.csv(
  ladder_disp_metrics_table,
  "TEMP_Fig1d_ladder_tree_dispersed_metrics_s8.csv",
  row.names = FALSE
)

ladder_disp_selected_tips <- data.frame(
  Experiment = "Ladder tree: dispersed subset",
  Subset_Type = "dispersed",
  Rank = seq_along(ladder_disp_names),
  tiplabel = ladder_disp_names,
  stringsAsFactors = FALSE
)

ladder_disp_selected_tips

write.csv(
  ladder_disp_selected_tips,
  "TEMP_Fig1d_ladder_tree_dispersed_selected_tips_s8.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. 单独输出 Fig. 1d
# 树枝：细黑线
# 选中物种：较大黑点
# 选中枝条：不加粗
# ------------------------------------------------------------

pdf("TEMP_Fig1d_ladder_tree_dispersed_subset_s8.pdf",
    width = 7, height = 7)

par(mar = c(1, 1, 3, 1))

plot_selected_tips_on_tree(
  tree = ladder_tree,
  selected_tips = ladder_disp_names,
  main_title = "d) Ladder tree: dispersed subset"
)

dev.off()

# ------------------------------------------------------------
# 4. 可选：生成临时 a-d 四联图
# 如果你前面已经有 balanced_disp_names, balanced_clust_names,
# ladder_clust_names，则这段可以直接运行
# ------------------------------------------------------------

pdf("TEMP_Figure1_idealized_tree_sanity_checks_with_1d_s8.pdf",
    width = 20, height = 5)

par(mfrow = c(1, 4), mar = c(1, 1, 3, 1))

plot_selected_tips_on_tree(
  tree = balanced_tree,
  selected_tips = balanced_disp_names,
  main_title = "a) Balanced: dispersed"
)

plot_selected_tips_on_tree(
  tree = balanced_tree,
  selected_tips = balanced_clust_names,
  main_title = "b) Balanced: clustered"
)

plot_selected_tips_on_tree(
  tree = ladder_tree,
  selected_tips = ladder_clust_names,
  main_title = "c) Ladder: clustered"
)

plot_selected_tips_on_tree(
  tree = ladder_tree,
  selected_tips = ladder_disp_names,
  main_title = "d) Ladder: dispersed"
)

dev.off()

# ------------------------------------------------------------
# 5. 控制台快速查看
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("TEMP Fig. 1d: Ladder tree dispersed subset\n")
cat("============================================================\n\n")

cat("Ladder tree ultrametric:",
    is.ultrametric(ladder_tree, tol = 1e-8), "\n\n")

cat("Selected tips:\n")
print(ladder_disp_names)

cat("\nDistance metrics:\n")
print(ladder_disp_metrics)

cat("\nStart species:\n")
print(ladder_disp_result$start_info$start_name)
cat("Mean distance of start species:",
    ladder_disp_result$start_info$mean_distance, "\n")

# 当前算法结果
ladder_disp_names
ladder_disp_metrics

# 人工构造你直觉中的 L121-L128
manual_l121_l128 <- paste0("L", sprintf("%03d", 121:128))

manual_l121_l128_idx <- match(manual_l121_l128, ladder_dist_obj$tip_labels)

manual_l121_l128_metrics <- calc_subset_metrics_extended_local(
  ladder_dist_obj$dist_mat,
  manual_l121_l128_idx
)

manual_l121_l128
manual_l121_l128_metrics

# 比较两者指标
rbind(
  Algorithm = ladder_disp_metrics,
  Manual_L121_L128 = manual_l121_l128_metrics
)

# ============================================================
# Figure 2 single-panel versions (density only)
# ------------------------------------------------------------
# Output:
#   Density versions only: 2c, 2d, 2e
#   PDF only
#
# Mapping:
#   Figure 2c = MinPD
#   Figure 2d = MeanPD
#   Figure 2e = MeanNND
#
# Requires existing objects:
#   random_dist_metrics
#   disp_metrics
#   clust_metrics
#   subset_size
# ============================================================

# ---------- 0. Basic checks ----------
stopifnot(exists("random_dist_metrics"))
stopifnot(exists("disp_metrics"))
stopifnot(exists("clust_metrics"))
stopifnot(exists("subset_size"))

# ---------- 1. Global size settings ----------
# Figure 2 单图比 Figure 3 单图略宽
single_panel_width_fig2 <- 4.0
single_panel_height_fig2 <- 4.4

# ---------- 2. Helper: safe range padding ----------
expand_range <- function(x, add_frac = 0.08) {
  xr <- range(x, na.rm = TRUE)
  dx <- diff(xr)
  
  if (!is.finite(dx) || dx == 0) {
    pad <- ifelse(abs(xr[1]) > 0, abs(xr[1]) * add_frac, 0.1)
  } else {
    pad <- dx * add_frac
  }
  
  c(xr[1] - pad, xr[2] + pad)
}

# ---------- 3. Helper: density panel ----------
# obs_disp 和 obs_clust 可以设为 NULL
plot_null_density_single_fig2 <- function(null_values,
                                          obs_disp = NULL,
                                          obs_clust = NULL,
                                          xlab,
                                          panel_label,
                                          ylab = "Density",
                                          adjust_density = 1,
                                          line_lwd = 2,
                                          density_lwd = 1.4,
                                          fill_col = "grey92",
                                          border_col = "grey35") {
  
  dens <- density(null_values, adjust = adjust_density, na.rm = TRUE)
  
  x_all <- dens$x
  if (!is.null(obs_disp)) {
    x_all <- c(x_all, obs_disp)
  }
  if (!is.null(obs_clust)) {
    x_all <- c(x_all, obs_clust)
  }
  
  xlim_use <- expand_range(x_all, add_frac = 0.10)
  ylim_use <- c(0, max(dens$y) * 1.10)
  
  plot(
    dens,
    type = "n",
    xlim = xlim_use,
    ylim = ylim_use,
    xlab = xlab,
    ylab = ylab,
    main = "",
    axes = FALSE,
    yaxs = "i"
  )
  
  axis(1, cex.axis = 0.9)
  axis(2, las = 1, cex.axis = 0.9)
  box()
  
  polygon(
    x = c(dens$x, rev(dens$x)),
    y = c(dens$y, rep(0, length(dens$y))),
    col = fill_col,
    border = NA
  )
  
  lines(dens$x, dens$y, lwd = density_lwd, col = border_col)
  
  if (!is.null(obs_disp)) {
    abline(v = obs_disp, lwd = line_lwd, lty = 1)
  }
  
  if (!is.null(obs_clust)) {
    abline(v = obs_clust, lwd = line_lwd, lty = 2)
  }
  
  mtext(panel_label, side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)
}

# ---------- 3.1 Helper -------------
plot_null_density_single_minpd_hardclip <- function(null_values,
                                                    obs_disp,
                                                    obs_clust = NULL,
                                                    xlab = "MinPD",
                                                    panel_label = "c",
                                                    ylab = "Density",
                                                    adjust_density = 1,
                                                    line_lwd = 2,
                                                    density_lwd = 1.4,
                                                    fill_col = "grey92",
                                                    border_col = "grey35") {
  
  dens_raw <- density(
    null_values,
    adjust = adjust_density,
    na.rm = TRUE
  )
  
  # 硬裁剪：只保留 0 到 30 之间的密度曲线
  keep <- dens_raw$x >= 0 & dens_raw$x <= 30
  dens_x <- dens_raw$x[keep]
  dens_y <- dens_raw$y[keep]
  
  if (length(dens_x) == 0) {
    stop("No density x values remain after clipping to [0, 30].")
  }
  
  # 补左端点 0
  if (min(dens_x) > 0) {
    y0 <- approx(dens_raw$x, dens_raw$y, xout = 0, rule = 2)$y
    dens_x <- c(0, dens_x)
    dens_y <- c(y0, dens_y)
  }
  
  # 补右端点 30
  if (max(dens_x) < 30) {
    y30 <- approx(dens_raw$x, dens_raw$y, xout = 30, rule = 2)$y
    dens_x <- c(dens_x, 30)
    dens_y <- c(dens_y, y30)
  }
  
  dens_x[1] <- 0
  dens_x[length(dens_x)] <- 30
  
  ylim_use <- c(0, max(dens_y, na.rm = TRUE) * 1.10)
  
  obs_disp_use <- min(max(obs_disp, 0), 30)
  
  if (!is.null(obs_clust)) {
    obs_clust_use <- min(max(obs_clust, 0), 30)
  }
  
  old_xpd <- par("xpd")
  par(xpd = FALSE)
  on.exit(par(xpd = old_xpd), add = TRUE)
  
  plot(
    NA,
    xlim = c(0, 30),
    ylim = ylim_use,
    xlab = xlab,
    ylab = ylab,
    main = "",
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  
  axis(1, at = seq(0, 30, by = 10), cex.axis = 0.9)
  axis(2, las = 1, cex.axis = 0.9)
  box()
  
  polygon(
    x = c(0, dens_x, 30),
    y = c(0, dens_y, 0),
    col = fill_col,
    border = NA
  )
  
  lines(
    dens_x,
    dens_y,
    lwd = density_lwd,
    col = border_col
  )
  
  segments(
    x0 = obs_disp_use,
    y0 = 0,
    x1 = obs_disp_use,
    y1 = ylim_use[2],
    lwd = line_lwd,
    lty = 1
  )
  
  if (!is.null(obs_clust)) {
    segments(
      x0 = obs_clust_use,
      y0 = 0,
      x1 = obs_clust_use,
      y1 = ylim_use[2],
      lwd = line_lwd,
      lty = 2
    )
  }
  
  mtext(panel_label, side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)
}
# ---------- 4. Common par settings ----------
set_single_panel_par_fig2 <- function() {
  par(
    mar = c(3.8, 4.0, 1.8, 0.8),
    mgp = c(2.0, 0.65, 0),
    tcl = -0.25,
    cex = 0.95
  )
}

# ============================================================
# 5. Density versions only
# ============================================================

# ---- Figure 2c: MinPD (dispersed only) ----
pdf(
  paste0("Figure2c_MinPD_density_s", subset_size, ".pdf"),
  width = single_panel_width_fig2,
  height = single_panel_height_fig2,
  pointsize = 9
)
set_single_panel_par_fig2()
plot_null_density_single_minpd_hardclip(
  null_values = random_dist_metrics$MinPD,
  obs_disp = disp_metrics$MinPD,
  obs_clust = NULL,
  xlab = "MinPD",
  panel_label = "c"
)
legend(
  "top",
  legend = c("Dispersed"),
  lty = c(1),
  lwd = 2,
  bty = "n",
  cex = 0.85,
  seg.len = 4
)
dev.off()

# ---- Figure 2d: MeanPD (dispersed + clustered) ----
pdf(
  paste0("Figure2d_MeanPD_density_s", subset_size, ".pdf"),
  width = single_panel_width_fig2,
  height = single_panel_height_fig2,
  pointsize = 9
)
set_single_panel_par_fig2()
plot_null_density_single_fig2(
  null_values = random_dist_metrics$MeanPD,
  obs_disp = disp_metrics$MeanPD,
  obs_clust = clust_metrics$MeanPD,
  xlab = "MeanPD",
  panel_label = "d"
)
legend(
  "top",
  legend = c("Dispersed", "Clustered"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n",
  cex = 0.85,
  seg.len = 4
)
dev.off()

# ---- Figure 2e: MeanNND (dispersed + clustered) ----
pdf(
  paste0("Figure2e_MeanNND_density_s", subset_size, ".pdf"),
  width = single_panel_width_fig2,
  height = single_panel_height_fig2,
  pointsize = 9
)
set_single_panel_par_fig2()
plot_null_density_single_fig2(
  null_values = random_dist_metrics$MeanNND,
  obs_disp = disp_metrics$MeanNND,
  obs_clust = clust_metrics$MeanNND,
  xlab = "MeanNND",
  panel_label = "e"
)
legend(
  "top",
  legend = c("Dispersed", "Clustered"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n",
  cex = 0.85,
  seg.len = 4
)
dev.off()

# ============================================================
# 6. Console summary
# ============================================================

cat("\n============================================================\n")
cat("Single-panel Figure 2 PDF files written.\n")
cat("Density versions only:\n")
cat("  ", paste0("Figure2c_MinPD_density_s", subset_size, ".pdf"), "\n")
cat("  ", paste0("Figure2d_MeanPD_density_s", subset_size, ".pdf"), "\n")
cat("  ", paste0("Figure2e_MeanNND_density_s", subset_size, ".pdf"), "\n")
cat("============================================================\n\n")

cat("Observed distance metrics:\n")
print(data.frame(
  Subset = c("Dispersed", "Clustered"),
  MinPD = c(disp_metrics$MinPD, NA),
  MeanPD = c(disp_metrics$MeanPD, clust_metrics$MeanPD),
  MeanNND = c(disp_metrics$MeanNND, clust_metrics$MeanNND)
))


# ============================================================
# Figure 3 single-panel versions (density only)
# ------------------------------------------------------------
# Output:
#   Density versions only: 3a, 3b, 3c
#   Panel 3d line plot
#   PDF only
#
# Requires existing objects:
#   random_dep_metrics
#   disp_dep
#   clust_dep
#   V_carn
#   disp_subset_names
#   subset_size
# ============================================================

# ---------- 0. Basic checks ----------
stopifnot(exists("random_dep_metrics"))
stopifnot(exists("disp_dep"))
stopifnot(exists("clust_dep"))
stopifnot(exists("V_carn"))
stopifnot(exists("disp_subset_names"))
stopifnot(exists("subset_size"))

# ---------- 1. Global size settings ----------
# 保持与原 Figure 3 子图相同的宽高比例（约 1.6 : 2.2）
single_panel_width <- 3.2
single_panel_height <- 4.4

# ---------- 2. Helper: safe range padding ----------
expand_range <- function(x, add_frac = 0.08) {
  xr <- range(x, na.rm = TRUE)
  dx <- diff(xr)
  
  if (!is.finite(dx) || dx == 0) {
    pad <- ifelse(abs(xr[1]) > 0, abs(xr[1]) * add_frac, 0.1)
  } else {
    pad <- dx * add_frac
  }
  
  c(xr[1] - pad, xr[2] + pad)
}

# ---------- 3. Helper: density panel ----------
plot_null_density_single <- function(null_values,
                                     obs_disp,
                                     obs_clust,
                                     xlab,
                                     panel_label,
                                     ylab = "Density",
                                     adjust_density = 1,
                                     line_lwd = 2,
                                     density_lwd = 1.4,
                                     fill_col = "grey92",
                                     border_col = "grey35",
                                     xlim_override = NULL,
                                     xaxs_style = "r",
                                     density_from = NULL,
                                     density_to = NULL,
                                     density_cut = 3) {
  
  # density() 不能显式传入 from = NULL 或 to = NULL
  # 所以这里分情况处理：
  if (is.null(density_from) && is.null(density_to)) {
    
    dens <- density(
      null_values,
      adjust = adjust_density,
      na.rm = TRUE,
      cut = density_cut
    )
    
  } else {
    
    from_use <- if (is.null(density_from)) {
      min(null_values, na.rm = TRUE)
    } else {
      density_from
    }
    
    to_use <- if (is.null(density_to)) {
      max(null_values, na.rm = TRUE)
    } else {
      density_to
    }
    
    dens <- density(
      null_values,
      adjust = adjust_density,
      na.rm = TRUE,
      from = from_use,
      to = to_use,
      cut = density_cut
    )
  }
  
  x_all <- c(dens$x, obs_disp, obs_clust)
  
  if (is.null(xlim_override)) {
    xlim_use <- expand_range(x_all, add_frac = 0.10)
  } else {
    xlim_use <- xlim_override
  }
  
  ylim_use <- c(0, max(dens$y) * 1.10)
  
  plot(
    dens,
    type = "n",
    xlim = xlim_use,
    ylim = ylim_use,
    xlab = xlab,
    ylab = ylab,
    main = "",
    axes = FALSE,
    xaxs = xaxs_style,
    yaxs = "i"
  )
  
  axis(1, cex.axis = 0.9)
  axis(2, las = 1, cex.axis = 0.9)
  box()
  
  polygon(
    x = c(dens$x, rev(dens$x)),
    y = c(dens$y, rep(0, length(dens$y))),
    col = fill_col,
    border = NA
  )
  
  lines(dens$x, dens$y, lwd = density_lwd, col = border_col)
  
  # dispersed: solid
  abline(v = obs_disp, lwd = line_lwd, lty = 1)
  
  # clustered: dashed
  abline(v = obs_clust, lwd = line_lwd, lty = 2)
  
  mtext(panel_label, side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)
}
# ---------- 3.1 Helper: MaxOffCor density panel with hard clipping to [0, 1] ----------
plot_null_density_single_rmax_hardclip <- function(null_values,
                                                   obs_disp,
                                                   obs_clust,
                                                   xlab = "MaxOffCor",
                                                   panel_label = "b",
                                                   ylab = "Density",
                                                   adjust_density = 1,
                                                   line_lwd = 2,
                                                   density_lwd = 1.4,
                                                   fill_col = "grey92",
                                                   border_col = "grey35") {
  
  # 先计算密度；允许 density 自己平滑
  dens_raw <- density(
    null_values,
    adjust = adjust_density,
    na.rm = TRUE
  )
  
  # 硬裁剪：只保留 0.6 到 1 之间的密度曲线
  keep <- dens_raw$x >= 0.6 & dens_raw$x <= 1
  dens_x <- dens_raw$x[keep]
  dens_y <- dens_raw$y[keep]
  
  # 如果裁剪后没有点，报错
  if (length(dens_x) == 0) {
    stop("No density x values remain after clipping to [0.6, 1].")
  }
  
  # 为了保证边界正好从 0.6 到 1，补上端点
  if (min(dens_x) > 0.6) {
    y0 <- approx(dens_raw$x, dens_raw$y, xout = 0.6, rule = 2)$y
    dens_x <- c(0.6, dens_x)
    dens_y <- c(y0, dens_y)
  }
  
  if (max(dens_x) < 1) {
    y1 <- approx(dens_raw$x, dens_raw$y, xout = 1, rule = 2)$y
    dens_x <- c(dens_x, 1)
    dens_y <- c(dens_y, y1)
  }
  
  # 强制端点精确为 0.6 和 1
  dens_x[1] <- 0.6
  dens_x[length(dens_x)] <- 1
  
  ylim_use <- c(0, max(dens_y, na.rm = TRUE) * 1.10)
  
  # obs 也限制在 [0.6, 1]，避免极小浮点误差造成越界
  obs_disp_use <- min(max(obs_disp, 0.6), 1)
  obs_clust_use <- min(max(obs_clust, 0.6), 1)
  
  old_xpd <- par("xpd")
  par(xpd = FALSE)
  on.exit(par(xpd = old_xpd), add = TRUE)
  
  plot(
    NA,
    xlim = c(0.6, 1),
    ylim = ylim_use,
    xlab = xlab,
    ylab = ylab,
    main = "",
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  
  axis(1, at = seq(0.6, 1, by = 0.1), labels = sprintf("%.1f", seq(0.6, 1, by = 0.1)), cex.axis = 0.9)
  axis(2, las = 1, cex.axis = 0.9)
  box()
  
  polygon(
    x = c(0.6, dens_x, 1),
    y = c(0, dens_y, 0),
    col = fill_col,
    border = NA
  )
  
  lines(
    dens_x,
    dens_y,
    lwd = density_lwd,
    col = border_col
  )
  
  # 用 segments 代替 abline，严格限制在绘图区内
  segments(
    x0 = obs_disp_use,
    y0 = 0,
    x1 = obs_disp_use,
    y1 = ylim_use[2],
    lwd = line_lwd,
    lty = 1
  )
  
  segments(
    x0 = obs_clust_use,
    y0 = 0,
    x1 = obs_clust_use,
    y1 = ylim_use[2],
    lwd = line_lwd,
    lty = 2
  )
  
  mtext(panel_label, side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)
}# ---------- 4. Helper: lambda-transformed covariance ----------
lambda_transform_V <- function(V, lambda) {
  V_lambda <- V
  off_diag <- row(V_lambda) != col(V_lambda)
  V_lambda[off_diag] <- lambda * V_lambda[off_diag]
  V_lambda
}

# ---------- 5. Helper: dependence metrics from covariance ----------
calc_dependence_from_V_simple <- function(V_full, subset_names) {
  
  V_sub <- V_full[subset_names, subset_names, drop = FALSE]
  R_sub <- cov2cor(V_sub)
  
  off_values <- R_sub[upper.tri(R_sub)]
  
  n <- length(subset_names)
  one_vec <- rep(1, n)
  
  neff_mean <- as.numeric(n^2 / (t(one_vec) %*% R_sub %*% one_vec))
  
  data.frame(
    off_mean = mean(off_values),
    rmax = max(off_values),
    neff_mean = neff_mean
  )
}

# ---------- 6. Lambda sensitivity for panel 3d ----------
lambda_grid <- seq(0, 1, by = 0.1)

lambda_meanESS <- sapply(lambda_grid, function(lam) {
  V_lam <- lambda_transform_V(V_carn, lam)
  calc_dependence_from_V_simple(V_lam, disp_subset_names)$neff_mean
})

lambda_meanESS_df <- data.frame(
  lambda = lambda_grid,
  MeanESS = as.numeric(lambda_meanESS)
)

write.csv(
  lambda_meanESS_df,
  paste0("Figure3d_lambda_meanESS_dispersed_s", subset_size, ".csv"),
  row.names = FALSE
)

# ---------- 7. Common par settings ----------
set_single_panel_par <- function() {
  par(
    mar = c(3.8, 4.0, 1.8, 0.8),
    mgp = c(2.0, 0.65, 0),
    tcl = -0.25,
    cex = 0.95
  )
}

# ============================================================
# 8. Density versions only: panels a, b, c
# ============================================================

# ---- 3a density: MeanOffCor ----
pdf(
  paste0("Figure3a_MeanOffCor_density_s", subset_size, ".pdf"),
  width = single_panel_width,
  height = single_panel_height,
  pointsize = 9
)
set_single_panel_par()
plot_null_density_single(
  null_values = random_dep_metrics$off_mean,
  obs_disp = disp_dep$off_mean,
  obs_clust = clust_dep$off_mean,
  xlab = "MeanOffCor",
  panel_label = "a"
)
legend(
  "top",
  legend = c("Dispersed", "Clustered"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n",
  cex = 0.85,
  seg.len = 4
)
dev.off()

# ---- 3b density: MaxOffCor ----
# 这里对核密度进行边界限制，并把横轴固定在 0 到 1
pdf(
  paste0("Figure3b_MaxOffCor_density_s", subset_size, ".pdf"),
  width = single_panel_width,
  height = single_panel_height,
  pointsize = 9
)
set_single_panel_par()
plot_null_density_single_rmax_hardclip(
  null_values = random_dep_metrics$rmax,
  obs_disp = disp_dep$rmax,
  obs_clust = clust_dep$rmax,
  xlab = "MaxOffCor",
  panel_label = "b"
)
legend(
  "top",
  legend = c("Dispersed", "Clustered"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n",
  cex = 0.85,
  seg.len = 4
)
dev.off()

# ---- 3c density: MeanESS ----
pdf(
  paste0("Figure3c_MeanESS_density_s", subset_size, ".pdf"),
  width = single_panel_width,
  height = single_panel_height,
  pointsize = 9
)
set_single_panel_par()
plot_null_density_single(
  null_values = random_dep_metrics$neff_mean,
  obs_disp = disp_dep$neff_mean,
  obs_clust = clust_dep$neff_mean,
  xlab = "MeanESS",
  panel_label = "c"
)
legend(
  "top",
  legend = c("Dispersed", "Clustered"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n",
  cex = 0.85,
  seg.len = 4
)
dev.off()

# ============================================================
# 9. Panel 3d single version
# ============================================================

ylim_d <- expand_range(c(lambda_meanESS_df$MeanESS, subset_size), add_frac = 0.08)

pdf(
  paste0("Figure3d_lambda_MeanESS_s", subset_size, ".pdf"),
  width = single_panel_width,
  height = single_panel_height,
  pointsize = 9
)
set_single_panel_par()
plot(
  lambda_meanESS_df$lambda,
  lambda_meanESS_df$MeanESS,
  type = "b",
  pch = 16,
  lwd = 1.4,
  xlim = c(0, 1),
  ylim = ylim_d,
  axes = FALSE,
  xlab = expression(lambda),
  ylab = "MeanESS",
  main = ""
)
axis(1, at = seq(0, 1, by = 0.1), labels = sprintf("%.1f", seq(0, 1, by = 0.1)), cex.axis = 0.9)
axis(2, las = 1, cex.axis = 0.9)
box()

# 如需预览参考线，取消注释即可
# abline(h = subset_size, lty = 3, lwd = 1.2)

mtext("d", side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)
dev.off()


# ============================================================
# 1. Panel wrapper functions
#    These wrappers reproduce the original single-panel plots
#    without changing the original plotting style.
# ============================================================

draw_fig2c_panel <- function() {
  plot_null_density_single_minpd_hardclip(
    null_values = random_dist_metrics$MinPD,
    obs_disp = disp_metrics$MinPD,
    obs_clust = NULL,
    xlab = "MinPD",
    panel_label = "c"
  )
  legend(
    "top",
    legend = c("Dispersed"),
    lty = c(1),
    lwd = 2,
    bty = "n",
    cex = 0.85,
    seg.len = 4
  )
}

draw_fig2d_panel <- function() {
  plot_null_density_single_fig2(
    null_values = random_dist_metrics$MeanPD,
    obs_disp = disp_metrics$MeanPD,
    obs_clust = clust_metrics$MeanPD,
    xlab = "MeanPD",
    panel_label = "d"
  )
  legend(
    "top",
    legend = c("Dispersed", "Clustered"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n",
    cex = 0.85,
    seg.len = 4
  )
}

draw_fig2e_panel <- function() {
  plot_null_density_single_fig2(
    null_values = random_dist_metrics$MeanNND,
    obs_disp = disp_metrics$MeanNND,
    obs_clust = clust_metrics$MeanNND,
    xlab = "MeanNND",
    panel_label = "e"
  )
  legend(
    "top",
    legend = c("Dispersed", "Clustered"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n",
    cex = 0.85,
    seg.len = 4
  )
}

draw_fig3a_panel <- function() {
  plot_null_density_single(
    null_values = random_dep_metrics$off_mean,
    obs_disp = disp_dep$off_mean,
    obs_clust = clust_dep$off_mean,
    xlab = "MeanOffCor",
    panel_label = "a"
  )
  legend(
    "top",
    legend = c("Dispersed", "Clustered"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n",
    cex = 0.85,
    seg.len = 4
  )
}

draw_fig3b_panel <- function() {
  plot_null_density_single_rmax_hardclip(
    null_values = random_dep_metrics$rmax,
    obs_disp = disp_dep$rmax,
    obs_clust = clust_dep$rmax,
    xlab = "MaxOffCor",
    panel_label = "b"
  )
  legend(
    "top",
    legend = c("Dispersed", "Clustered"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n",
    cex = 0.85,
    seg.len = 4
  )
}

draw_fig3c_panel <- function() {
  plot_null_density_single(
    null_values = random_dep_metrics$neff_mean,
    obs_disp = disp_dep$neff_mean,
    obs_clust = clust_dep$neff_mean,
    xlab = "MeanESS",
    panel_label = "c"
  )
  legend(
    "top",
    legend = c("Dispersed", "Clustered"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n",
    cex = 0.85,
    seg.len = 4
  )
}

draw_fig3d_panel <- function() {
  lambda_grid <- seq(0, 1, by = 0.1)
  
  lambda_meanESS <- sapply(lambda_grid, function(lam) {
    V_lam <- lambda_transform_V(V_carn, lam)
    calc_dependence_from_V_simple(V_lam, disp_subset_names)$neff_mean
  })
  
  lambda_meanESS_df <- data.frame(
    lambda = lambda_grid,
    MeanESS = as.numeric(lambda_meanESS)
  )
  
  ylim_d <- expand_range(c(lambda_meanESS_df$MeanESS, subset_size), add_frac = 0.08)
  
  plot(
    lambda_meanESS_df$lambda,
    lambda_meanESS_df$MeanESS,
    type = "b",
    pch = 16,
    lwd = 1.4,
    xlim = c(0, 1),
    ylim = ylim_d,
    axes = FALSE,
    xlab = expression(lambda),
    ylab = "MeanESS",
    main = ""
  )
  axis(
    1,
    at = seq(0, 1, by = 0.1),
    labels = sprintf("%.1f", seq(0, 1, by = 0.1)),
    cex.axis = 0.9
  )
  axis(2, las = 1, cex.axis = 0.9)
  box()
  
  # 如需参考线，取消注释
  # abline(h = subset_size, lty = 3, lwd = 1.2)
  
  mtext("d", side = 3, adj = 0, line = 0.2, font = 2, cex = 1.0)
}

# ============================================================
# 2. Common par settings for combined figures
#    (basically same style as your single-panel figures)
# ============================================================

set_combined_par <- function() {
  par(
    mar = c(3.8, 4.0, 1.8, 0.8),
    mgp = c(2.0, 0.65, 0),
    tcl = -0.25,
    cex = 0.95
  )
}

# ============================================================
# 3. Combined Figure 2: c + d + e
#    Layout: 1 row x 3 columns
# ============================================================

pdf(
  paste0("Figure2_cde_combined_direct_s", subset_size, ".pdf"),
  width = 12,
  height = 4.4,
  pointsize = 9
)

par(mfrow = c(1, 3))
set_combined_par()

draw_fig2c_panel()
draw_fig2d_panel()
draw_fig2e_panel()

dev.off()

# ============================================================
# 4. Combined Figure 3: a + b + c + d
#    Layout: 1 row x 4 columns
# ============================================================

pdf(
  paste0("Figure3_abcd_combined_direct_s", subset_size, ".pdf"),
  width = 14.4,
  height = 4.4,
  pointsize = 9
)

par(mfrow = c(1, 4))
set_combined_par()

draw_fig3a_panel()
draw_fig3b_panel()
draw_fig3c_panel()
draw_fig3d_panel()

dev.off()

# 用 layout() 来实现“上面 3 张、下面 4 张”的不等列排版
# 这里把整页想象成 12 列：
# 上排每张占 4 列，共 3 张
# 下排每张占 3 列，共 4 张

layout_mat <- matrix(
  c(
    1,1,1,1,  2,2,2,2,  3,3,3,3,
    4,4,4,    5,5,5,    6,6,6,    7,7,7
  ),
  nrow = 2,
  byrow = TRUE
)

pdf(
  paste0("Figure2cde_Figure3abcd_combined_direct_s", subset_size, ".pdf"),
  width = 15,
  height = 8.6,
  pointsize = 9
)

layout(layout_mat)
set_combined_par()

# ---- top row: Figure 2 c, d, e ----
draw_fig2c_panel()
draw_fig2d_panel()
draw_fig2e_panel()

# ---- bottom row: Figure 3 a, b, c, d ----
draw_fig3a_panel()
draw_fig3b_panel()
draw_fig3c_panel()
draw_fig3d_panel()

dev.off()
