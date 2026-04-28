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

source("D:/R/distance_metrics.R")
source("D:/R/objective_compare.R")
source("D:/R/subset_greedy.R")
source("D:/R/subset_exchange.R")
source("D:/R/subset_random.R")

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
# ============================================================

set.seed(2)

disp_result <- run_complete_algorithm(
  dist_obj = dist_obj,
  subset_size = subset_size,
  maximize = TRUE,
  n_greedy_starts = n_greedy_starts
)

disp_subset_idx <- disp_result$final_subset
disp_subset_names <- disp_result$final_subset_names
disp_metrics <- disp_result$final_metrics

length(disp_subset_idx)
disp_metrics
disp_subset_names

# ============================================================
# 9. 构造 phylogenetically clustered subset
# ============================================================

set.seed(1)

clust_result <- run_complete_algorithm(
  dist_obj = dist_obj,
  subset_size = subset_size,
  maximize = FALSE,
  n_greedy_starts = n_greedy_starts
)

clust_subset_idx <- clust_result$final_subset
clust_subset_names <- clust_result$final_subset_names
clust_metrics <- clust_result$final_metrics

length(clust_subset_idx)
clust_metrics
clust_subset_names

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
# ============================================================

random_dist_metrics <- calc_multiple_subsets_metrics(
  dist_mat = dist_obj$dist_mat,
  subsets = random_subsets_idx
)

head(random_dist_metrics)
summary(random_dist_metrics)

# ============================================================
# 12. 检查 observed distance metrics
# ============================================================

disp_metrics <- calc_subset_metrics(dist_obj$dist_mat, disp_subset_idx)
clust_metrics <- calc_subset_metrics(dist_obj$dist_mat, clust_subset_idx)

disp_metrics
clust_metrics

# ============================================================
# 13. 距离指标的 p-value 和 SES
# ============================================================

calc_ses <- function(obs, null_values) {
  (obs - mean(null_values)) / sd(null_values)
}

calc_p_high <- function(obs, null_values) {
  (1 + sum(null_values >= obs)) / (1 + length(null_values))
}

calc_p_low <- function(obs, null_values) {
  (1 + sum(null_values <= obs)) / (1 + length(null_values))
}
# dispersed: high-tail p-values
disp_p_MinPD <- calc_p_high(disp_metrics$MinPD, random_dist_metrics$MinPD)
disp_p_MeanPD <- calc_p_high(disp_metrics$MeanPD, random_dist_metrics$MeanPD)
disp_p_MeanNND <- calc_p_high(disp_metrics$MeanNND, random_dist_metrics$MeanNND)

# clustered: low-tail p-values
clust_p_MinPD <- calc_p_low(clust_metrics$MinPD, random_dist_metrics$MinPD)
clust_p_MeanPD <- calc_p_low(clust_metrics$MeanPD, random_dist_metrics$MeanPD)
clust_p_MeanNND <- calc_p_low(clust_metrics$MeanNND, random_dist_metrics$MeanNND)

# SES: 统一用 (observed - random mean) / random sd
disp_ses_MinPD <- calc_ses(disp_metrics$MinPD, random_dist_metrics$MinPD)
disp_ses_MeanPD <- calc_ses(disp_metrics$MeanPD, random_dist_metrics$MeanPD)
disp_ses_MeanNND <- calc_ses(disp_metrics$MeanNND, random_dist_metrics$MeanNND)

clust_ses_MinPD <- calc_ses(clust_metrics$MinPD, random_dist_metrics$MinPD)
clust_ses_MeanPD <- calc_ses(clust_metrics$MeanPD, random_dist_metrics$MeanPD)
clust_ses_MeanNND <- calc_ses(clust_metrics$MeanNND, random_dist_metrics$MeanNND)

# 汇总距离结果
distance_summary <- data.frame(
  Metric = c("MinPD", "MeanPD", "MeanNND"),
  
  Dispersed = c(
    disp_metrics$MinPD,
    disp_metrics$MeanPD,
    disp_metrics$MeanNND
  ),
  
  Clustered = c(
    clust_metrics$MinPD,
    clust_metrics$MeanPD,
    clust_metrics$MeanNND
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
  
  Dispersed_SES = c(
    disp_ses_MinPD,
    disp_ses_MeanPD,
    disp_ses_MeanNND
  ),
  
  Clustered_SES = c(
    clust_ses_MinPD,
    clust_ses_MeanPD,
    clust_ses_MeanNND
  ),
  
  Dispersed_P_high = c(
    disp_p_MinPD,
    disp_p_MeanPD,
    disp_p_MeanNND
  ),
  
  Clustered_P_low = c(
    clust_p_MinPD,
    clust_p_MeanPD,
    clust_p_MeanNND
  )
)

distance_summary

write.csv(distance_summary,
          "Carnivora_distance_summary_s30.csv",
          row.names = FALSE)

distance_summary

write.csv(distance_summary,
          "Carnivora_distance_summary_s30_2.csv",
          row.names = FALSE)

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
  
  neff_mean <- as.numeric(n^2 / (t(one_vec) %*% R_sub %*% one_vec))
  
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
          "Carnivora_random_dependence_metrics_s30.csv",
          row.names = FALSE)

# ============================================================
# 18. dependence diagnostics 的 p-value 和 SES
# ============================================================

# dispersed
disp_p_off <- calc_p_low(disp_dep$off_mean, random_dep_metrics$off_mean)
disp_p_rmax <- calc_p_low(disp_dep$rmax, random_dep_metrics$rmax)
disp_p_neff <- calc_p_high(disp_dep$neff_mean, random_dep_metrics$neff_mean)

# clustered
clust_p_off <- calc_p_high(clust_dep$off_mean, random_dep_metrics$off_mean)
clust_p_rmax <- calc_p_high(clust_dep$rmax, random_dep_metrics$rmax)
clust_p_neff <- calc_p_low(clust_dep$neff_mean, random_dep_metrics$neff_mean)

# SES: 仍然统一用 observed - random mean / random sd
disp_ses_off <- calc_ses(disp_dep$off_mean, random_dep_metrics$off_mean)
disp_ses_rmax <- calc_ses(disp_dep$rmax, random_dep_metrics$rmax)
disp_ses_neff <- calc_ses(disp_dep$neff_mean, random_dep_metrics$neff_mean)

clust_ses_off <- calc_ses(clust_dep$off_mean, random_dep_metrics$off_mean)
clust_ses_rmax <- calc_ses(clust_dep$rmax, random_dep_metrics$rmax)
clust_ses_neff <- calc_ses(clust_dep$neff_mean, random_dep_metrics$neff_mean)

dependence_summary <- data.frame(
  Metric = c("off_mean", "rmax", "neff_mean"),
  
  Dispersed = c(
    disp_dep$off_mean,
    disp_dep$rmax,
    disp_dep$neff_mean
  ),
  
  Clustered = c(
    clust_dep$off_mean,
    clust_dep$rmax,
    clust_dep$neff_mean
  ),
  
  Random_Mean = c(
    mean(random_dep_metrics$off_mean),
    mean(random_dep_metrics$rmax),
    mean(random_dep_metrics$neff_mean)
  ),
  
  Random_SD = c(
    sd(random_dep_metrics$off_mean),
    sd(random_dep_metrics$rmax),
    sd(random_dep_metrics$neff_mean)
  ),
  
  Dispersed_SES = c(
    disp_ses_off,
    disp_ses_rmax,
    disp_ses_neff
  ),
  
  Clustered_SES = c(
    clust_ses_off,
    clust_ses_rmax,
    clust_ses_neff
  ),
  
  Dispersed_P_favorable = c(
    disp_p_off,
    disp_p_rmax,
    disp_p_neff
  ),
  
  Clustered_P_favorable = c(
    clust_p_off,
    clust_p_rmax,
    clust_p_neff
  )
)

dependence_summary

write.csv(dependence_summary,
          "Carnivora_dependence_summary_s30.csv",
          row.names = FALSE)

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
          "Carnivora_selected_species_taxonomy_s30.csv",
          row.names = FALSE)
# ============================================================
# 24. 树上标出 dispersed subset
# ============================================================

svg("Carnivora_tree_dispersed_subset_s20_2.svg",
    width = 9, height = 9)

plot.phylo(carn_tree,
           type = "fan",
           show.tip.label = FALSE,
           no.margin = FALSE,
           cex = 0.3)

tiplabels(
  pch = 21,
  bg = "black",
  cex = 0.9,
  tip = which(carn_tree$tip.label %in% disp_subset_names)
)

title("Phylogenetically dispersed subset, s = 20")

dev.off()

# ============================================================
# 25. 树上标出 clustered subset
# ============================================================

svg("Carnivora_tree_clustered_subset_s20_1.svg",
    width = 9, height = 9)

plot.phylo(carn_tree,
           type = "fan",
           show.tip.label = FALSE,
           no.margin = FALSE,
           cex = 0.3)

tiplabels(
  pch = 21,
  bg = "black",
  cex = 0.9,
  tip = which(carn_tree$tip.label %in% clust_subset_names)
)

title("Phylogenetically clustered subset, s = 20")

dev.off()


# ============================================================
# TEMP analysis:
# 从 1000 个 random subsets 中找出 4 个目标子集
# ============================================================

# 给 random metrics 加一个 Subset_ID，便于回溯
random_summary <- data.frame(
  Subset_ID = seq_len(nrow(random_dist_metrics)),
  random_dist_metrics,
  stringsAsFactors = FALSE
)

head(random_summary)

# ============================================================
# 1. 三个单指标最小子集
# ============================================================

min_MinPD_id   <- which.min(random_summary$MinPD)
min_MeanPD_id  <- which.min(random_summary$MeanPD)
min_MeanNND_id <- which.min(random_summary$MeanNND)

min_MinPD_id
min_MeanPD_id
min_MeanNND_id

# ============================================================
# 2. 字典序最优子集（clustered 方向：越小越好）
# 先比 MinPD，再比 MeanPD，再比 MeanNND
# ============================================================

lexico_order_ids <- order(
  random_summary$MinPD,
  random_summary$MeanPD,
  random_summary$MeanNND
)

# 字典序最优的 random subset
lexico_best_id <- lexico_order_ids[1]

lexico_best_id

# 看一下字典序前 10 名
random_summary[lexico_order_ids[1:10], ]

# ============================================================
# 3. 提取子集信息的辅助函数
# ============================================================

extract_random_subset_info <- function(subset_id, label) {
  
  subset_idx <- random_subsets_idx[[subset_id]]
  subset_names <- dist_obj$tip_labels[subset_idx]
  subset_metrics <- random_summary[random_summary$Subset_ID == subset_id, ]
  
  list(
    Label = label,
    Subset_ID = subset_id,
    subset_idx = subset_idx,
    subset_names = subset_names,
    metrics = subset_metrics
  )
}

# ============================================================
# 4. 提取四个目标子集
# ============================================================

subset_min_MinPD <- extract_random_subset_info(
  subset_id = min_MinPD_id,
  label = "Random subset with minimum MinPD"
)

subset_min_MeanPD <- extract_random_subset_info(
  subset_id = min_MeanPD_id,
  label = "Random subset with minimum MeanPD"
)

subset_min_MeanNND <- extract_random_subset_info(
  subset_id = min_MeanNND_id,
  label = "Random subset with minimum MeanNND"
)

subset_lexico_best <- extract_random_subset_info(
  subset_id = lexico_best_id,
  label = "Random subset with lexicographic optimum"
)

# ============================================================
# 5. 汇总成总表
# ============================================================

selected_random_subsets_summary <- rbind(
  data.frame(
    Type = subset_min_MinPD$Label,
    subset_min_MinPD$metrics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Type = subset_min_MeanPD$Label,
    subset_min_MeanPD$metrics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Type = subset_min_MeanNND$Label,
    subset_min_MeanNND$metrics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Type = subset_lexico_best$Label,
    subset_lexico_best$metrics,
    stringsAsFactors = FALSE
  )
)

selected_random_subsets_summary

# ============================================================
# 8. 画图函数
# ============================================================

plot_subset_on_tree <- function(tree, subset_names, file_name, main_title) {
  
  svg(file_name, width = 9, height = 9)
  
  par(mar = c(1, 1, 3, 1))
  
  plot.phylo(tree,
             type = "fan",
             show.tip.label = FALSE,
             no.margin = FALSE,
             cex = 0.3)
  
  tiplabels(
    pch = 21,
    bg = "black",
    cex = 0.9,
    tip = which(tree$tip.label %in% subset_names)
  )
  
  title(main_title, line = 1)
  
  dev.off()
}

# ============================================================
# 9. 分别画四个子集
# ============================================================

plot_subset_on_tree(
  tree = carn_tree,
  subset_names = subset_min_MinPD$subset_names,
  file_name = "TEMP_random_subset_min_MinPD_s20.svg",
  main_title = "Random subset with minimum MinPD, s = 20"
)

plot_subset_on_tree(
  tree = carn_tree,
  subset_names = subset_min_MeanPD$subset_names,
  file_name = "TEMP_random_subset_min_MeanPD_s20.svg",
  main_title = "Random subset with minimum MeanPD, s = 20"
)

plot_subset_on_tree(
  tree = carn_tree,
  subset_names = subset_min_MeanNND$subset_names,
  file_name = "TEMP_random_subset_min_MeanNND_s20.svg",
  main_title = "Random subset with minimum MeanNND, s = 20"
)

plot_subset_on_tree(
  tree = carn_tree,
  subset_names = subset_lexico_best$subset_names,
  file_name = "TEMP_random_subset_lexico_best_s20.svg",
  main_title = "Random subset with lexicographic optimum, s = 20"
)


# ============================================================
# TEMP function 1:
# 找出全树中 patristic distance 最小的一对物种
# ============================================================

find_closest_pair <- function(dist_obj) {
  
  d <- dist_obj$dist_mat
  
  # 排除自己到自己的距离
  diag(d) <- Inf
  
  # 只看上三角，避免 i-j 和 j-i 重复
  d_upper <- d
  d_upper[lower.tri(d_upper)] <- Inf
  
  min_dist <- min(d_upper, na.rm = TRUE)
  
  pair_pos <- which(d_upper == min_dist, arr.ind = TRUE)
  
  # 如果有并列最近 pair，默认取第一对
  pair_pos_first <- pair_pos[1, ]
  
  pair_idx <- as.integer(pair_pos_first)
  pair_names <- dist_obj$tip_labels[pair_idx]
  
  return(list(
    pair_idx = pair_idx,
    pair_names = pair_names,
    distance = min_dist,
    n_ties = nrow(pair_pos),
    all_tied_pairs = pair_pos
  ))
}

# ============================================================
# TEMP function 2:
# 从固定初始 subset 开始 greedy 添加物种
# 后续比较逻辑与 build_subset_greedy() 保持一致
# ============================================================

build_subset_greedy_from_initial <- function(dist_obj,
                                             subset_size,
                                             initial_subset,
                                             maximize = FALSE,
                                             single_objective = NULL) {
  
  all_tips <- seq_along(dist_obj$tip_labels)
  
  # 允许 initial_subset 是物种名或 index
  if (is.character(initial_subset)) {
    current_subset <- match(initial_subset, dist_obj$tip_labels)
  } else {
    current_subset <- initial_subset
  }
  
  current_subset <- unique(as.integer(current_subset))
  
  if (any(is.na(current_subset))) {
    stop("Some initial species are not found in dist_obj$tip_labels.")
  }
  
  if (length(current_subset) > subset_size) {
    stop("Initial subset is larger than target subset_size.")
  }
  
  if (subset_size > length(all_tips)) {
    stop("subset_size cannot be larger than total number of tips.")
  }
  
  available <- setdiff(all_tips, current_subset)
  
  while (length(current_subset) < subset_size && length(available) > 0) {
    
    best_candidate <- NULL
    best_metrics <- NULL
    best_value <- NULL
    
    for (candidate in available) {
      
      temp_subset <- c(current_subset, candidate)
      temp_metrics <- calc_subset_metrics(dist_obj$dist_mat, temp_subset)
      
      if (is.null(single_objective)) {
        
        if (is.null(best_candidate)) {
          best_candidate <- candidate
          best_metrics <- temp_metrics
        } else {
          if (maximize) {
            if (is_better_lexico_max(temp_metrics, best_metrics)) {
              best_candidate <- candidate
              best_metrics <- temp_metrics
            }
          } else {
            if (is_better_lexico_min(temp_metrics, best_metrics)) {
              best_candidate <- candidate
              best_metrics <- temp_metrics
            }
          }
        }
        
      } else {
        
        temp_value <- temp_metrics[[single_objective]]
        
        if (is.null(best_candidate)) {
          best_candidate <- candidate
          best_value <- temp_value
          best_metrics <- temp_metrics
        } else {
          if (maximize) {
            if (temp_value > best_value) {
              best_candidate <- candidate
              best_value <- temp_value
              best_metrics <- temp_metrics
            }
          } else {
            if (temp_value < best_value) {
              best_candidate <- candidate
              best_value <- temp_value
              best_metrics <- temp_metrics
            }
          }
        }
      }
    }
    
    current_subset <- c(current_subset, best_candidate)
    available <- setdiff(all_tips, current_subset)
  }
  
  final_metrics <- calc_subset_metrics(dist_obj$dist_mat, current_subset)
  
  return(list(
    subset = current_subset,
    subset_names = dist_obj$tip_labels[current_subset],
    metrics = final_metrics,
    initial_subset = initial_subset,
    algorithm = ifelse(
      maximize,
      "greedy_fixed_initial_max",
      "greedy_fixed_initial_min"
    )
  ))
}

# ============================================================
# TEMP function 3:
# 固定最近 pair 起始的 clustered complete algorithm
# ============================================================

run_complete_algorithm_fixed_closest_pair <- function(dist_obj,
                                                      subset_size,
                                                      single_objective = NULL) {
  
  if (subset_size < 2) {
    stop("subset_size must be at least 2 for closest-pair initialization.")
  }
  
  cat("Running fixed-closest-pair clustered algorithm...\n")
  
  # Phase 0: 找全树最近 pair
  cat("  Phase 0: Finding closest pair\n")
  
  closest_pair <- find_closest_pair(dist_obj)
  
  cat("    Closest pair:\n")
  cat("      ", closest_pair$pair_names[1], "\n")
  cat("      ", closest_pair$pair_names[2], "\n")
  cat("    Distance =", closest_pair$distance, "\n")
  cat("    Number of tied closest pairs =", closest_pair$n_ties, "\n")
  
  # Phase 1: 从最近 pair 开始 greedy construction
  cat("  Phase 1: Greedy construction from fixed closest pair\n")
  
  greedy_result <- build_subset_greedy_from_initial(
    dist_obj = dist_obj,
    subset_size = subset_size,
    initial_subset = closest_pair$pair_idx,
    maximize = FALSE,
    single_objective = single_objective
  )
  
  cat("    Greedy result: MinPD =", greedy_result$metrics$MinPD,
      "MeanPD =", greedy_result$metrics$MeanPD,
      "MeanNND =", greedy_result$metrics$MeanNND, "\n")
  
  # Phase 2: exchange refinement，沿用原来的函数
  cat("  Phase 2: Exchange refinement\n")
  
  exchange_result <- refine_subset_exchange(
    dist_obj = dist_obj,
    current_subset = greedy_result$subset,
    maximize = FALSE,
    single_objective = single_objective
  )
  
  cat("    Exchange result: MinPD =", exchange_result$metrics$MinPD,
      "MeanPD =", exchange_result$metrics$MeanPD,
      "MeanNND =", exchange_result$metrics$MeanNND, "\n")
  cat("    Iterations:", exchange_result$iterations, "\n")
  cat("    Improvements:", length(exchange_result$improvements), "\n")
  
  improvement <- list(
    MinPD = greedy_result$metrics$MinPD - exchange_result$metrics$MinPD,
    MeanPD = greedy_result$metrics$MeanPD - exchange_result$metrics$MeanPD,
    MeanNND = greedy_result$metrics$MeanNND - exchange_result$metrics$MeanNND
  )
  
  return(list(
    closest_pair = closest_pair,
    greedy_result = greedy_result,
    exchange_result = exchange_result,
    final_subset = exchange_result$subset,
    final_subset_names = exchange_result$subset_names,
    final_metrics = exchange_result$metrics,
    improvement = improvement,
    algorithm = "fixed_closest_pair_greedy_exchange_min"
  ))
}

# ============================================================
# TEMP analysis:
# 固定最近 pair 起始的 clustered subset
# ============================================================

fixed_clust_result <- run_complete_algorithm_fixed_closest_pair(
  dist_obj = dist_obj,
  subset_size = subset_size,
  single_objective = NULL
)

fixed_clust_subset_idx <- fixed_clust_result$final_subset
fixed_clust_subset_names <- fixed_clust_result$final_subset_names
fixed_clust_metrics <- fixed_clust_result$final_metrics

length(fixed_clust_subset_idx)
fixed_clust_metrics
fixed_clust_subset_names

# ============================================================
# TEMP plot:
# 固定最近 pair 起始 clustered subset
# ============================================================

svg("TEMP_Carnivora_tree_fixed_closest_pair_clustered_s20.svg",
    width = 9, height = 9)

par(mar = c(1, 1, 3, 1))

plot.phylo(carn_tree,
           type = "fan",
           show.tip.label = FALSE,
           no.margin = FALSE,
           cex = 0.3)

# 所有 clustered subset 物种
tiplabels(
  pch = 21,
  bg = "black",
  cex = 0.9,
  tip = which(carn_tree$tip.label %in% fixed_clust_subset_names)
)


title("Fixed closest-pair clustered subset, s = 20",
      line = 1)

dev.off()

# clustered: low-tail p-values
clust_2_MinPD <- calc_p_low(fixed_clust_metrics$MinPD, random_dist_metrics$MinPD)
clust_2_MeanPD <- calc_p_low(fixed_clust_metrics$MeanPD, random_dist_metrics$MeanPD)
clust_2_MeanNND <- calc_p_low(fixed_clust_metrics$MeanNND, random_dist_metrics$MeanNND)
