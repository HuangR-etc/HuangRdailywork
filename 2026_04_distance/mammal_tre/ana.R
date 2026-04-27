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



#临时测试

# ============================================================
# TEMPORARY TEST:
# 多个随机种子重复构造 clustered subset
# 并与固定 random baseline 比较 p 值
# ============================================================

# 这部分是临时测试代码，正式分析完成后可以删除

# ------------------------------------------------------------
# 0. 确认前面已经存在这些对象
# ------------------------------------------------------------

exists("dist_obj")
exists("subset_size")
exists("n_greedy_starts")
exists("random_dist_metrics")

dim(random_dist_metrics)
head(random_dist_metrics)


# ------------------------------------------------------------
# 2. 设置测试参数
# ------------------------------------------------------------

# 你可以改成 1:50、1:100，或者自定义一组种子
cluster_test_seeds <- 1:100

n_cluster_tests <- length(cluster_test_seeds)

alpha <- 0.05

# 用于保存每次 clustered 结果的总表
cluster_test_summary <- data.frame()

# 用于保存每次选中的物种
cluster_test_species_list <- vector("list", n_cluster_tests)

# ------------------------------------------------------------
# 3. 多个 seed 重复运行 clustered subset construction
# ------------------------------------------------------------

for (i in seq_along(cluster_test_seeds)) {
  
  seed_now <- cluster_test_seeds[i]
  
  cat("\n========================================\n")
  cat("Clustered test", i, "of", n_cluster_tests, "\n")
  cat("Seed:", seed_now, "\n")
  cat("========================================\n")
  
  set.seed(seed_now)
  
  test_clust_result <- run_complete_algorithm(
    dist_obj = dist_obj,
    subset_size = subset_size,
    maximize = FALSE,
    n_greedy_starts = n_greedy_starts
  )
  
  test_clust_idx <- test_clust_result$final_subset
  test_clust_names <- test_clust_result$final_subset_names
  
  # 重新计算一次 metrics，避免直接依赖 result 内部字段
  test_clust_metrics <- calc_subset_metrics(
    dist_mat = dist_obj$dist_mat,
    subset = test_clust_idx
  )
  
  # clustered subset 是 lower-tail comparison:
  # p = P(random <= observed)
  p_MinPD <- calc_p_low(
    obs = test_clust_metrics$MinPD,
    null_values = random_dist_metrics$MinPD
  )
  
  p_MeanPD <- calc_p_low(
    obs = test_clust_metrics$MeanPD,
    null_values = random_dist_metrics$MeanPD
  )
  
  p_MeanNND <- calc_p_low(
    obs = test_clust_metrics$MeanNND,
    null_values = random_dist_metrics$MeanNND
  )
  
  # 三个距离指标是否都显著
  all_three_significant <- all(
    c(p_MinPD, p_MeanPD, p_MeanNND) < alpha
  )
  
  # 至少一个不显著
  any_not_significant <- any(
    c(p_MinPD, p_MeanPD, p_MeanNND) >= alpha
  )
  
  # 记录该次结果
  temp_summary <- data.frame(
    Test_ID = i,
    Seed = seed_now,
    Subset_Size = length(test_clust_idx),
    
    MinPD = test_clust_metrics$MinPD,
    MeanPD = test_clust_metrics$MeanPD,
    MeanNND = test_clust_metrics$MeanNND,
    
    P_MinPD_low = p_MinPD,
    P_MeanPD_low = p_MeanPD,
    P_MeanNND_low = p_MeanNND,
    
    Sig_MinPD = p_MinPD < alpha,
    Sig_MeanPD = p_MeanPD < alpha,
    Sig_MeanNND = p_MeanNND < alpha,
    
    All_Three_Significant = all_three_significant,
    Any_Not_Significant = any_not_significant,
    
    stringsAsFactors = FALSE
  )
  
  cluster_test_summary <- rbind(cluster_test_summary, temp_summary)
  
  # 保存该次选中的物种
  cluster_test_species_list[[i]] <- data.frame(
    Test_ID = i,
    Seed = seed_now,
    Rank_In_Subset = seq_along(test_clust_names),
    tiplabel = test_clust_names,
    stringsAsFactors = FALSE
  )
  
  # 简单打印当前结果，方便你边跑边看
  print(temp_summary)
}
# ============================================================
# 4. 查看多次 clustered 测试结果
# ============================================================

cluster_test_summary

summary(cluster_test_summary[, c(
  "MinPD", "MeanPD", "MeanNND",
  "P_MinPD_low", "P_MeanPD_low", "P_MeanNND_low"
)])

# 三个 p 值是否每次都显著
table(cluster_test_summary$All_Three_Significant)

# 哪些 seed 没有三个都显著
cluster_test_summary[
  cluster_test_summary$All_Three_Significant == FALSE,
]

# 分别看每个指标的显著次数
colSums(cluster_test_summary[, c(
  "Sig_MinPD",
  "Sig_MeanPD",
  "Sig_MeanNND"
)])

# p 值最大的是哪几次
cluster_test_summary[
  order(cluster_test_summary$P_MinPD_low, decreasing = TRUE),
][1:10, ]

cluster_test_summary[
  order(cluster_test_summary$P_MeanPD_low, decreasing = TRUE),
][1:10, ]

cluster_test_summary[
  order(cluster_test_summary$P_MeanNND_low, decreasing = TRUE),
][1:10, ]