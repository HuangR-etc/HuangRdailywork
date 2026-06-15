# scripts/08_diagnose_s64_random_maxoffcor_pairs.R
# Diagnose random MaxOffCor distribution and the species pair producing rmax
# for Cricetidae s = 64 cases.

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_file[grepl("^--file=", args_file)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}

source("R/01_load_modules.R")
load_project_modules()

cat("=== 08_diagnose_s64_random_maxoffcor_pairs ===\n")

raw_dir <- file.path(RESULTS_DIR, "sensitivity", "raw_cases")
out_dir <- file.path(RESULTS_DIR, "sensitivity", "diagnostics_s64_maxoffcor")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# s = 64 is only feasible when N > 64
N_values <- CRICETIDAE_POOL_SIZES[CRICETIDAE_POOL_SIZES > 64]
s_i <- 64

# ------------------------------------------------------------
# Helper: normalize random subset object to a list of character vectors
# ------------------------------------------------------------
as_subset_list <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    return(x)
  }
  
  if (is.matrix(x)) {
    return(split(x, row(x)))
  }
  
  if (is.data.frame(x)) {
    return(split(as.matrix(x), row(as.matrix(x))))
  }
  
  stop("Unsupported random_names structure: ", paste(class(x), collapse = ", "))
}

# ------------------------------------------------------------
# Helper: find MaxOffCor and the species pair producing it
# ------------------------------------------------------------
find_maxoffcor_pair <- function(V, subset_names) {
  subset_names <- as.character(subset_names)
  subset_names <- subset_names[!is.na(subset_names)]
  
  R <- cov2cor(V)
  R_sub <- R[subset_names, subset_names, drop = FALSE]
  
  # Use upper triangle only; exclude diagonal
  upper_idx <- upper.tri(R_sub)
  vals <- R_sub[upper_idx]
  
  max_val <- max(vals, na.rm = TRUE)
  
  # If there are ties, keep the first pair but also report number of ties
  pair_idx_all <- which(upper_idx & R_sub == max_val, arr.ind = TRUE)
  pair_idx <- pair_idx_all[1, ]
  
  species_1 <- rownames(R_sub)[pair_idx[1]]
  species_2 <- colnames(R_sub)[pair_idx[2]]
  
  data.frame(
    MaxOffCor = max_val,
    Species_1 = species_1,
    Species_2 = species_2,
    N_tied_pairs = nrow(pair_idx_all),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------
all_random_pairs <- list()
all_observed_pairs <- list()
all_threshold_summaries <- list()

for (N_i in N_values) {
  pool_label <- paste0("Cricetidae_C", N_i)
  case_id <- paste0(pool_label, "_s", s_i)
  rds_file <- file.path(raw_dir, paste0(case_id, ".rds"))
  
  cat("Processing:", case_id, "\n")
  
  if (!file.exists(rds_file)) {
    cat("  WARNING: missing raw RDS:", rds_file, "\n")
    next
  }
  
  result <- readRDS(rds_file)
  
  pool_tree <- result$pool_tree
  V_bm <- make_bm_covariance(pool_tree)
  
  random_list <- as_subset_list(result$random_names)
  
  # ---- Random subsets: MaxOffCor distribution + max pair per replicate ----
  random_pair_df <- do.call(rbind, lapply(seq_along(random_list), function(rep_i) {
    one <- find_maxoffcor_pair(V_bm, random_list[[rep_i]])
    
    data.frame(
      N = N_i,
      s = s_i,
      Replicate = rep_i,
      one,
      stringsAsFactors = FALSE
    )
  }))
  
  all_random_pairs[[length(all_random_pairs) + 1]] <- random_pair_df
  
  # ---- Observed dispersed and clustered subsets, for comparison ----
  disp_pair <- find_maxoffcor_pair(
    V_bm,
    result$dispersed$final_subset_names
  )
  clust_pair <- find_maxoffcor_pair(
    V_bm,
    result$clustered$final_subset_names
  )
  
  observed_pair_df <- rbind(
    data.frame(
      N = N_i,
      s = s_i,
      Subset_Type = "dispersed",
      disp_pair,
      stringsAsFactors = FALSE
    ),
    data.frame(
      N = N_i,
      s = s_i,
      Subset_Type = "clustered",
      clust_pair,
      stringsAsFactors = FALSE
    )
  )
  
  all_observed_pairs[[length(all_observed_pairs) + 1]] <- observed_pair_df
  
  # ---- Optional: summarize how extreme the observed clustered value is ----
  random_rmax <- random_pair_df$MaxOffCor
  clust_rmax <- clust_pair$MaxOffCor[1]
  
  threshold_summary <- data.frame(
    N = N_i,
    s = s_i,
    Random_mean = mean(random_rmax, na.rm = TRUE),
    Random_sd = sd(random_rmax, na.rm = TRUE),
    Random_min = min(random_rmax, na.rm = TRUE),
    Random_q50 = as.numeric(quantile(random_rmax, 0.50, na.rm = TRUE)),
    Random_q95 = as.numeric(quantile(random_rmax, 0.95, na.rm = TRUE)),
    Random_q99 = as.numeric(quantile(random_rmax, 0.99, na.rm = TRUE)),
    Random_max = max(random_rmax, na.rm = TRUE),
    Clustered_observed = clust_rmax,
    Clustered_percentile_vs_random = mean(random_rmax <= clust_rmax, na.rm = TRUE),
    P_one_sided_clustered_greater = mean(random_rmax >= clust_rmax, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  all_threshold_summaries[[length(all_threshold_summaries) + 1]] <- threshold_summary
  
  # ---- Save one histogram per N ----
  pdf(
    file.path(out_dir, paste0("Random_MaxOffCor_distribution_N", N_i, "_s64.pdf")),
    width = 6,
    height = 4
  )
  
  hist(
    random_rmax,
    breaks = 30,
    main = paste0("Random MaxOffCor distribution: N = ", N_i, ", s = 64"),
    xlab = "Random subset MaxOffCor"
  )
  
  abline(v = clust_rmax, lwd = 2, lty = 2)
  abline(v = disp_pair$MaxOffCor[1], lwd = 2, lty = 3)
  
  legend(
    "topleft",
    legend = c("Clustered observed", "Dispersed observed"),
    lty = c(2, 3),
    lwd = 2,
    bty = "n"
  )
  
  dev.off()
}

# ------------------------------------------------------------
# Save combined outputs
# ------------------------------------------------------------
random_pairs_combined <- do.call(rbind, all_random_pairs)
observed_pairs_combined <- do.call(rbind, all_observed_pairs)
threshold_summaries_combined <- do.call(rbind, all_threshold_summaries)

write.csv(
  random_pairs_combined,
  file.path(out_dir, "s64_random_MaxOffCor_pairs_by_replicate.csv"),
  row.names = FALSE
)

write.csv(
  observed_pairs_combined,
  file.path(out_dir, "s64_observed_MaxOffCor_pairs.csv"),
  row.names = FALSE
)

write.csv(
  threshold_summaries_combined,
  file.path(out_dir, "s64_random_MaxOffCor_distribution_summary.csv"),
  row.names = FALSE
)

cat("Done. Outputs saved to:\n")
cat(out_dir, "\n")