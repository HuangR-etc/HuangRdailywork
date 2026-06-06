# 09_summary_tables.R
# Summary table construction and the main run_empirical_case() runner

# ============================================================
# Summary table helpers
# ============================================================

#' Summarize distance metrics: observed vs random baseline
#'
#' @param pool_label Label for the pool (e.g. "Carnivora_N261", "Cricetidae_C128")
#' @param N Pool size
#' @param subset_size Subset size s
#' @param disp_metrics Extended metrics list for dispersed subset
#' @param clust_metrics Extended metrics list for clustered subset
#' @param random_dist_metrics Data frame of extended metrics for random subsets
#' @return Data frame with columns: Analysis, Pool, N, s, Subset_Type, Metric,
#'         Observed, Baseline_Mean, Baseline_SD, Baseline_Q025, Baseline_Q975,
#'         SES, P_value, Tail, Covariance_Model, Covariance_Param
summarize_distance_observed_vs_random <- function(pool_label, N, subset_size,
                                                  disp_metrics, clust_metrics,
                                                  random_dist_metrics) {
  metrics_names <- c("MinPD", "MeanPD", "MeanNND", "MaxPD")
  
  out_list <- list()
  
  for (stype in c("dispersed", "clustered")) {
    obs <- if (stype == "dispersed") disp_metrics else clust_metrics
    
    for (mname in metrics_names) {
      obs_val <- obs[[mname]]
      null_vals <- random_dist_metrics[[mname]]
      
      ses_val <- calc_ses(obs_val, null_vals)
      
      # Dispersed: upper tail (observed > random)
      # Clustered: lower tail (observed < random)
      if (stype == "dispersed") {
        p_val <- calc_p_high(obs_val, null_vals)
        tail <- "upper"
      } else {
        p_val <- calc_p_low(obs_val, null_vals)
        tail <- "lower"
      }
      
      out_list[[length(out_list) + 1]] <- data.frame(
        Analysis = pool_label,
        Pool = pool_label,
        N = N,
        s = subset_size,
        Subset_Type = stype,
        Metric = mname,
        Observed = obs_val,
        Baseline_Mean = mean(null_vals, na.rm = TRUE),
        Baseline_SD = sd(null_vals, na.rm = TRUE),
        Baseline_Q025 = quantile(null_vals, 0.025, na.rm = TRUE),
        Baseline_Q975 = quantile(null_vals, 0.975, na.rm = TRUE),
        SES = ses_val,
        P_value = p_val,
        Tail = tail,
        Covariance_Model = NA_character_,
        Covariance_Param = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, out_list)
}

#' Summarize dependence diagnostics: observed vs random baseline
#'
#' @param pool_label Label for the pool
#' @param N Pool size
#' @param subset_size Subset size s
#' @param disp_dep Data frame from calc_dependence_from_V for dispersed
#' @param clust_dep Data frame from calc_dependence_from_V for clustered
#' @param random_dep_metrics Data frame from calc_multiple_dependence_from_V
#' @param covariance_model Name of covariance model ("BM", "lambda_BM", "OU")
#' @param covariance_param Parameter value (NA for BM, lambda value, or half_life_frac)
#' @param covariance_param_label Human-readable label for the parameter
#' @return Data frame
summarize_dependence_observed_vs_random <- function(pool_label, N, subset_size,
                                                    disp_dep, clust_dep,
                                                    random_dep_metrics,
                                                    covariance_model,
                                                    covariance_param = NA_real_,
                                                    covariance_param_label = NA_character_) {
  metrics_names <- c("off_mean", "rmax", "neff_mean")
  
  out_list <- list()
  
  for (stype in c("dispersed", "clustered")) {
    obs <- if (stype == "dispersed") disp_dep else clust_dep
    
    for (mname in metrics_names) {
      obs_val <- obs[[mname]]
      null_vals <- random_dep_metrics[[mname]]
      
      ses_val <- calc_ses(obs_val, null_vals)
      
      # Dispersed: lower tail for off_mean/rmax (less dependence), upper for neff_mean
      # Clustered: upper tail for off_mean/rmax (more dependence), lower for neff_mean
      if (mname %in% c("off_mean", "rmax")) {
        if (stype == "dispersed") {
          p_val <- calc_p_low(obs_val, null_vals)
          tail <- "lower"
        } else {
          p_val <- calc_p_high(obs_val, null_vals)
          tail <- "upper"
        }
      } else {
        # neff_mean: dispersed should have higher neff, clustered lower
        if (stype == "dispersed") {
          p_val <- calc_p_high(obs_val, null_vals)
          tail <- "upper"
        } else {
          p_val <- calc_p_low(obs_val, null_vals)
          tail <- "lower"
        }
      }
      
      out_list[[length(out_list) + 1]] <- data.frame(
        Analysis = pool_label,
        Pool = pool_label,
        N = N,
        s = subset_size,
        Subset_Type = stype,
        Metric = mname,
        Observed = obs_val,
        Baseline_Mean = mean(null_vals, na.rm = TRUE),
        Baseline_SD = sd(null_vals, na.rm = TRUE),
        Baseline_Q025 = quantile(null_vals, 0.025, na.rm = TRUE),
        Baseline_Q975 = quantile(null_vals, 0.975, na.rm = TRUE),
        SES = ses_val,
        P_value = p_val,
        Tail = tail,
        Covariance_Model = covariance_model,
        Covariance_Param = covariance_param,
        Covariance_Param_Label = covariance_param_label,
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, out_list)
}

#' Make selected species table
#'
#' @param pool_label Label for the pool
#' @param N Pool size
#' @param subset_size Subset size s
#' @param disp_names Character vector of dispersed subset tip names
#' @param clust_names Character vector of clustered subset tip names
#' @return Data frame
make_selected_species_table <- function(pool_label, N, subset_size,
                                        disp_names, clust_names) {
  max_len <- max(length(disp_names), length(clust_names))
  disp_names_padded <- c(disp_names, rep(NA, max_len - length(disp_names)))
  clust_names_padded <- c(clust_names, rep(NA, max_len - length(clust_names)))
  
  data.frame(
    Analysis = pool_label,
    N = N,
    s = subset_size,
    Species_Index = seq_len(max_len),
    Dispersed = disp_names_padded,
    Clustered = clust_names_padded,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Main empirical case runner
# ============================================================

#' Run a complete empirical analysis case
#'
#' This is the main workhorse function. It:
#' 1. Creates distance object
#' 2. Selects dispersed and clustered subsets
#' 3. Samples random baseline
#' 4. Calculates distance metrics for all
#' 5. Calculates BM dependence diagnostics
#' 6. Saves result as RDS
#'
#' @param pool_tree Phylo object for the pool
#' @param pool_label Label (e.g. "Carnivora_N261", "Cricetidae_C128")
#' @param subset_size Desired subset size s
#' @param n_null_reps Number of random subsets for baseline
#' @param seed Random seed
#' @param out_dir Output directory for RDS
#' @param save_raw If TRUE, save RDS
#' @param overwrite If TRUE, overwrite existing RDS
#' @return List with all results
run_empirical_case <- function(pool_tree,
                               pool_label,
                               subset_size,
                               n_null_reps = 1000,
                               seed = GLOBAL_SEED,
                               out_dir,
                               save_raw = TRUE,
                               overwrite = FALSE) {
  
  N <- ape::Ntip(pool_tree)
  case_id <- paste0(pool_label, "_s", subset_size)
  rds_file <- file.path(out_dir, paste0(case_id, ".rds"))
  
  if (file.exists(rds_file) && !overwrite) {
    message("Loading existing result: ", rds_file)
    return(readRDS(rds_file))
  }
  
  if (subset_size >= N) {
    stop("subset_size must be smaller than N.")
  }
  
  dist_obj <- create_distance_object(pool_tree)
  
  # Selected subsets
  disp_result <- run_dispersed_algorithm(dist_obj, subset_size)
  clust_result <- select_clustered_greedy_exchange(
    dist_obj = dist_obj,
    subset_size = subset_size,
    max_iterations = CLUSTERED_MAX_EXCHANGE_ITERATIONS,
    tol = CLUSTERED_EXCHANGE_TOL,
    verbose = FALSE
  )
  
  # Random baseline
  set.seed(seed)
  random_idx <- sample_random_subsets(
    dist_obj = dist_obj,
    subset_size = subset_size,
    n_reps = n_null_reps,
    replace = FALSE
  )
  
  random_names <- lapply(random_idx, function(x) dist_obj$tip_labels[x])
  
  # Distance metrics
  random_dist <- calc_multiple_subsets_metrics_extended(
    dist_obj$dist_mat,
    random_idx
  )
  
  disp_dist <- calc_subset_metrics_extended(
    dist_obj$dist_mat,
    disp_result$final_subset
  )
  
  clust_dist <- calc_subset_metrics_extended(
    dist_obj$dist_mat,
    clust_result$final_subset
  )
  
  distance_summary <- summarize_distance_observed_vs_random(
    pool_label = pool_label,
    N = N,
    subset_size = subset_size,
    disp_metrics = disp_dist,
    clust_metrics = clust_dist,
    random_dist_metrics = random_dist
  )
  
  # BM dependence
  V_bm <- make_bm_covariance(pool_tree)
  
  disp_dep <- calc_dependence_from_V(V_bm, disp_result$final_subset_names)
  clust_dep <- calc_dependence_from_V(V_bm, clust_result$final_subset_names)
  random_dep <- calc_multiple_dependence_from_V(V_bm, random_names)
  
  dependence_summary <- summarize_dependence_observed_vs_random(
    pool_label = pool_label,
    N = N,
    subset_size = subset_size,
    disp_dep = disp_dep,
    clust_dep = clust_dep,
    random_dep_metrics = random_dep,
    covariance_model = "BM",
    covariance_param = NA_real_,
    covariance_param_label = "BM"
  )
  
  selected_species <- make_selected_species_table(
    pool_label = pool_label,
    N = N,
    subset_size = subset_size,
    disp_names = disp_result$final_subset_names,
    clust_names = clust_result$final_subset_names
  )
  
  result <- list(
    case_id = case_id,
    pool_label = pool_label,
    N = N,
    s = subset_size,
    seed = seed,
    n_null_reps = n_null_reps,
    pool_tree = pool_tree,
    dispersed = disp_result,
    clustered = clust_result,
    random_idx = random_idx,
    random_names = random_names,
    random_dist = random_dist,
    dispersed_distance = disp_dist,
    clustered_distance = clust_dist,
    random_dep_BM = random_dep,
    distance_summary = distance_summary,
    dependence_summary_BM = dependence_summary,
    selected_species = selected_species
  )
  
  if (save_raw) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, rds_file)
  }
  
  result
}
