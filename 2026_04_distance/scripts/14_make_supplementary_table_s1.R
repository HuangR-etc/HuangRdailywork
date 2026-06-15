# 14_make_supplementary_table_s1.R
# Build formal long-format Supplementary Table S1 for selected Cricetidae species.
#
# Usage: Rscript scripts/14_make_supplementary_table_s1.R
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(file_arg) > 0) {
  setwd(normalizePath(file.path(dirname(file_arg[1]), "..")))
}
source("R/01_load_modules.R")
load_project_modules()

cat("=== 14_make_supplementary_table_s1 ===\n")

in_file <- file.path(
  RESULTS_DIR, "sensitivity", "selected_subsets",
  "cricetidae_sensitivity_selected_species.csv"
)
out_dir <- file.path(RESULTS_DIR, "supplementary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file)) {
  stop("Missing selected species file: ", in_file,
       "\nRun scripts/04_run_cricetidae_sensitivity_grid.R first.")
}

selected <- read.csv(in_file, stringsAsFactors = FALSE)
selected <- selected[selected$N %in% CRICETIDAE_POOL_SIZES &
                       selected$s %in% SENSITIVITY_SUBSET_SIZES &
                       selected$s < selected$N, ]

make_long <- function(df, subset_col, subset_type) {
  out <- df[, c("N", "s", "Species_Index", subset_col)]
  names(out)[names(out) == subset_col] <- "Species"
  out <- out[!is.na(out$Species) & out$Species != "", ]
  out$Subset_Type <- subset_type
  out[, c("N", "s", "Subset_Type", "Species_Index", "Species")]
}

s1 <- rbind(
  make_long(selected, "Dispersed", "dispersed"),
  make_long(selected, "Clustered", "clustered")
)
s1 <- s1[order(s1$N, s1$s, s1$Subset_Type, s1$Species_Index), ]

out_file <- file.path(out_dir, "Supplementary_Table_S1_selected_species.csv")
write.csv(s1, out_file, row.names = FALSE)

cat("Saved:", out_file, "\n")
cat("Rows:", nrow(s1), "\n")
