# 10_plotting_trees.R
# Tree plotting functions for figures

#' Plot a tree with highlighted subset tips
#'
#' @param tree Phylo object
#' @param subset_names Character vector of tip names to highlight
#' @param main Title for the plot
#' @param highlight_col Color for highlighted tips
#' @param cex Tip label size
plot_tree_with_subset <- function(tree, subset_names, main = "",
                                  highlight_col = "red", cex = 0.7) {
  tip_colors <- ifelse(tree$tip.label %in% subset_names, highlight_col, "black")
  tip_cex <- ifelse(tree$tip.label %in% subset_names, cex * 1.2, cex)
  
  ape::plot.phylo(tree, main = main, cex = tip_cex,
                  tip.color = tip_colors, label.offset = 0.01)
  
  if (length(subset_names) > 0) {
    legend("topleft", legend = paste("Selected (n =", length(subset_names), ")"),
           text.col = highlight_col, bty = "n")
  }
}

#' Save a tree plot to PDF
#'
#' @param tree Phylo object
#' @param subset_names Character vector of tip names to highlight
#' @param file_path Output PDF path
#' @param main Title
#' @param highlight_col Color
#' @param width PDF width
#' @param height PDF height
save_tree_plot_pdf <- function(tree, subset_names, file_path, main = "",
                               highlight_col = "red",
                               width = 10, height = 8, cex = 0.5) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  pdf(file_path, width = width, height = height)
  plot_tree_with_subset(tree, subset_names, main = main,
                        highlight_col = highlight_col)
  dev.off()
}
