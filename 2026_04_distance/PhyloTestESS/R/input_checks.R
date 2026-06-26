#' Input validation and tree preprocessing
#'
#' Lightweight helpers for checking tree inputs and building patristic
#' distance matrices for candidate pools.
NULL

#' @export
check_phylo_input <- function(tree, candidates = NULL, stop_on_error = TRUE) {
  issues <- character()

  if (!inherits(tree, "phylo")) {
    issues <- c(issues, "tree must be a phylo object")
  }

  if (inherits(tree, "phylo")) {
    if (is.null(tree$edge.length)) {
      issues <- c(issues, "tree must contain branch lengths")
    }
    if (anyDuplicated(tree$tip.label) > 0) {
      issues <- c(issues, "tree contains duplicate tip labels")
    }
    if (anyNA(tree$tip.label) || any(tree$tip.label == "")) {
      issues <- c(issues, "tree contains missing or empty tip labels")
    }
  }

  if (!is.null(candidates)) {
    if (!is.character(candidates) || length(candidates) == 0) {
      issues <- c(issues, "candidates must be a non-empty character vector")
    }
    if (anyDuplicated(candidates) > 0) {
      issues <- c(issues, "candidates contains duplicate species names")
    }
    if (inherits(tree, "phylo")) {
      missing_candidates <- setdiff(candidates, tree$tip.label)
      if (length(missing_candidates) > 0) {
        issues <- c(
          issues,
          paste0(length(missing_candidates), " candidate species not found in tree")
        )
      }
    }
  }

  if (stop_on_error && length(issues) > 0) {
    stop(paste(issues, collapse = "; "), call. = FALSE)
  }

  if (!stop_on_error) {
    return(list(valid = length(issues) == 0, issues = issues))
  }

  invisible(TRUE)
}

#' @export
match_species <- function(tree, species, fuzzy = FALSE) {
  check_phylo_input(tree)

  matched <- character()
  unmatched <- character()
  duplicated <- unique(species[duplicated(species)])

  for (sp in unique(species)) {
    if (sp %in% tree$tip.label) {
      matched <- c(matched, sp)
    } else if (fuzzy) {
      alt <- if (grepl("_", sp)) gsub("_", " ", sp) else gsub(" ", "_", sp)
      if (alt %in% tree$tip.label) {
        matched <- c(matched, alt)
      } else {
        unmatched <- c(unmatched, sp)
      }
    } else {
      unmatched <- c(unmatched, sp)
    }
  }

  list(
    matched = matched,
    unmatched = unmatched,
    duplicated = duplicated
  )
}

#' @export
prune_to_candidates <- function(tree, candidates) {
  check_phylo_input(tree, candidates)

  matched <- intersect(candidates, tree$tip.label)
  if (length(matched) == 0) {
    stop("No candidate species found in tree.", call. = FALSE)
  }

  list(
    pruned_tree = ape::drop.tip(tree, setdiff(tree$tip.label, matched)),
    removed = setdiff(tree$tip.label, matched),
    unmatched = setdiff(candidates, tree$tip.label)
  )
}

#' @export
patristic_matrix <- function(tree, tips = NULL) {
  check_phylo_input(tree)
  d <- ape::cophenetic.phylo(tree)

  if (!is.null(tips)) {
    missing_tips <- setdiff(tips, rownames(d))
    if (length(missing_tips) > 0) {
      stop(
        paste0(length(missing_tips), " specified tips not found in the tree."),
        call. = FALSE
      )
    }
    d <- d[tips, tips, drop = FALSE]
  }

  d
}
