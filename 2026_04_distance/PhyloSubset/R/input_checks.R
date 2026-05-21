#' Input validation and preprocessing functions
#'
#' Functions for checking phylogenetic input data, matching species names,
#' pruning trees to candidate pools, and computing patristic distance matrices.
#'
#' @name input_checks
NULL

#' @rdname input_checks
#' @export
#' @title Check phylogenetic input validity
#' @description Validates a phylogenetic tree and optional candidate species list.
#'   Checks that \code{tree} is a \code{phylo} object, has branch lengths,
#'   has no duplicate tip labels, and that all candidate species are present
#'   in the tree.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param candidates Optional character vector of candidate species names.
#' @param stop_on_error If \code{TRUE}, stops on first error; otherwise returns
#'   a list of issues.
#' @return If \code{stop_on_error = TRUE}, returns \code{TRUE} invisibly on success.
#'   If \code{stop_on_error = FALSE}, returns a list with components \code{valid}
#'   (logical) and \code{issues} (character vector).
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(10)
#' check_phylo_input(tree)
#' }
check_phylo_input <- function(tree, candidates = NULL, stop_on_error = TRUE) {
  issues <- character()

  if (!inherits(tree, "phylo")) {
    issues <- c(issues, "tree must be a phylo object")
  }

  if (inherits(tree, "phylo")) {
    if (is.null(tree$edge.length)) {
      issues <- c(issues, "tree has no branch lengths")
    }

    if (anyDuplicated(tree$tip.label) > 0) {
      issues <- c(issues, "tree has duplicate tip labels")
    }

    if (any(is.na(tree$tip.label)) || any(tree$tip.label == "")) {
      issues <- c(issues, "tree has NA or empty tip labels")
    }
  }

  if (!is.null(candidates)) {
    if (any(is.na(candidates)) || any(candidates == "")) {
      issues <- c(issues, "candidates contains NA or empty strings")
    }

    if (anyDuplicated(candidates) > 0) {
      issues <- c(issues, "candidates contains duplicate species names")
    }

    if (inherits(tree, "phylo")) {
      missing_candidates <- setdiff(candidates, tree$tip.label)
      if (length(missing_candidates) > 0) {
        issues <- c(issues, paste0(
          length(missing_candidates),
          " candidate species not found in tree"
        ))
      }
    }
  }

  if (stop_on_error && length(issues) > 0) {
    stop(paste(issues, collapse = "; "))
  }

  if (!stop_on_error) {
    return(list(valid = length(issues) == 0, issues = issues))
  }

  invisible(TRUE)
}

#' @rdname input_checks
#' @export
#' @title Match species names between tree and query
#' @description Matches a vector of species names against tree tip labels.
#'   Optionally handles underscore/space conversions.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param species Character vector of species names to match.
#' @param fuzzy If \code{TRUE}, also try replacing underscores with spaces
#'   and vice versa when exact matches fail.
#' @return A list with components:
#'   \item{matched}{Character vector of successfully matched species.}
#'   \item{unmatched}{Character vector of species not found in tree.}
#'   \item{duplicated}{Character vector of species appearing multiple times
#'     in the input.}
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(5)
#' tree$tip.label <- gsub("t", "Species_", tree$tip.label)
#' match_species(tree, c("Species_1", "Species_X"))
#' }
match_species <- function(tree, species, fuzzy = FALSE) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }

  tips <- tree$tip.label
  dup_input <- species[duplicated(species)]

  matched <- character()
  unmatched <- character()

  for (sp in unique(species)) {
    if (sp %in% tips) {
      matched <- c(matched, sp)
    } else if (fuzzy) {
      sp_alt <- if (grepl("_", sp)) gsub("_", " ", sp) else gsub(" ", "_", sp)
      if (sp_alt %in% tips) {
        matched <- c(matched, sp_alt)
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
    duplicated = unique(dup_input)
  )
}

#' @rdname input_checks
#' @export
#' @title Prune tree to candidate species
#' @description Prunes a phylogenetic tree to retain only the specified
#'   candidate species. Returns the pruned tree and lists of removed and
#'   unmatched species.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param candidates Character vector of candidate species to retain.
#' @return A list with components:
#'   \item{pruned_tree}{The pruned \code{phylo} object.}
#'   \item{removed}{Character vector of species removed from the tree.}
#'   \item{unmatched}{Character vector of candidates not found in the tree.}
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(10)
#' candidates <- tree$tip.label[1:5]
#' result <- prune_to_candidates(tree, candidates)
#' }
prune_to_candidates <- function(tree, candidates) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }

  matched <- intersect(candidates, tree$tip.label)
  unmatched <- setdiff(candidates, tree$tip.label)
  removed <- setdiff(tree$tip.label, matched)

  if (length(matched) == 0) {
    stop("No candidate species found in tree")
  }

  pruned <- ape::drop.tip(tree, setdiff(tree$tip.label, matched))

  list(
    pruned_tree = pruned,
    removed = removed,
    unmatched = unmatched
  )
}

#' @rdname input_checks
#' @export
#' @title Compute patristic distance matrix
#' @description Computes the patristic (phylogenetic) distance matrix for a
#'   tree, optionally restricted to a subset of tips.
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param tips Optional character vector of tip names to include. If NULL,
#'   all tips are used.
#' @return A numeric distance matrix with tip names as row and column names.
#' @examples
#' \donttest{
#' library(ape)
#' tree <- rtree(10)
#' D <- patristic_matrix(tree)
#' }
patristic_matrix <- function(tree, tips = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }

  D <- ape::cophenetic.phylo(tree)

  if (!is.null(tips)) {
    missing <- setdiff(tips, rownames(D))
    if (length(missing) > 0) {
      stop(paste0(
        length(missing),
        " specified tips not found in tree"
      ))
    }
    D <- D[tips, tips, drop = FALSE]
  }

  D
}
