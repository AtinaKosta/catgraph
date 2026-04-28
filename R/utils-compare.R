#' @keywords internal
#' @noRd
# -----------------------------------------------------------------------------
# Shared helpers used by compare_catgraphs() and compare_modality_graphs().
# Kept deliberately minimal: only the pieces that are genuinely identical
# across the two functions are extracted here. Aesthetics, pruning pipelines,
# and vertex-set policies remain in the caller because they legitimately
# differ.
# -----------------------------------------------------------------------------

# Validate a named list of objects of a single S3 class.
# Stops with a tailored message; returns invisible(TRUE) on success.
.validate_compare_input <- function(x, class_name,
                                    min_len = 2L,
                                    arg_name = "x") {
  # English-number strings for the error message. The existing test suite
  # matches on literal substrings like "at least two", so we preserve the
  # exact wording for the common cases (2, 3) and fall back to the numeral.
  len_word <- switch(as.character(min_len),
                     "1" = "one", "2" = "two", "3" = "three",
                     as.character(min_len))
  if (!is.list(x) || length(x) < min_len) {
    stop(sprintf("`%s` must be a list of at least %s %s objects.",
                 arg_name, len_word, class_name),
         call. = FALSE)
  }
  if (is.null(names(x)) || any(names(x) == "") || anyDuplicated(names(x))) {
    stop(sprintf("`%s` must be a named list with unique, non-empty names.",
                 arg_name),
         call. = FALSE)
  }
  if (!all(vapply(x, inherits, logical(1L), class_name))) {
    stop(sprintf("All elements of `%s` must be %s objects.",
                 arg_name, class_name),
         call. = FALSE)
  }
  invisible(TRUE)
}


# Build a union igraph from a list of igraph objects that share the same
# vertex namespace. Isolated vertices from each input are preserved; edges
# from later graphs are added if both endpoints already exist in the running
# union. Caller is responsible for ensuring the vertex sets are reconciled
# upstream (compare_modality_graphs does this via `restrict`).
.union_graph <- function(graphs) {
  if (length(graphs) == 0L) {
    stop(".union_graph(): empty input.", call. = FALSE)
  }
  u <- graphs[[1L]]
  if (length(graphs) == 1L) return(igraph::simplify(u))
  for (i in seq.int(2L, length(graphs))) {
    el_i <- igraph::as_edgelist(graphs[[i]])
    if (nrow(el_i) == 0L) next
    mask <- el_i[, 1L] %in% igraph::V(u)$name &
            el_i[, 2L] %in% igraph::V(u)$name
    if (any(mask)) {
      u <- igraph::add_edges(u, t(el_i[mask, , drop = FALSE]))
    }
  }
  igraph::simplify(u)
}


# Compute a layout on the union graph and return a name-indexed matrix
# suitable for realigning per-panel subgraphs. If `rescale_to` is supplied
# (two-element numeric), coordinates are rescaled into that range per axis.
.shared_layout <- function(union_g, layout_fn, rescale_to = NULL) {
  lay <- layout_fn(union_g)
  rownames(lay) <- igraph::V(union_g)$name
  if (!is.null(rescale_to)) {
    if (!is.numeric(rescale_to) || length(rescale_to) != 2L) {
      stop(".shared_layout(): `rescale_to` must be a length-2 numeric.",
           call. = FALSE)
    }
    lay[, 1L] <- scales::rescale(lay[, 1L], to = rescale_to)
    lay[, 2L] <- scales::rescale(lay[, 2L], to = rescale_to)
  }
  lay
}


# Return the sub-matrix of `lay_ref` in the order of V(g)$name. Used inside
# the panel loop to align each panel's graph to the shared layout.
.align_layout <- function(g, lay_ref) {
  lay_ref[igraph::V(g)$name, , drop = FALSE]
}
