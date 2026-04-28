#' Prune edges from a modality graph
#'
#' Removes weak or non-significant edges from a \code{catmodgraph} object.
#' This is the modality-level analogue of \code{\link{prune_edges}()} for
#' variable graphs.
#'
#' @param x A \code{catmodgraph} object.
#' @param min_weight Numeric, non-negative. Edges with phi strictly below this
#'   value are removed. Since phi lies in [0, 1], meaningful thresholds are
#'   usually in that range. Default is \code{0}.
#' @param max_p Numeric in [0, 1]. Edges with p-value strictly above this
#'   value are removed. Default is \code{1}.
#' @param remove_isolates Logical. If \code{TRUE}, vertices with degree 0 after
#'   pruning are removed. Default is \code{FALSE}.
#'
#' @return A pruned \code{catmodgraph} object.
#'
#' @examples
#' df <- expand_table(Titanic)
#' mg <- build_modality_graph(df)
#' mg2 <- prune_modality_edges(mg, min_weight = 0.1, max_p = 0.05)
#' mg2
#'
#' @importFrom igraph delete_edges delete_vertices E degree ecount vcount edge_attr_names
#' @export
prune_modality_edges <- function(x,
                                 min_weight = 0,
                                 max_p = 1,
                                 remove_isolates = FALSE) {
  
  if (!inherits(x, "catmodgraph")) {
    stop("`x` must be a catmodgraph object.", call. = FALSE)
  }
  
  if (!is.numeric(min_weight) || length(min_weight) != 1L || is.na(min_weight) ||
      min_weight < 0) {
    stop("`min_weight` must be a single non-negative number.", call. = FALSE)
  }
  
  if (!is.numeric(max_p) || length(max_p) != 1L || is.na(max_p) ||
      max_p < 0 || max_p > 1) {
    stop("`max_p` must be a single number in [0, 1].", call. = FALSE)
  }
  
  if (!is.logical(remove_isolates) || length(remove_isolates) != 1L || is.na(remove_isolates)) {
    stop("`remove_isolates` must be TRUE or FALSE.", call. = FALSE)
  }
  
  g <- x$graph
  
  if (igraph::ecount(g) == 0L) {
    return(x)
  }
  
  w <- igraph::E(g)$weight
  p <- igraph::E(g)$p_value
  
  keep <- !is.na(w) & !is.na(p) & (w >= min_weight) & (p <= max_p)
  
  if (!all(keep)) {
    g <- igraph::delete_edges(g, igraph::E(g)[!keep])
  }
  
  if (remove_isolates && igraph::vcount(g) > 0L) {
    deg <- igraph::degree(g)
    if (any(deg == 0L)) {
      g <- igraph::delete_vertices(g, which(deg == 0L))
    }
  }
  
  # Keep metadata aligned with surviving vertices
  verts <- if (igraph::vcount(g) > 0L) igraph::V(g)$name else character(0)
  
  x$graph <- g
  x$modalities <- x$modalities[match(verts, x$modalities$node, nomatch = 0L), , drop = FALSE]
  
  # Keep only indicator columns corresponding to surviving modalities
  if (length(verts) > 0L) {
    keep_cols <- colnames(x$indicator_matrix) %in% verts
    x$indicator_matrix <- x$indicator_matrix[, keep_cols, drop = FALSE]
  } else {
    x$indicator_matrix <- x$indicator_matrix[, FALSE, drop = FALSE]
  }
  
  # Drop membership if it no longer matches the graph
  if (!is.null(x$membership)) {
    x$membership <- x$membership[verts]
    if (length(x$membership) == 0L) {
      x$membership <- NULL
    } else {
      igraph::V(x$graph)$cluster <- unname(x$membership)
    }
  }
  
  x
}