#' Detect communities of co-associated modalities
#'
#' Applies graph community detection to a \code{catmodgraph} object and writes
#' the resulting community membership onto the graph vertices. Communities
#' here are groups of modalities (factor levels) that tend to co-associate
#' across different variables. This is a category co-association analysis,
#' not a respondent-segmentation method; for the latter, see the
#' \pkg{poLCA} or \pkg{FactoMineR} packages.
#'
#' @param x A \code{catmodgraph} object.
#' @param method Character. Community detection method. One of
#'   \code{"louvain"} (default) or \code{"walktrap"}.
#' @param signed Logical. If \code{TRUE}, edges with negative standardised
#'   Pearson residual (repulsion: modalities co-occurring less than
#'   expected under independence) are temporarily dropped before
#'   community detection, so communities are defined solely by positive
#'   co-association. Only supported for \code{method = "louvain"}.
#'   Default \code{FALSE}.
#'
#' @return The input \code{catmodgraph} object with additional vertex attribute
#'   \code{cluster}, and a new component \code{membership} giving the cluster
#'   label for each modality node.
#'
#' @details
#' By default, community detection uses the \code{weight} edge attribute
#' (phi), which is non-negative. This groups modalities by association
#' \emph{magnitude} and does not distinguish attraction from repulsion.
#' Two modalities with large |phi| but opposite signs of the standardised
#' Pearson residual will be pulled together by unsigned Louvain, which
#' can lump opposing substantive patterns into the same community.
#'
#' When \code{signed = TRUE}, edges with negative standardised residual
#' (repulsion) are dropped from a local copy of the graph before Louvain
#' runs. Communities are then defined purely by positive co-association.
#' Repulsion information is \strong{not discarded} — the
#' \code{std_resid} attribute is retained on the original graph, so
#' downstream plotting with \code{plot(mg, signed = TRUE)} still shows
#' attraction (green) and repulsion (red) edges together.
#'
#' This is a pragmatic adaptation that avoids the need for a proper
#' signed community-detection algorithm, which standard Louvain is not
#' (modern \pkg{igraph} rejects negative edge weights in
#' \code{cluster_louvain()}). For principled signed community detection
#' (Traag & Bruggeman, 2009), see the \pkg{signnet} package.
#'
#' \code{signed = TRUE} requires the \code{std_resid} edge attribute,
#' which is attached by \code{\link{build_modality_graph}} and preserved
#' through \code{\link{prune_modality_edges}}.
#'
#' @references
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#'   Fast unfolding of communities in large networks.
#'   \emph{Journal of Statistical Mechanics}, 2008(10), P10008.
#'   \doi{10.1088/1742-5468/2008/10/P10008}
#'
#' Traag, V. A., & Bruggeman, J. (2009). Community detection in networks
#'   with positive and negative links. \emph{Physical Review E}, 80(3),
#'   036115. \doi{10.1103/PhysRevE.80.036115}
#'
#' @examples
#' df <- expand_table(Titanic)
#' mg <- build_modality_graph(df)
#'
#' # Default: unsigned Louvain on phi weights
#' mg <- cluster_modalities(mg)
#' table(mg$membership)
#'
#' # Positive-residual-only Louvain: repulsion edges excluded from clustering
#' mg_pos <- cluster_modalities(mg, signed = TRUE)
#' table(mg_pos$membership)
#'
#' @importFrom igraph cluster_louvain cluster_walktrap membership V E
#'   edge_attr_names delete_edges
#' @export
cluster_modalities <- function(x,
                               method = c("louvain", "walktrap"),
                               signed = FALSE) {
  
  if (!inherits(x, "catmodgraph")) {
    stop("`x` must be a catmodgraph object.", call. = FALSE)
  }
  
  method <- match.arg(method)
  
  if (!is.logical(signed) || length(signed) != 1L || is.na(signed)) {
    stop("`signed` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (signed && method != "louvain") {
    stop("`signed = TRUE` is currently supported only for method = 'louvain'.",
         call. = FALSE)
  }
  
  g <- x$graph
  
  if (igraph::vcount(g) < 2L) {
    stop("Modality graph must contain at least two nodes.", call. = FALSE)
  }
  
  if (igraph::ecount(g) == 0L) {
    stop("Modality graph contains no edges; clustering is not possible.",
         call. = FALSE)
  }
  
  # For signed clustering, drop repulsion edges (std_resid < 0) on a local
  # copy of the graph. The original x$graph is not modified; std_resid is
  # retained on all edges for downstream plot(signed = TRUE).
  g_clust <- g
  if (signed) {
    if (!"std_resid" %in% igraph::edge_attr_names(g_clust)) {
      stop(
        "`signed = TRUE` requires the `std_resid` edge attribute, which is ",
        "attached by build_modality_graph(). This graph appears to have been ",
        "constructed without it.",
        call. = FALSE
      )
    }
    std_res <- igraph::E(g_clust)$std_resid
    
    if (all(is.na(std_res))) {
      stop("All `std_resid` values are NA; signed clustering cannot proceed.",
           call. = FALSE)
    }
    
    drop <- is.na(std_res) | std_res <= 0
    if (any(drop)) {
      g_clust <- igraph::delete_edges(g_clust, which(drop))
    }
    
    if (igraph::ecount(g_clust) == 0L) {
      stop(
        "No positive-residual edges remain after filtering; signed community ",
        "detection cannot proceed. The graph contains only repulsion edges ",
        "or none at all.",
        call. = FALSE
      )
    }
  }
  
  comm <- switch(
    method,
    louvain  = igraph::cluster_louvain(g_clust,
                                       weights = igraph::E(g_clust)$weight),
    walktrap = igraph::cluster_walktrap(g_clust,
                                        weights = igraph::E(g_clust)$weight)
  )
  
  memb <- igraph::membership(comm)
  names(memb) <- igraph::V(g_clust)$name
  
  # Map membership back onto the ORIGINAL graph (which still has all edges)
  # so downstream plotting shows attraction/repulsion edges even in nodes
  # whose community was determined only by positive residuals.
  full_memb <- rep(NA_integer_, igraph::vcount(g))
  names(full_memb) <- igraph::V(g)$name
  full_memb[names(memb)] <- as.integer(memb)
  
  # Nodes that got isolated by the filter but survive in the original graph
  # need a community label; assign each to its own singleton community so
  # downstream code that expects no NA in membership still works.
  if (any(is.na(full_memb))) {
    next_id <- max(full_memb, na.rm = TRUE) + 1L
    miss_idx <- which(is.na(full_memb))
    for (i in miss_idx) {
      full_memb[i] <- next_id
      next_id <- next_id + 1L
    }
  }
  
  igraph::V(g)$cluster <- unname(full_memb)
  x$graph <- g
  x$membership <- full_memb
  x
}