#' Detect variable communities in a catgraph using graph clustering algorithms
#'
#' Applies a community detection algorithm to the weighted undirected graph
#' of a \code{catgraph} object. Assigns cluster membership as a vertex
#' attribute and returns the updated object. The Louvain algorithm
#' (Blondel et al., 2008) is the default as it is fast, handles weighted
#' edges, and requires no pre-specification of the number of clusters.
#'
#' @param x A \code{catgraph} object.
#' @param method Character. Community detection algorithm to use. One of
#'   \code{"louvain"} (default), \code{"walktrap"}, \code{"fast_greedy"},
#'   \code{"label_prop"}, or \code{"edge_betweenness"}. All methods use
#'   edge weights when available.
#' @param resolution Numeric. Resolution parameter for the Louvain method
#'   (\code{igraph::cluster_louvain} \code{resolution} argument). Higher
#'   values favour smaller communities. Default \code{1}.
#' @param steps Integer. Number of random walk steps for the Walktrap method.
#'   Default \code{4L}.
#'
#' @return The input \code{catgraph} object with two additions:
#' \describe{
#'   \item{\code{graph}}{Vertex attribute \code{cluster} (integer membership
#'     vector) and \code{cluster_method} (character) are added to the graph.}
#'   \item{\code{clustering}}{A named list with the raw \pkg{igraph}
#'     communities object (\code{communities}), the membership vector
#'     (\code{membership}), the number of communities (\code{n_clusters}),
#'     and the modularity score (\code{modularity}).}
#' }
#'
#' @details
#' Cluster detection is performed on the graph at the time of the call. If
#' you want to cluster a pruned graph, call \code{\link{prune_edges}} first.
#' This function clusters vertices of the \emph{variable-level} association
#' graph. Since vertices represent variables rather than factor levels or
#' respondents, the detected communities should be read as blocks of
#' variables that co-vary pairwise. For community detection at the modality
#' level (factor-level co-association), see \code{\link{cluster_modalities}}.
#'
#' The Louvain algorithm is non-deterministic; set a seed with
#' \code{set.seed()} before calling if reproducibility is needed.
#'
#' @references
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
#'   Fast unfolding of communities in large networks.
#'   \emph{Journal of Statistical Mechanics: Theory and Experiment}, 2008(10),
#'   P10008. \doi{10.1088/1742-5468/2008/10/P10008}
#'
#' Pons, P., & Latapy, M. (2005). Computing communities in large networks
#'   using random walks. In \emph{Computer and Information Sciences – ISCIS
#'   2005} (pp. 284--293). Springer. \doi{10.1007/11569596_31}
#'
#' @examples
#' df <- as.data.frame(Titanic)
#' df_exp <- df[rep(seq_len(nrow(df)), df$Freq), -5]
#' cg <- catgraph(df_exp)
#' set.seed(42)
#' cg <- detect_clusters(cg)
#' cg$clustering$n_clusters
#' cg$clustering$modularity
#'
#' @seealso \code{\link{catgraph}}, \code{\link{prune_edges}}
#' @importFrom igraph cluster_louvain cluster_walktrap cluster_fast_greedy
#'   cluster_label_prop cluster_edge_betweenness membership modularity V
#' @export
detect_clusters <- function(x,
                             method     = c("louvain", "walktrap",
                                            "fast_greedy", "label_prop",
                                            "edge_betweenness"),
                             resolution = 1,
                             steps      = 4L) {

  if (!inherits(x, "catgraph")) {
    stop("`x` must be a catgraph object.")
  }

  method <- match.arg(method)
  g      <- x$graph

  comm <- switch(
    method,
    louvain = igraph::cluster_louvain(
      g,
      weights    = igraph::E(g)$weight,
      resolution = resolution
    ),
    walktrap = igraph::cluster_walktrap(
      g,
      weights = igraph::E(g)$weight,
      steps   = steps
    ),
    fast_greedy = igraph::cluster_fast_greedy(
      g,
      weights = igraph::E(g)$weight
    ),
    label_prop = igraph::cluster_label_prop(
      g,
      weights = igraph::E(g)$weight
    ),
    edge_betweenness = igraph::cluster_edge_betweenness(
      g,
      weights  = igraph::E(g)$weight,
      directed = FALSE
    )
  )

  mem <- igraph::membership(comm)
  igraph::V(g)$cluster        <- as.integer(mem)
  igraph::V(g)$cluster_method <- method

  x$graph      <- g
  x$clustering <- list(
    communities = comm,
    membership  = mem,
    n_clusters  = max(mem),
    modularity  = igraph::modularity(comm)
  )

  x
}
