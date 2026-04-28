#' Compute weighted centrality indices for all variables in a catgraph
#'
#' Returns a data frame of five weighted centrality measures for every node
#' (variable) in a \code{catgraph} object. All measures use the phi /
#' Cramer's V edge weights stored in the graph, so they reflect association
#' strength rather than mere connectivity.
#'
#' @param x A \code{catgraph} object.
#'
#' @return A \code{data.frame} with one row per variable and the following
#'   columns:
#' \describe{
#'   \item{\code{variable}}{Variable name (character).}
#'   \item{\code{strength}}{Weighted degree / strength centrality: the sum
#'     of all edge weights incident to a node. The primary measure of how
#'     strongly a variable is globally associated with all others
#'     (Barrat et al., 2004).}
#'   \item{\code{w_betweenness}}{Weighted betweenness centrality: the
#'     number of weighted shortest paths passing through a node, where path
#'     lengths are \eqn{1/w_{ij}} so that stronger associations correspond
#'     to shorter distances (Brandes, 2001).}
#'   \item{\code{w_closeness}}{Weighted closeness centrality: the inverse
#'     of the mean weighted distance to all other reachable nodes.}
#'   \item{\code{w_eigenvector}}{Weighted eigenvector centrality: a node
#'     scores highly if it is connected to other high-scoring nodes
#'     (Bonacich, 1987).}
#'   \item{\code{w_pagerank}}{Weighted PageRank: a random-walk based
#'     centrality robust to isolated subgraphs (Brin & Page, 1998).
#'     Uses the \pkg{igraph} default damping factor of 0.85.}
#' }
#' The data frame is returned sorted by \code{strength} in descending order.
#'
#' @details
#' Edges with \code{NA} or non-positive weights are dropped from a local
#' copy of the graph before centrality is computed; the input
#' \code{catgraph} object is never modified. If the graph has no edges
#' after this cleanup, a zero-valued result is returned (with a uniform
#' \code{1/p} PageRank, the stationary distribution on an edgeless graph).
#'
#' \strong{Strength centrality (Barrat et al., 2004):}
#'
#' \deqn{s_i = \sum_{j \in \mathcal{N}(i)} w_{ij}}
#'
#' where \eqn{\mathcal{N}(i)} is the set of neighbours of node \eqn{i} and
#' \eqn{w_{ij}} is the edge weight (phi or Cramer's V). In a categorical
#' association graph this is the most interpretable measure: it is simply
#' the total association strength of a variable with all other variables.
#'
#' \strong{Weighted betweenness (Brandes, 2001):}
#'
#' Path lengths are defined as \eqn{1/w_{ij}} so that stronger edges are
#' traversed preferentially. A variable with high weighted betweenness sits
#' on many shortest paths and therefore acts as a statistical bridge between
#' otherwise weakly associated variable groups.
#'
#' \strong{Weighted closeness:}
#'
#' \deqn{C_i = \frac{n-1}{\sum_{j \neq i} d(i,j)}}
#'
#' where \eqn{d(i,j)} is the weighted shortest-path distance. A high value
#' indicates a variable that reaches all other variables through short
#' (strong-association) paths.
#'
#' \strong{Eigenvector centrality (Bonacich, 1987):}
#'
#' The leading eigenvector of the weighted adjacency matrix. A variable
#' scores highly when it is associated with other variables that are
#' themselves highly associated variables — a second-order hub measure.
#'
#' \strong{PageRank (Brin & Page, 1998):}
#'
#' A random-walk measure that down-weights the contribution of nodes with
#' many weak connections. Robust when the graph is sparse or has near-zero
#' weight edges.
#'
#' @references
#' Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A. (2004).
#'   The architecture of complex weighted networks.
#'   \emph{Proceedings of the National Academy of Sciences}, 101(11),
#'   3747--3752. \doi{10.1073/pnas.0400087101}
#'
#' Bonacich, P. (1987). Power and centrality: A family of measures.
#'   \emph{American Journal of Sociology}, 92(5), 1170--1182.
#'   \doi{10.1086/228631}
#'
#' Brandes, U. (2001). A faster algorithm for betweenness centrality.
#'   \emph{Journal of Mathematical Sociology}, 25(2), 163--177.
#'   \doi{10.1080/0022250X.2001.9990249}
#'
#' Brin, S., & Page, L. (1998). The anatomy of a large-scale hypertextual
#'   web search engine. \emph{Computer Networks}, 30(1-7), 107--117.
#'   \doi{10.1016/S0169-7552(98)00110-X}
#'
#' @examples
#' df <- expand_table(Titanic)
#' cg <- catgraph(df)
#' nc <- node_centrality(cg)
#' nc
#'
#' # Identify the most central variable
#' nc$variable[1]
#'
#' # Raw (unnormalised) values
#' node_centrality(cg, normalize = FALSE)
#'
#' @seealso \code{\link{plot_centrality}}, \code{\link{catgraph}},
#'   \code{\link{clustering_coef}}
#' @importFrom igraph strength betweenness closeness eigen_centrality page_rank
#' @importFrom graphics par barplot
# S3 generic: dispatches to node_centrality.catgraph (catgraph objects)
# or node_centrality.catmodgraph (modality networks, in modality_gravity.R)
#' @rdname node_centrality
#' @param x A \code{catgraph} or \code{catmodgraph} object.
#' @param ... Additional arguments passed to the method.
#' @export
node_centrality <- function(x, ...) UseMethod("node_centrality")

#' @export
node_centrality.catgraph <- function(x, normalize = TRUE, ...) {
  
  if (!inherits(x, "catgraph")) stop("`x` must be a catgraph object.")
  
  # Work on a local copy so we never mutate the caller's graph.
  g <- x$graph
  
  vars <- igraph::V(g)$name
  p    <- length(vars)
  
  # Empty-graph short circuit: return zero-valued centrality for every vertex
  # rather than propagating NaN/NA from igraph on a no-edge graph.
  if (igraph::ecount(g) == 0L) {
    out <- data.frame(
      variable      = vars,
      strength      = rep(0, p),
      w_betweenness = rep(0, p),
      w_closeness   = rep(0, p),
      w_eigenvector = rep(0, p),
      w_pagerank    = rep(1 / p, p),   # uniform stationary distribution
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    return(out)
  }
  
  # Drop edges with NA or non-positive weights from the LOCAL copy only.
  # After v0.4.0, build_graph() does not create zero-weight edges, so this
  # branch is defensive: it handles weights that became NA downstream (e.g.
  # after a user-supplied transformation) without polluting shortest-path
  # distances with a synthetic epsilon.
  w <- igraph::E(g)$weight
  bad <- is.na(w) | w <= 0
  if (any(bad)) {
    g <- igraph::delete_edges(g, which(bad))
  }
  
  # Recompute after potential edge deletion.
  if (igraph::ecount(g) == 0L) {
    out <- data.frame(
      variable      = vars,
      strength      = rep(0, p),
      w_betweenness = rep(0, p),
      w_closeness   = rep(0, p),
      w_eigenvector = rep(0, p),
      w_pagerank    = rep(1 / p, p),
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    return(out)
  }
  
  w <- igraph::E(g)$weight
  
  # Strength: sum of incident edge weights (Barrat et al., 2004)
  s <- igraph::strength(g, weights = w)
  
  # Weighted betweenness / closeness: path length is 1/w so stronger edges
  # are shorter distances (Newman, 2001). No epsilon floor needed because
  # we deleted non-positive edges above.
  inv_w <- 1 / w
  
  btw <- igraph::betweenness(g, weights = inv_w, normalized = FALSE)
  clo <- igraph::closeness(g,  weights = inv_w, normalized = FALSE)
  
  eig <- igraph::eigen_centrality(g, weights = w)$vector
  pgr <- igraph::page_rank(g,        weights = w)$vector
  
  # Align to the original vertex order (delete_edges preserves vertex order,
  # but be defensive in case future igraph versions reorder).
  out <- data.frame(
    variable      = igraph::V(g)$name,
    strength      = as.numeric(s),
    w_betweenness = as.numeric(btw),
    w_closeness   = as.numeric(clo),
    w_eigenvector = as.numeric(eig),
    w_pagerank    = as.numeric(pgr),
    stringsAsFactors = FALSE
  )
  
  # Normalise to [0,1] per measure
  if (normalize) {
    num_cols <- c("strength", "w_betweenness", "w_closeness",
                  "w_eigenvector", "w_pagerank")
    for (col in num_cols) {
      vals <- out[[col]]
      if (length(vals) == 0L) next
      if (all(is.na(vals))) next
      mx <- suppressWarnings(max(vals, na.rm = TRUE))
      if (is.finite(mx) && mx > 0) {
        out[[col]] <- vals / mx
      }
    }
  }
  
  out <- out[order(out$strength, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}


#' Plot weighted centrality indices for a catgraph
#'
#' Produces a ranked horizontal bar chart of one or all centrality measures
#' from a \code{catgraph} object. When a single measure is selected the chart
#' shows bars ranked by that measure. When \code{measure = "all"}, a faceted
#' or overlaid comparison chart is produced.
#'
#' @param x A \code{catgraph} object.
#' @param measure Character. One of \code{"strength"}, \code{"w_betweenness"},
#'   \code{"w_closeness"}, \code{"w_eigenvector"}, \code{"w_pagerank"}, or
#'   \code{"all"}. Default \code{"strength"}.
#' @param normalize Logical. Passed to \code{\link{node_centrality}}.
#'   Default \code{TRUE}.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param color Character. Bar fill colour. Default \code{"#7F77DD"}
#'   (purple-400).
#' @param engine Character. \code{"ggplot2"} (default) or \code{"base"}.
#'
#' @return For \code{engine = "ggplot2"}: a \code{ggplot} object.
#'   For \code{engine = "base"}: \code{NULL}, invisibly.
#'
#' @examples
#' df <- expand_table(Titanic)
#' cg <- catgraph(df)
#' plot_centrality(cg)
#' plot_centrality(cg, measure = "w_betweenness")
#' # plot_centrality(cg, measure = "all")  # requires ggplot2
#'
#' @seealso \code{\link{node_centrality}}
#' @importFrom graphics par barplot
#' @importFrom stats reshape
#' @export
plot_centrality <- function(x,
                            measure   = "strength",
                            normalize = TRUE,
                            title     = NULL,
                            color     = "#7F77DD",
                            engine    = c("ggplot2", "base")) {
  
  engine  <- match.arg(engine)
  nc      <- node_centrality(x, normalize = normalize)
  measures <- c("strength", "w_betweenness", "w_closeness",
                "w_eigenvector", "w_pagerank")
  
  if (measure != "all" && !measure %in% measures) {
    stop("`measure` must be one of: ", paste(c(measures, "all"), collapse = ", "))
  }
  
  if (engine == "base") {
    if (measure == "all") measure <- "strength"
    vals <- nc[[measure]]
    ord  <- order(vals)
    labs <- nc$variable[ord]
    graphics::par(mar = c(4, 8, 3, 2))
    graphics::barplot(
      vals[ord],
      names.arg = labs,
      horiz     = TRUE,
      las       = 1,
      col       = color,
      border    = NA,
      xlab      = if (normalize) paste(measure, "(normalised)") else measure,
      main      = if (!is.null(title)) title else paste("Node centrality:", measure)
    )
    return(invisible(NULL))
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }
  
  if (measure == "all") {
    # Long format for faceting
    long <- stats::reshape(
      nc,
      varying   = measures,
      v.names   = "value",
      timevar   = "measure",
      times     = measures,
      direction = "long"
    )
    long$measure  <- factor(long$measure,  levels = measures)
    long$variable <- factor(long$variable,
                            levels = nc$variable[order(nc$strength)])
    
    p <- ggplot2::ggplot(long,
                         ggplot2::aes(x = .data$value, y = .data$variable)) +
      ggplot2::geom_col(fill = color, width = 0.7) +
      ggplot2::facet_wrap(~ .data$measure, scales = "free_x", nrow = 2) +
      ggplot2::labs(
        x     = if (normalize) "Normalised value" else "Value",
        y     = NULL,
        title = title
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
    
  } else {
    nc$variable <- factor(nc$variable,
                          levels = nc$variable[order(nc[[measure]])])
    
    label_map <- c(
      strength      = "Strength (sum of weights)",
      w_betweenness = "Weighted betweenness",
      w_closeness   = "Weighted closeness",
      w_eigenvector = "Weighted eigenvector",
      w_pagerank    = "Weighted PageRank"
    )
    
    p <- ggplot2::ggplot(nc,
                         ggplot2::aes(x = .data[[measure]], y = .data$variable)) +
      ggplot2::geom_col(fill = color, width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = formatC(.data[[measure]], digits = 3, format = "f")),
        hjust = -0.1, size = 3, color = "#444441"
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.15))
      ) +
      ggplot2::labs(
        x     = if (normalize) paste(label_map[measure], "(normalised)")
        else label_map[measure],
        y     = NULL,
        title = title
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  }
  
  p
}