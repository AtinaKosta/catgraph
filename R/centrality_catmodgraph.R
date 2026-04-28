# =============================================================================
# centrality_catmodgraph.R
# S3 method: node_centrality() for catmodgraph objects
#
# Extends node_centrality() to return traditional centrality measures
# augmented with gravity indices (MGI+, MGI-, dMGI, OS, role) in a
# single unified table.
# =============================================================================

# -----------------------------------------------------------------------------
#' node_centrality method for catmodgraph objects
#'
#' Extends \code{\link{node_centrality}} to accept a \code{catmodgraph}
#' object (modality-level network).  Returns the standard centrality
#' measures - strength, betweenness, closeness, eigenvector, PageRank -
#' augmented with the gravity indices MGI+, MGI-, dMGI, OS, and role from
#' \code{\link{modality_gravity}}.
#'
#' @param x A \code{catmodgraph} object as returned by
#'   \code{\link{build_modality_graph}} or \code{\link{prune_modality_edges}}.
#' @param normalize Logical. If \code{TRUE} (default), the five traditional
#'   centrality columns are normalised to \eqn{[0, 1]}.  Gravity indices are
#'   never normalised here (use \code{mgi_plus_norm} for a normalised MGI+).
#' @param ... Additional arguments passed to \code{\link{modality_gravity}}.
#'
#' @return A data frame with one row per modality node and columns:
#'   \code{node}, \code{variable}, \code{modality},
#'   \code{prevalence}, \code{degree},
#'   \code{strength}, \code{w_betweenness}, \code{w_closeness},
#'   \code{w_eigenvector}, \code{w_pagerank}
#'   (traditional, optionally normalised), plus
#'   \code{mgi_plus}, \code{mgi_minus}, \code{delta_mgi},
#'   \code{mgi_plus_norm}, \code{os}, \code{role}
#'   (gravity indices, always on their natural scale).
#'   Rows are ordered by \code{delta_mgi} descending.
#'
#' @seealso \code{\link{modality_gravity}}, \code{\link{node_centrality}},
#'   \code{\link{plot_gravity}}
#'
#' @examples
#' data(survey_health)
#' mg  <- build_modality_graph(survey_health)
#' mg  <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
#' nc  <- node_centrality(mg)
#' print(nc[, c("node", "strength", "w_eigenvector", "delta_mgi", "os", "role")])
#'
#' @export
node_centrality.catmodgraph <- function(x, normalize = TRUE, ...) {
  
  if (!inherits(x, "catmodgraph")) {
    stop("`x` must be a catmodgraph object.", call. = FALSE)
  }
  
  g     <- x$graph
  nodes <- igraph::V(g)$name
  p     <- length(nodes)
  
  # ---- Traditional centrality on the modality graph -----------------------
  if (igraph::ecount(g) == 0L) {
    trad <- data.frame(
      node          = nodes,
      strength      = rep(0, p),
      w_betweenness = rep(0, p),
      w_closeness   = rep(0, p),
      w_eigenvector = rep(0, p),
      w_pagerank    = rep(1 / p, p),
      stringsAsFactors = FALSE
    )
  } else {
    w     <- abs(igraph::edge_attr(g, "weight"))
    inv_w <- 1 / pmax(w, 1e-9)
    
    trad <- data.frame(
      node          = nodes,
      strength      = as.numeric(igraph::strength(g, weights = w)),
      w_betweenness = as.numeric(igraph::betweenness(g,
                                                     weights = inv_w, normalized = FALSE)),
      w_closeness   = as.numeric(igraph::closeness(g,
                                                   weights = inv_w, normalized = FALSE)),
      w_eigenvector = as.numeric(
        igraph::eigen_centrality(g, weights = w)$vector),
      w_pagerank    = as.numeric(igraph::page_rank(g, weights = w)$vector),
      stringsAsFactors = FALSE
    )
    
    # closeness returns NA for nodes in disconnected components - replace with 0
    trad$w_closeness[is.na(trad$w_closeness)] <- 0
  }
  
  # Normalise traditional measures to [0, 1]
  if (normalize) {
    num_cols <- c("strength", "w_betweenness", "w_closeness",
                  "w_eigenvector", "w_pagerank")
    for (col in num_cols) {
      v  <- trad[[col]]
      v[is.na(v)] <- 0          # defensive NA replacement before normalising
      mx <- max(v, na.rm = TRUE)
      trad[[col]] <- if (is.finite(mx) && mx > 0) v / mx else v
    }
  }
  
  # ---- Gravity indices -----------------------------------------------------
  grav <- modality_gravity(x, ...)
  
  # Align gravity to vertex order
  grav_aligned <- grav[match(nodes, grav$node), ]
  
  # ---- Merge ---------------------------------------------------------------
  out <- data.frame(
    node          = nodes,
    variable      = grav_aligned$variable,
    modality      = grav_aligned$modality,
    prevalence    = grav_aligned$prevalence,
    degree        = grav_aligned$degree,
    strength      = trad$strength,
    w_betweenness = trad$w_betweenness,
    w_closeness   = trad$w_closeness,
    w_eigenvector = trad$w_eigenvector,
    w_pagerank    = trad$w_pagerank,
    mgi_plus      = grav_aligned$mgi_plus,
    mgi_minus     = grav_aligned$mgi_minus,
    delta_mgi     = grav_aligned$delta_mgi,
    mgi_plus_norm = grav_aligned$mgi_plus_norm,
    os            = grav_aligned$os,
    role          = grav_aligned$role,
    stringsAsFactors = FALSE
  )
  
  out <- out[order(-out$delta_mgi), ]
  rownames(out) <- NULL
  out
  
}