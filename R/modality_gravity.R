# =============================================================================
# modality_gravity.R
# Modality Gravity Index (MGI) and Orbital Score (OS) for catmodgraph objects
#
# Exported functions:
#   modality_gravity()   -- compute gravity indices per modality node
#   plot_gravity()       -- 6-panel comparison: traditional vs gravity indices
#   compare_gravity()    -- compare gravity profiles across two subgroups
#
# Theoretical background:
#   Standard centrality indices treat nodes as exchangeable except for their
#   connectivity structure. For modality networks, nodes are categorical
#   response options with empirical prevalences; a rare modality connected to
#   a common one occupies a fundamentally different structural position than
#   two equally prevalent modalities with the same edge weight.
#
#   MGI corrects for this by weighting each edge contribution by the
#   prevalence ratio p_i / p_j:
#
#     MGI+(i)  = sum_{j in N(i), p_j < p_i}  w_ij * (p_i / p_j)
#     MGI-(i)  = sum_{j in N(i), p_j > p_i}  w_ij * (p_i / p_j)
#     dMGI(i)  = MGI+(i) - MGI-(i)    [signed; + = attractor, - = satellite]
#     OS(i)    = sum_{j in N(i)}       w_ij * (p_j / p_i)  [orbital pull]
#
#   When all nodes are equally prevalent, p_i/p_j = 1 for all edges and
#   MGI+(i) + MGI-(i) reduces to weighted degree (strength). MGI is
#   therefore a generalisation of strength to prevalence-heterogeneous graphs.
# =============================================================================


# -----------------------------------------------------------------------------
#' Modality Gravity Index for catmodgraph objects
#'
#' Computes four prevalence-weighted structural indices for every node in a
#' modality-level network: gravitational mass (MGI+), orbital drag (MGI-),
#' net gravity (dMGI), and orbital score (OS).  These indices extend standard
#' centrality measures by incorporating the empirical prevalence of each
#' modality, allowing direction-sensitive identification of attractor and
#' satellite modalities.
#'
#' @param x A \code{catmodgraph} object as returned by
#'   \code{\link{build_modality_graph}} or \code{\link{prune_modality_edges}}.
#'   The object must contain a \code{$data} slot with the original data frame
#'   from which prevalences are computed.
#' @param weight_attr Character. Name of the edge attribute to use as
#'   association weight. Defaults to \code{"weight"} (absolute phi /
#'   Cramer's V stored by \code{build_modality_graph}).
#' @param min_prevalence Numeric in (0, 1). Modalities with prevalence below
#'   this value are retained in the output but their prevalence-ratio
#'   contributions are flagged. Default \code{0.01} (1\%).
#'
#' @return A data frame with one row per modality node and columns:
#'   \describe{
#'     \item{node}{Full node label in \code{variable=modality} format.}
#'     \item{variable}{Variable name (left of \code{=}).}
#'     \item{modality}{Modality label (right of \code{=}).}
#'     \item{prevalence}{Empirical relative frequency computed from
#'       \code{x$data}, ignoring \code{NA}s.}
#'     \item{degree}{Number of edges incident to the node.}
#'     \item{strength}{Weighted degree: sum of incident edge weights.}
#'     \item{mgi_plus}{MGI+: sum of \eqn{w_{ij} \cdot p_i/p_j} over
#'       neighbours \eqn{j} with \eqn{p_j < p_i}.}
#'     \item{mgi_minus}{MGI-: sum of \eqn{w_{ij} \cdot p_i/p_j} over
#'       neighbours \eqn{j} with \eqn{p_j > p_i}.}
#'     \item{delta_mgi}{Net gravity: \code{mgi_plus - mgi_minus}.
#'       Positive = net attractor; negative = net satellite.}
#'     \item{mgi_plus_norm}{MGI+ normalised to \eqn{[0, 1]} within the
#'       graph, for cross-graph comparisons.}
#'     \item{os}{Orbital Score: sum of \eqn{w_{ij} \cdot p_j/p_i} over
#'       all neighbours. High OS indicates strong gravitational capture.}
#'     \item{role}{Character: \code{"attractor"} if \code{delta_mgi > 0},
#'       \code{"satellite"} if \code{delta_mgi < 0}, \code{"neutral"} if
#'       \code{delta_mgi == 0} (isolated nodes or perfectly balanced).}
#'   }
#'   Rows are ordered by \code{delta_mgi} descending (strongest attractors
#'   first).
#'
#' @details
#' \strong{Prevalence computation.}  Prevalences are computed from
#' \code{x$data} as \eqn{n_{ij} / n_j} where \eqn{n_{ij}} is the number of
#' non-missing observations of variable \eqn{j} taking value \eqn{i}, and
#' \eqn{n_j} is the total non-missing count for that variable.  Each variable
#' is normalised independently, so prevalences sum to 1 within each variable.
#'
#' \strong{Relationship to standard centrality.}
#' When all nodes in the graph have equal prevalence, \eqn{p_i/p_j = 1} for
#' every edge, and \code{mgi_plus + mgi_minus} equals weighted degree
#' (strength).  \code{delta_mgi} is then determined solely by whether a
#' node's neighbours are uniformly more or less prevalent, which is
#' uninformative.  The metric is most useful when prevalences vary
#' substantially across modalities.
#'
#' \strong{Edge weights.}  The function takes the absolute value of edge
#' weights before computing ratios so that signed phi values (stored for
#' the modality graph) do not cancel positive and negative contributions.
#' The sign of the association is a separate dimension captured in
#' \code{std_resid}; \code{delta_mgi} captures prevalence-weighted
#' directional dominance.
#'
#' @seealso \code{\link{plot_gravity}}, \code{\link{compare_gravity}},
#'   \code{\link{node_centrality}}, \code{\link{build_modality_graph}}
#'
#' @examples
#' data(survey_health)
#' mg <- build_modality_graph(survey_health)
#' mg <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
#' grav <- modality_gravity(mg)
#' print(grav)
#'
#' # Top attractors
#' head(grav[grav$role == "attractor", ], 5)
#'
#' # Top satellites (highest orbital score)
#' head(grav[order(-grav$os), ], 5)
#'
#' @export
modality_gravity <- function(x,
                             weight_attr   = "weight",
                             min_prevalence = 0.01) {
  
  # ---- Input validation ----------------------------------------------------
  if (!inherits(x, "catmodgraph")) {
    stop("`x` must be a catmodgraph object (from build_modality_graph()).",
         call. = FALSE)
  }
  
  g    <- x$graph
  data <- x$data
  
  if (is.null(data) || !is.data.frame(data)) {
    stop(
      "The catmodgraph object has no '$data' slot. ",
      "Rebuild with build_modality_graph() to attach the source data.",
      call. = FALSE
    )
  }
  
  if (!weight_attr %in% igraph::edge_attr_names(g)) {
    stop(sprintf(
      "Edge attribute '%s' not found. Available: %s",
      weight_attr,
      paste(igraph::edge_attr_names(g), collapse = ", ")
    ), call. = FALSE)
  }
  
  if (igraph::ecount(g) == 0L) {
    warning("Graph has no edges; all gravity indices will be zero.",
            call. = FALSE)
  }
  
  # ---- Prevalence computation ----------------------------------------------
  nodes <- igraph::V(g)$name
  
  prev <- vapply(nodes, function(nd) {
    var <- sub("=.*", "", nd)
    mod <- sub(".*=", "", nd)
    if (!var %in% names(data)) return(NA_real_)
    col <- as.character(data[[var]])
    n_total <- sum(!is.na(col))
    if (n_total == 0L) return(NA_real_)
    sum(col == mod, na.rm = TRUE) / n_total
  }, numeric(1L))
  
  low_prev <- !is.na(prev) & prev > 0 & prev < min_prevalence
  if (any(low_prev)) {
    warning(sprintf(
      "%d modality/modalities have prevalence below %.2f: %s\n",
      sum(low_prev), min_prevalence,
      paste(nodes[low_prev], collapse = ", ")
    ), call. = FALSE)
  }
  
  # ---- Initialise accumulators ---------------------------------------------
  mgi_plus  <- setNames(numeric(length(nodes)), nodes)
  mgi_minus <- setNames(numeric(length(nodes)), nodes)
  os_vec    <- setNames(numeric(length(nodes)), nodes)
  
  # ---- Edge loop -----------------------------------------------------------
  el <- igraph::as_data_frame(g, what = "edges")
  
  for (k in seq_len(nrow(el))) {
    i   <- el$from[k]
    j   <- el$to[k]
    w   <- abs(el[[weight_attr]][k])
    p_i <- prev[i]
    p_j <- prev[j]
    
    if (anyNA(c(p_i, p_j)) || p_i == 0 || p_j == 0) next
    
    # --- node i ---
    r_ij <- p_i / p_j
    if (p_i > p_j) {
      mgi_plus[i]  <- mgi_plus[i]  + w * r_ij
    } else if (p_i < p_j) {
      mgi_minus[i] <- mgi_minus[i] + w * r_ij
    }
    
    # --- node j ---
    r_ji <- p_j / p_i
    if (p_j > p_i) {
      mgi_plus[j]  <- mgi_plus[j]  + w * r_ji
    } else if (p_j < p_i) {
      mgi_minus[j] <- mgi_minus[j] + w * r_ji
    }
    
    # --- OS (both nodes) ---
    os_vec[i] <- os_vec[i] + w * (p_j / p_i)
    os_vec[j] <- os_vec[j] + w * (p_i / p_j)
  }
  
  # ---- Derived quantities --------------------------------------------------
  delta_mgi     <- mgi_plus - mgi_minus
  max_plus      <- max(mgi_plus, na.rm = TRUE)
  mgi_plus_norm <- if (max_plus > 0) mgi_plus / max_plus else mgi_plus
  
  role <- ifelse(delta_mgi > 0,  "attractor",
                 ifelse(delta_mgi < 0,  "satellite", "neutral"))
  
  # ---- Standard centrality for reference -----------------------------------
  deg      <- igraph::degree(g)
  strength <- igraph::strength(g,
                               weights = abs(igraph::edge_attr(g, weight_attr)))
  
  # ---- Assemble output -----------------------------------------------------
  out <- data.frame(
    node          = nodes,
    variable      = sub("=.*", "", nodes),
    modality      = sub(".*=", "", nodes),
    prevalence    = round(prev,          4L),
    degree        = as.integer(deg[nodes]),
    strength      = round(strength[nodes], 4L),
    mgi_plus      = round(mgi_plus,       4L),
    mgi_minus     = round(mgi_minus,      4L),
    delta_mgi     = round(delta_mgi,      4L),
    mgi_plus_norm = round(mgi_plus_norm,  4L),
    os            = round(os_vec,         4L),
    role          = role,
    row.names     = NULL,
    stringsAsFactors = FALSE
  )
  
  out <- out[order(-out$delta_mgi), ]
  rownames(out) <- NULL
  class(out) <- c("modality_gravity", "data.frame")
  out
}