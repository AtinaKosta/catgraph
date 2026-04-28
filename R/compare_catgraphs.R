#' Compare multiple catgraph networks on one panel
#'
#' Visualises two or more \code{catgraph} objects side-by-side with a
#' shared node layout, for exploring how variable-level association
#' structure differs across populations, sites, time points, or other
#' grouping variables. All graphs must be built from the same variable
#' set (same column names).
#'
#' @param x A named list of \code{catgraph} objects. The list names
#'   become panel titles. A minimum of two graphs is required.
#' @param pruning Character. How to filter edges before plotting. One of:
#'   \describe{
#'     \item{\code{"pooled"} (default)}{Build a reference catgraph on
#'       all rows combined, apply BH-adjusted pruning at
#'       \code{max_p} and \code{min_weight}, and retain only those
#'       edges in every panel. This prevents sample-size artifacts:
#'       panels never differ in edge \emph{presence} solely because
#'       one group had fewer observations.}
#'     \item{\code{"individual"}}{Apply \code{\link{prune_edges}} to
#'       each graph independently at the specified thresholds.
#'       Differences in edge presence across panels are then shown
#'       directly, but may reflect power differences as well as
#'       substantive differences.}
#'     \item{\code{"overlay"}}{Show all edges in every panel,
#'       with thin grey edges for non-significant pairs and thick
#'       coloured edges for those surviving \code{"individual"}
#'       pruning. Most information-dense; busiest visually.}
#'     \item{\code{"none"}}{No filtering; every pair is shown.}
#'   }
#' @param min_weight Numeric. Effect-size threshold used by the
#'   \code{"pooled"}, \code{"individual"}, and \code{"overlay"} modes.
#'   Default \code{0.05}.
#' @param max_p Numeric. Adjusted p-value threshold used by the same
#'   modes. Default \code{0.05}.
#' @param p_adjust Character. Multiple-testing correction method, passed
#'   to \code{\link{prune_edges}}. Default \code{"BH"}.
#' @param layout_fn Function. An \pkg{igraph} layout function applied
#'   to the union graph to produce the shared node coordinates. Default
#'   \code{igraph::layout_with_fr}.
#' @param edge_width_range Numeric vector of length 2. Min and max edge
#'   widths when rescaling edge weights for display. Default
#'   \code{c(0.5, 4)}.
#' @param vertex_size Numeric. Vertex size. Default \code{28}.
#' @param ... Further arguments passed to \code{\link[igraph]{plot.igraph}}.
#'
#' @return Invisibly returns the reference union graph used to compute
#'   the shared layout (useful for further inspection). Called for its
#'   side effect: drawing the multi-panel comparison.
#'
#' @details
#' The default \code{"pooled"} mode is the most statistically conservative
#' choice for cross-group comparison. It ensures every panel shows the
#' same edge set, so weight differences across panels are interpretable
#' as substantive differences rather than power differences. Use
#' \code{"individual"} only when power differences themselves are the
#' scientific object of interest (e.g., to document that group A has
#' enough data to detect an association that group B does not).
#'
#' \strong{Note on \code{"pooled"} and mixed estimators.} When pooled
#' pruning is used, the reference catgraph is rebuilt with the
#' \code{corrected} flag copied from the \emph{first} element of
#' \code{x}. If your list mixes corrected and uncorrected graphs, pool
#' manually and call \code{compare_catgraphs(..., pruning = "none")}
#' on the pre-pruned objects instead.
#'
#' @section Formal testing:
#' This function visualises differences; it does not test them. For
#' permutation-based inferential comparison of modality-level
#' networks, see \code{\link{test_modality_graph_equality}} and
#' \code{\link{test_modality_edge_differences}}.
#'
#' @examples
#' # Split HairEyeColor into two populations and compare
#' df <- expand_table(HairEyeColor)
#' df_f <- df[df$Sex == "Female", c("Hair", "Eye")]
#' df_m <- df[df$Sex == "Male",   c("Hair", "Eye")]
#'
#' cg_f <- catgraph(df_f, corrected = TRUE)
#' cg_m <- catgraph(df_m, corrected = TRUE)
#'
#' compare_catgraphs(list(Female = cg_f, Male = cg_m))
#'
#' @seealso \code{\link{catgraph}}, \code{\link{prune_edges}},
#'   \code{\link{compare_modality_graphs}}
#' @importFrom graphics par
#' @importFrom igraph layout_with_fr as_edgelist add_edges simplify
#'   E V vcount ecount delete_edges induced_subgraph edge_attr_names
#' @export
compare_catgraphs <- function(x,
                              pruning          = c("pooled", "individual",
                                                   "overlay", "none"),
                              min_weight       = 0.05,
                              max_p            = 0.05,
                              p_adjust         = "BH",
                              layout_fn        = igraph::layout_with_fr,
                              edge_width_range = c(0.5, 4),
                              vertex_size      = 28,
                              ...) {

  .validate_compare_input(x, class_name = "catgraph")
  pruning <- match.arg(pruning)

  # Check all graphs share the same variable set
  var_sets <- lapply(x, function(cg) sort(igraph::V(cg$graph)$name))
  if (!all(vapply(var_sets, identical, logical(1L), var_sets[[1L]]))) {
    stop("All catgraph objects must share the same variable set (node names).",
         call. = FALSE)
  }

  # ------- Pooled pruning: compute reference edge set from pooled data -----
  retained_edges <- NULL
  if (pruning == "pooled") {
    pooled_data <- do.call(rbind, lapply(x, function(cg) cg$data))
    pooled_cg <- suppressWarnings(
      catgraph(pooled_data, corrected = x[[1L]]$corrected)
    )
    pooled_pruned <- prune_edges(
      pooled_cg, min_weight = min_weight, max_p = max_p, p_adjust = p_adjust
    )
    el <- igraph::as_edgelist(pooled_pruned$graph)
    if (nrow(el) > 0L) {
      retained_edges <- apply(el, 1, function(e) paste(sort(e), collapse = "\r"))
    } else {
      retained_edges <- character(0)
    }
  }

  # ------- Build per-panel graphs according to pruning choice --------------
  panel_graphs <- lapply(x, function(cg) {
    g <- cg$graph
    if (pruning == "none") return(g)

    if (pruning == "individual") {
      cg2 <- prune_edges(cg, min_weight = min_weight, max_p = max_p,
                         p_adjust = p_adjust)
      return(cg2$graph)
    }

    if (pruning == "pooled") {
      el <- igraph::as_edgelist(g)
      if (nrow(el) == 0L) return(g)
      keys <- apply(el, 1, function(e) paste(sort(e), collapse = "\r"))
      drop <- !(keys %in% retained_edges)
      if (any(drop)) g <- igraph::delete_edges(g, which(drop))
      return(g)
    }

    # overlay: keep full graph, flag survivors via edge attribute
    cg2 <- prune_edges(cg, min_weight = min_weight, max_p = max_p,
                       p_adjust = p_adjust)
    el_p <- igraph::as_edgelist(cg2$graph)
    keep_keys <- if (nrow(el_p) > 0L) {
      apply(el_p, 1, function(e) paste(sort(e), collapse = "\r"))
    } else character(0)

    el <- igraph::as_edgelist(g)
    keys <- if (nrow(el) > 0L) {
      apply(el, 1, function(e) paste(sort(e), collapse = "\r"))
    } else character(0)
    igraph::E(g)$significant <- keys %in% keep_keys
    g
  })

  # ------- Shared layout on the union graph --------------------------------
  union_g <- .union_graph(panel_graphs)
  lay     <- .shared_layout(union_g, layout_fn)   # no rescale: matches prior behaviour

  scale_w <- function(w, r = edge_width_range) {
    if (length(w) == 0L) return(numeric(0))
    if (diff(range(w, na.rm = TRUE)) < .Machine$double.eps) {
      return(rep(mean(r), length(w)))
    }
    r[1L] + (w - min(w, na.rm = TRUE)) /
      (max(w, na.rm = TRUE) - min(w, na.rm = TRUE)) * diff(r)
  }

  # ------- Plot multi-panel layout -----------------------------------------
  n <- length(panel_graphs)
  op <- graphics::par(mfrow = c(1, n), mar = c(2, 2, 3, 2))
  on.exit(graphics::par(op), add = TRUE)

  for (i in seq_len(n)) {
    g <- panel_graphs[[i]]
    title_i <- sprintf("%s (%d edges)", names(panel_graphs)[i],
                       igraph::ecount(g))

    if (pruning == "overlay" &&
        "significant" %in% igraph::edge_attr_names(g)) {
      sig <- igraph::E(g)$significant
      edge_cols <- ifelse(sig, "#534AB7", "#D3D1C7")
      edge_w    <- ifelse(sig, scale_w(igraph::E(g)$weight), 0.3)
    } else {
      edge_cols <- "#534AB7"
      edge_w    <- scale_w(igraph::E(g)$weight)
    }

    plot(
      g,
      layout             = .align_layout(g, lay),
      vertex.size        = vertex_size,
      vertex.color       = "#AFA9EC",
      vertex.label.cex   = 0.75,
      vertex.label.color = "#2C2C2A",
      edge.width         = edge_w,
      edge.color         = edge_cols,
      main               = title_i,
      ...
    )
  }

  invisible(union_g)
}
