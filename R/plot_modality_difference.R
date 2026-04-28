#' Plot modality-network differences on a single graph
#'
#' Renders the edge-wise differences between two \code{catmodgraph}
#' objects (from \code{\link{test_modality_edge_differences}}) on one
#' shared-layout graph, so the reader can see \emph{where} in the joint
#' categorical structure two samples disagree, not just \emph{that} they
#' disagree.
#'
#' Edge colour encodes the sign of the difference \code{weight_x - weight_y}:
#' edges stronger in \code{x} are drawn in one colour, edges stronger in
#' \code{y} in another. Edge width scales with \eqn{|weight_x - weight_y|}.
#' Edge opacity scales with \code{-log10(p_adjusted)} so that edges with
#' smaller adjusted p-values dominate the visual field; non-significant
#' edges remain visible but faded.
#'
#' This is a complement to the \code{plot()} method for
#' \code{catmodedgetest} objects, which
#' renders a bar chart of the top-n edges by adjusted p-value. Use the bar
#' chart for a ranked-list read and this function for a network-structural
#' read.
#'
#' @param x A \code{catmodedgetest} object, typically from
#'   \code{\link{test_modality_edge_differences}}.
#' @param reference Optional. One of the two source \code{catmodgraph}
#'   objects, or a list of both. Used only to derive the node set and
#'   (by default) the layout. If \code{NULL} (default), a graph is built
#'   from the union of edges in \code{x$edge_table} and laid out
#'   directly; node-variable colouring is still recovered from the
#'   \code{"Variable=level"} naming convention.
#' @param alpha_fdr Numeric in (0, 1]. Adjusted p-value cutoff for
#'   "significant" edges. Non-significant edges are drawn with a floor
#'   opacity set by \code{alpha_floor}. Default \code{0.05}.
#' @param alpha_floor Numeric in [0, 1]. Minimum edge opacity, applied
#'   to non-significant edges and as a lower bound on the
#'   \code{-log10(p_adjusted)} alpha mapping. Default \code{0.15}.
#' @param show_nonsig Logical. If \code{FALSE}, edges with
#'   \code{p_adjusted >= alpha_fdr} are omitted entirely. If \code{TRUE}
#'   (default), they are drawn at \code{alpha_floor}.
#' @param group_labels Character vector of length 2. Legend labels for
#'   the two groups; the first is the direction \code{obs_diff > 0}
#'   (stronger in \code{x}), the second \code{obs_diff < 0} (stronger
#'   in \code{y}). Default \code{c("stronger in x", "stronger in y")}.
#' @param color_pos,color_neg Character. Hex colours for the two
#'   directions. Defaults are a teal / coral pair chosen to match the
#'   signed-residual palette in \code{\link{plot.catmodgraph}}.
#' @param layout_fn Function. An \pkg{igraph} layout function applied
#'   to the reference graph. Default \code{igraph::layout_with_fr}.
#' @param vertex_size Numeric. Vertex size. Default \code{14}.
#' @param edge_scale Numeric. Multiplier for edge widths after scaling
#'   \eqn{|obs_diff|} into \code{[0, 1]}. Default \code{8}.
#' @param title Character. Plot title. Default \code{NULL} (auto).
#' @param ... Further arguments passed to \code{\link[igraph]{plot.igraph}}.
#'
#' @return Invisibly returns the igraph object used for plotting,
#'   with edge attributes \code{obs_diff} and \code{p_adjusted} set.
#'
#' @section Interpretation caveats:
#' \itemize{
#'   \item The sign of the difference is group-order-dependent: if
#'     \code{x} and \code{y} are swapped at the testing step, every
#'     colour flips. The \code{group_labels} argument lets you write the
#'     legend in plain language rather than relying on the reader to
#'     remember which group was which.
#'   \item Edges absent from both input graphs (weight 0 in each) will
#'     have \code{obs_diff == 0} and are omitted by
#'     \code{edges = "union"} in
#'     \code{\link{test_modality_edge_differences}}. If you called the
#'     test with \code{edges = "all"}, consider pre-filtering the
#'     edge table before plotting to avoid a dense low-magnitude
#'     background.
#'   \item The \code{-log10(p_adjusted)} alpha mapping is clamped at
#'     4 (i.e., \code{p_adjusted <= 1e-4} all plot at full opacity).
#'     This prevents a single ultra-significant edge from visually
#'     dominating the entire panel.
#' }
#'
#' @examples
#' \donttest{
#' data(survey_health)
#' df_f <- subset(survey_health, sex == "female")[, -1]
#' df_m <- subset(survey_health, sex == "male")[, -1]
#'
#' mg_f <- build_modality_graph(df_f)
#' mg_m <- build_modality_graph(df_m)
#'
#' edge_test <- test_modality_edge_differences(
#'   mg_f, mg_m, n_perm = 200, edges = "union",
#'   seed = 1, verbose = FALSE
#' )
#'
#' plot_modality_difference(
#'   edge_test,
#'   reference    = list(mg_f, mg_m),
#'   group_labels = c("stronger in women", "stronger in men")
#' )
#' }
#'
#' @seealso
#' \code{\link{test_modality_edge_differences}},
#' \code{plot.catmodedgetest()},
#' \code{\link{compare_modality_graphs}}
#' @importFrom igraph graph_from_data_frame layout_with_fr V E vcount
#'   ecount edge_attr_names
#' @importFrom grDevices adjustcolor hcl.colors
#' @importFrom graphics legend par
#' @importFrom scales rescale
#' @export
plot_modality_difference <- function(x,
                                     reference    = NULL,
                                     alpha_fdr    = 0.05,
                                     alpha_floor  = 0.15,
                                     show_nonsig  = TRUE,
                                     group_labels = c("stronger in x",
                                                      "stronger in y"),
                                     color_pos    = "#1D9E75",
                                     color_neg    = "#D85A30",
                                     layout_fn    = igraph::layout_with_fr,
                                     vertex_size  = 14,
                                     edge_scale   = 8,
                                     title        = NULL,
                                     ...) {
  
  if (!inherits(x, "catmodedgetest")) {
    stop("`x` must be a catmodedgetest object (from ",
         "test_modality_edge_differences()).", call. = FALSE)
  }
  if (!is.numeric(alpha_fdr) || length(alpha_fdr) != 1L ||
      is.na(alpha_fdr) || alpha_fdr <= 0 || alpha_fdr > 1) {
    stop("`alpha_fdr` must be a single number in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(alpha_floor) || length(alpha_floor) != 1L ||
      is.na(alpha_floor) || alpha_floor < 0 || alpha_floor > 1) {
    stop("`alpha_floor` must be a single number in [0, 1].", call. = FALSE)
  }
  if (!is.logical(show_nonsig) || length(show_nonsig) != 1L ||
      is.na(show_nonsig)) {
    stop("`show_nonsig` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(group_labels) || length(group_labels) != 2L) {
    stop("`group_labels` must be a length-2 character vector.",
         call. = FALSE)
  }
  
  et <- x$edge_table
  
  # -- Optionally filter non-significant edges ---------------------------------
  if (!show_nonsig) {
    et <- et[!is.na(et$p_adjusted) & et$p_adjusted < alpha_fdr, , drop = FALSE]
    if (nrow(et) == 0L) {
      stop("No edges survive alpha_fdr = ", alpha_fdr,
           "; nothing to plot. Try show_nonsig = TRUE.", call. = FALSE)
    }
  }
  
  # Drop exact-zero differences (no visible signal, and obs_diff = 0 has
  # ambiguous sign)
  et <- et[!is.na(et$obs_diff) & et$obs_diff != 0, , drop = FALSE]
  if (nrow(et) == 0L) {
    stop("All edge differences are zero or NA; nothing to plot.",
         call. = FALSE)
  }
  
  # -- Resolve the vertex set --------------------------------------------------
  # Priority: reference graph(s) > union of endpoints in et
  if (is.null(reference)) {
    node_names <- unique(c(et$from, et$to))
  } else {
    ref_list <- if (inherits(reference, "catmodgraph")) list(reference)
    else if (is.list(reference) &&
             all(vapply(reference, inherits, logical(1L),
                        "catmodgraph"))) reference
    else stop("`reference` must be a catmodgraph or a list of ",
              "catmodgraph objects.", call. = FALSE)
    node_names <- unique(unlist(lapply(ref_list, function(g)
      igraph::V(g$graph)$name)))
    # Include any edge-table endpoints missing from the reference (e.g.,
    # modalities that exist only in one group under restrict = "union")
    extra <- setdiff(unique(c(et$from, et$to)), node_names)
    if (length(extra) > 0L) node_names <- c(node_names, extra)
  }
  
  vertices <- data.frame(
    name     = node_names,
    variable = sub("=.*", "", node_names),
    stringsAsFactors = FALSE
  )
  
  # -- Build the difference graph ---------------------------------------------
  edges_df <- data.frame(
    from       = et$from,
    to         = et$to,
    obs_diff   = et$obs_diff,
    abs_diff   = abs(et$obs_diff),
    p_adjusted = et$p_adjusted,
    stringsAsFactors = FALSE
  )
  
  g <- igraph::graph_from_data_frame(
    d        = edges_df,
    directed = FALSE,
    vertices = vertices
  )
  
  # -- Edge aesthetics ---------------------------------------------------------
  obs_d  <- igraph::E(g)$obs_diff
  p_adj  <- igraph::E(g)$p_adjusted
  is_sig <- !is.na(p_adj) & p_adj < alpha_fdr
  
  # Width: scale |diff| into [0, 1] then multiply. Guard constant-width case.
  abs_d <- igraph::E(g)$abs_diff
  if (length(abs_d) > 0L && diff(range(abs_d)) < .Machine$double.eps) {
    width_scaled <- rep(0.5, length(abs_d))
  } else {
    width_scaled <- scales::rescale(abs_d, to = c(0, 1))
  }
  edge_w <- pmax(width_scaled, 0.05) * edge_scale
  
  # Alpha: -log10(p_adj) clamped to [alpha_floor, 1]; non-sig -> alpha_floor
  safe_p <- ifelse(is.na(p_adj) | p_adj <= 0, 1e-16, p_adj)
  alpha_sig <- pmin(pmax(-log10(safe_p) / 4, alpha_floor), 1)
  edge_alpha <- ifelse(is_sig, alpha_sig, alpha_floor)
  
  base_cols <- ifelse(obs_d > 0, color_pos, color_neg)
  edge_cols <- vapply(seq_along(edge_alpha),
                      function(i) grDevices::adjustcolor(base_cols[i],
                                                         alpha.f = edge_alpha[i]),
                      FUN.VALUE = character(1))
  
  # -- Node colouring by originating variable --------------------------------
  v_var    <- igraph::V(g)$variable
  var_lev  <- sort(unique(v_var))
  pal      <- grDevices::hcl.colors(length(var_lev), palette = "Dark 3")
  v_cols   <- pal[match(v_var, var_lev)]
  
  # -- Layout ----------------------------------------------------------------
  lay <- layout_fn(g)
  
  # -- Title -----------------------------------------------------------------
  if (is.null(title)) {
    n_sig <- sum(is_sig)
    title <- sprintf("Modality-network differences (%d of %d edges at FDR %.2f)",
                     n_sig, length(is_sig), alpha_fdr)
  }
  
  plot(
    g,
    layout             = lay,
    vertex.color       = v_cols,
    vertex.size        = vertex_size,
    vertex.label       = igraph::V(g)$name,
    vertex.label.cex   = 0.65,
    vertex.label.color = "black",
    edge.width         = edge_w,
    edge.color         = edge_cols,
    main               = title,
    ...
  )
  
  graphics::legend(
    "topleft",
    legend = var_lev,
    col    = pal,
    pch    = 19, pt.cex = 1.1, bty = "n", cex = 0.65, title = "variable"
  )
  graphics::legend(
    "bottomleft",
    legend = group_labels,
    col    = c(color_pos, color_neg),
    lwd    = 3, bty = "n", cex = 0.7, title = "edge direction"
  )
  
  invisible(g)
}
