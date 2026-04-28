#' Plot a modality graph
#'
#' Visualises a \code{catmodgraph} object. Nodes represent modalities and
#' edges represent cross-variable modality associations. By default, node
#' colours indicate the originating variable. If modality communities have
#' been detected, colours can instead reflect community membership.
#'
#' When \code{signed = TRUE}, edges are coloured by the sign of the stored
#' standardised Pearson residual (\code{std_resid} edge attribute): green
#' for positive (modalities co-occurring more than expected under
#' independence) and red for negative (co-occurring less than expected).
#' Edge alpha transparency then scales with \code{|std_resid|}.
#'
#' @param x A \code{catmodgraph} object.
#' @param color_by Character. One of \code{"variable"} (default) or
#'   \code{"cluster"}.
#' @param signed Logical. If \code{TRUE}, edge colour encodes the sign of
#'   the stored standardised Pearson residual: green for positive
#'   (attraction), red for negative (repulsion). Default is \code{FALSE},
#'   which uses a uniform grey edge colour. Requires the \code{std_resid}
#'   edge attribute (present by default in graphs built with
#'   \code{\link{build_modality_graph}}).
#' @param show_labels Logical. If \code{TRUE}, node labels are drawn.
#'   Default is \code{TRUE}.
#' @param layout Character. Graph layout passed to
#'   \code{igraph::layout_with_fr()} or \code{igraph::layout_with_kk()}.
#'   One of \code{"fr"} (default) or \code{"kk"}.
#' @param vertex_size Numeric. Node size. Default is \code{24}.
#' @param edge_scale Numeric. Multiplicative factor applied to edge
#'   widths. Default is \code{8}.
#' @param remove_isolates Logical. If \code{TRUE} (default), vertices with
#'   degree 0 are hidden from the plot. Isolated modalities typically
#'   arise after pruning and carry no community structure information,
#'   so removing them improves interpretability. Set to \code{FALSE} to
#'   see every vertex in the graph. Does not modify the input object.
#' @param ... Further arguments passed to \code{plot.igraph()}.
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' df <- expand_table(Titanic)
#' mg <- build_modality_graph(df)
#' mg <- cluster_modalities(mg)
#'
#' # Base plotting
#' plot(mg)
#'
#' # Colour nodes by detected modality community
#' plot(mg, color_by = "cluster")
#'
#' # Signed edges: green = attraction, red = repulsion
#' plot(mg, signed = TRUE)
#' 
#' @examples
#' df <- expand_table(Titanic)
#' mg <- build_modality_graph(df)
#' mg <- cluster_modalities(mg)
#'
#' # Base plotting
#' plot(mg)
#'
#' # Colour nodes by detected modality community
#' plot(mg, color_by = "cluster")
#'
#' # Signed edges: green = attraction, red = repulsion
#' plot(mg, signed = TRUE)
#'
#' # Show isolated vertices (not hidden by default)
#' plot(mg, remove_isolates = FALSE)
#'
#'
#' @importFrom graphics plot legend par
#' @importFrom grDevices hcl.colors adjustcolor
#' @importFrom igraph V E layout_with_fr layout_with_kk degree delete_vertices
#' @export
plot.catmodgraph <- function(x,
                             color_by        = c("variable", "cluster"),
                             signed          = FALSE,
                             show_labels     = TRUE,
                             layout          = c("fr", "kk"),
                             vertex_size     = 24,
                             edge_scale      = 8,
                             remove_isolates = TRUE,
                             ...) {
  
  if (!inherits(x, "catmodgraph")) {
    stop("`x` must be a catmodgraph object.", call. = FALSE)
  }
  
  color_by <- match.arg(color_by)
  layout <- match.arg(layout)
  
  if (!is.logical(signed) || length(signed) != 1L || is.na(signed)) {
    stop("`signed` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(show_labels) || length(show_labels) != 1L || is.na(show_labels)) {
    stop("`show_labels` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(remove_isolates) || length(remove_isolates) != 1L || is.na(remove_isolates)) {
    stop("`remove_isolates` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.numeric(vertex_size) || length(vertex_size) != 1L || is.na(vertex_size) || vertex_size <= 0) {
    stop("`vertex_size` must be a single positive number.", call. = FALSE)
  }
  
  if (!is.numeric(edge_scale) || length(edge_scale) != 1L || is.na(edge_scale) || edge_scale <= 0) {
    stop("`edge_scale` must be a single positive number.", call. = FALSE)
  }
  
  g <- x$graph
  
  if (igraph::vcount(g) == 0L) {
    stop("Cannot plot an empty modality graph.", call. = FALSE)
  }
  
  if (remove_isolates) {
    iso <- which(igraph::degree(g) == 0L)
    if (length(iso) > 0L) {
      g <- igraph::delete_vertices(g, iso)
    }
    if (igraph::vcount(g) == 0L) {
      stop(
        "All vertices are isolates; nothing to plot. Try `remove_isolates = FALSE` ",
        "to see the isolated vertices, or adjust the pruning threshold.",
        call. = FALSE
      )
    }
  }
  
  verts <- igraph::as_data_frame(g, what = "vertices")
  
  if (color_by == "cluster") {
    if (!"cluster" %in% names(verts)) {
      stop("`color_by = \"cluster\"` requires `cluster_modalities()` to be run first.", call. = FALSE)
    }
    grp <- as.factor(verts$cluster)
  } else {
    grp <- as.factor(verts$variable)
  }
  
  palette_vals <- grDevices::hcl.colors(n = nlevels(grp), palette = "Dark 3")
  vertex_cols  <- palette_vals[as.integer(grp)]
  
  # Edge widths always scale with absolute weight (phi).
  edge_w <- igraph::E(g)$weight
  edge_w <- ifelse(is.na(edge_w), 0, edge_w)
  edge_w <- pmax(edge_w, 0.01) * edge_scale
  
  # Edge colours: either signed (green/red by std_resid) or neutral grey.
  if (signed) {
    has_attr <- "std_resid" %in% igraph::edge_attr_names(g)
    std_res  <- if (has_attr) igraph::E(g)$std_resid else NULL
    
    if (!has_attr || length(std_res) == 0L || all(is.na(std_res))) {
      warning(
        "`signed = TRUE` requested but `std_resid` edge attribute is missing ",
        "or entirely NA; falling back to unsigned edges.",
        call. = FALSE
      )
      signed <- FALSE
    }
  }
  
  if (signed) {
    # Per-edge mapping: sign -> hue, |std_resid| -> alpha (clamped to [0, 1]).
    # adjustcolor() does not vectorise over alpha.f, so call per edge.
    edge_cols <- vapply(
      seq_along(std_res),
      function(i) {
        r <- std_res[i]
        if (is.na(r)) return(grDevices::adjustcolor("grey60", alpha.f = 0.4))
        base_col <- if (r > 0) "#1D9E75" else "#D85A30"
        a <- pmin(abs(r) / 10, 1)
        a <- pmax(a, 0.15)  # floor so very weak edges remain visible
        grDevices::adjustcolor(base_col, alpha.f = a)
      },
      FUN.VALUE = character(1)
    )
  } else {
    edge_alpha <- igraph::E(g)$weight
    edge_alpha <- ifelse(is.na(edge_alpha), 0.2, pmin(pmax(edge_alpha, 0.1), 1))
    edge_cols <- vapply(
      edge_alpha,
      function(a) grDevices::adjustcolor("grey40", alpha.f = a),
      FUN.VALUE = character(1)
    )
  }
  
  lay <- switch(
    layout,
    fr = igraph::layout_with_fr(g),
    kk = igraph::layout_with_kk(g)
  )
  
  plot(
    g,
    layout             = lay,
    vertex.color       = vertex_cols,
    vertex.size        = vertex_size,
    vertex.label       = if (show_labels) igraph::V(g)$name else NA,
    vertex.label.cex   = 0.8,
    vertex.label.color = "black",
    edge.width         = edge_w,
    edge.color         = edge_cols,
    ...
  )
  
  legend(
    "topleft",
    legend = levels(grp),
    col    = palette_vals,
    pch    = 19,
    pt.cex = 1.4,
    bty    = "n",
    title  = color_by
  )
  
  if (signed) {
    legend(
      "bottomleft",
      legend = c("positive std. residual (attraction)",
                 "negative std. residual (repulsion)"),
      col    = c("#1D9E75", "#D85A30"),
      lwd    = 3,
      bty    = "n",
      cex    = 0.75,
      title  = "edge sign"
    )
  }
  
  invisible(x)
}