#' Construct a categorical association network
#'
#' The primary user-facing constructor for \pkg{catgraph}. It computes
#' pairwise effect sizes (phi or Cramer's V) for all categorical variable
#' pairs, stores the resulting weighted \pkg{igraph} network, and preserves
#' processed data and metadata for downstream analysis. Use this function for
#' standard workflows; use \code{\link{build_graph}} only when a raw
#' \pkg{igraph} object is required.
#'
#' \strong{Scope.} A \code{catgraph} is a \emph{pairwise association network},
#' not a conditional-independence graphical model. Edges encode bivariate
#' dependence between two variables and do not imply that the two variables
#' remain dependent after controlling for the remaining variables. Interpret
#' centrality, community, and bridge measures accordingly. See the package
#' vignette for a full discussion.
#'
#' @param data A data frame or tibble whose columns represent categorical
#'   variables. Factor, character, and logical columns are supported.
#'   Numeric columns are coerced to character with a message.
#' @param corrected Logical. If \code{FALSE} (default), classical phi and
#'   Cramer's V are computed. If \code{TRUE}, the bias-corrected estimators
#'   of Bergsma (2013) are used.
#' @param correct Logical. Yates' continuity correction for the chi-square
#'   test. Default \code{FALSE}.
#' @param simulate_p Logical. Monte Carlo p-value simulation. Default
#'   \code{FALSE}.
#' @param B Integer. Monte Carlo resamples when \code{simulate_p = TRUE}.
#'   Default \code{2000L}.
#'
#' @return An S3 object of class \code{catgraph} containing:
#' \describe{
#'   \item{\code{graph}}{An undirected weighted \pkg{igraph} object. True
#'     zero associations are absent edges, not near-zero edges.}
#'   \item{\code{data}}{The \strong{processed} data frame actually used for
#'     estimation (after non-categorical coercion and constant-column
#'     removal). Downstream functions such as \code{\link{catgraph_ci}}
#'     resample from this object. Changed from \code{raw_data} in v0.4.0 to
#'     fix an internal-consistency bug.}
#'   \item{\code{raw_data}}{The original input data frame, for reference.}
#'   \item{\code{corrected}}{Logical flag indicating which estimator is active.}
#'   \item{\code{n_vars}}{Number of variables (graph vertices).}
#'   \item{\code{n_pairs_total}}{Number of variable pairs evaluated.}
#'   \item{\code{n_pairs}}{Number of retained graph edges (pairs with non-zero
#'     effect size).}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' All variable pairs with non-zero effect size are retained by default (no
#' thresholding at construction time). To remove weak or non-significant
#' edges, pass the object to \code{\link{prune_edges}}.
#'
#' @references
#' Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
#'   \emph{Journal of the Korean Statistical Society}, 42(3), 323--328.
#'   \doi{10.1016/j.jkss.2012.10.002}
#'
#' @examples
#' df <- expand_table(Titanic)
#' cg <- catgraph(df)
#' cg
#' summary(cg)
#'
#' cg_bc <- catgraph(df, corrected = TRUE)
#'
#' @seealso \code{\link{prune_edges}}, \code{\link{detect_clusters}},
#'   \code{\link{plot.catgraph}}, \code{\link{assoc_matrix}},
#'   \code{\link{assoc_similarity}}
#' @export
catgraph <- function(data,
                     corrected  = FALSE,
                     correct    = FALSE,
                     simulate_p = FALSE,
                     B          = 2000L) {
  
  mc <- match.call()
  
  g <- build_graph(
    data       = data,
    corrected  = corrected,
    correct    = correct,
    simulate_p = simulate_p,
    B          = B
  )
  
  processed <- attr(g, "processed_data")
  pair_results <- attr(g, "pair_results")
  attr(g, "processed_data") <- NULL
  attr(g, "pair_results") <- NULL
  
  structure(
    list(
      graph         = g,
      data          = processed,
      raw_data      = data,
      pair_results  = pair_results,
      corrected     = corrected,
      n_vars        = igraph::vcount(g),
      n_pairs_total = if (!is.null(pair_results)) nrow(pair_results) else NA_integer_,
      n_pairs       = igraph::ecount(g),
      call          = mc
    ),
    class = "catgraph"
  )
}


# ---------------------------------------------------------- S3: print ------

#' @describeIn catgraph Print a concise summary of a \code{catgraph} object.
#' @param x A \code{catgraph} object.
#' @param ... Ignored.
#' @importFrom stats median
#' @export
print.catgraph <- function(x, ...) {
  cat("catgraph object (pairwise association network)\n")
  cat("  Variables :", x$n_vars, "\n")
  cat("  Edges     :", x$n_pairs, "\n")
  estimator <- if (x$corrected) "bias-corrected (Bergsma 2013)" else "classical"
  cat("  Estimator :", estimator, "\n")
  
  w <- igraph::E(x$graph)$weight
  w_nna <- w[!is.na(w)]
  if (length(w_nna) > 0L) {
    cat(sprintf(
      "  Weights   : min = %.4f  median = %.4f  max = %.4f\n",
      min(w_nna), stats::median(w_nna), max(w_nna)
    ))
  }
  
  metrics <- igraph::E(x$graph)$metric
  if (!is.null(metrics) && length(metrics) > 0L) {
    tab <- table(metrics)
    parts <- paste(names(tab), tab, sep = " = ")
    cat("  Metric mix:", paste(parts, collapse = ", "), "\n")
  }
  
  cat("  Note      : edges encode pairwise marginal association, not\n")
  cat("              conditional independence. Edge weights use phi\n")
  cat("              (2x2) and Cramer's V (RxC); both lie on [0, 1],\n")
  cat("              but are not strictly exchangeable across table\n")
  cat("              dimensions. Interpret mixed-metric graphs with care.\n")
  cat("              See vignette 'Methodological caveats', item 2.\n")
  invisible(x)
}


# -------------------------------------------------------- S3: summary ------

#' @describeIn catgraph Summarise a \code{catgraph} object, listing edges
#'   sorted by effect size.
#' @param object A \code{catgraph} object.
#' @param top Integer. Number of strongest edges to display. Use \code{Inf}
#'   for all edges. Default \code{10L}.
#' @param ... Ignored.
#' @export
summary.catgraph <- function(object, top = 10L, ...) {
  if (!is.numeric(top) || length(top) != 1L || is.na(top) || top < 1) {
    stop("`top` must be a single positive number.", call. = FALSE)
  }
  
  g <- object$graph
  
  if (igraph::ecount(g) == 0L) {
    cat("catgraph summary: no edges (all pairs had zero effect size).\n")
    return(invisible(data.frame()))
  }
  
  el <- igraph::as_edgelist(g)
  tbl <- data.frame(
    var1        = el[, 1],
    var2        = el[, 2],
    effect_size = igraph::E(g)$weight,
    metric      = igraph::E(g)$metric,
    p_value     = igraph::E(g)$p_value,
    n           = igraph::E(g)$n,
    type        = igraph::E(g)$type,
    stringsAsFactors = FALSE
  )
  
  tbl <- tbl[order(tbl$effect_size, decreasing = TRUE), ]
  rownames(tbl) <- NULL
  n_show <- min(top, nrow(tbl))
  
  cat("catgraph summary\n")
  cat("  Variables       :", object$n_vars, "\n")
  if (!is.null(object$n_pairs_total)) {
    cat("  Pairs evaluated :", object$n_pairs_total, "\n")
  }
  cat("  Edges retained  :", object$n_pairs, "\n\n")
  
  estimator <- if (object$corrected) "bias-corrected (Bergsma 2013)" else "classical"
  cat("  Estimator       :", estimator, "\n\n")
  cat(sprintf("  Top %d edges by effect size:\n\n", n_show))
  print(tbl[seq_len(n_show), ], digits = 4)
  
  invisible(tbl)
}


# ---------------------------------------------------------- S3: plot -------

#' Plot a \code{catgraph} object
#'
#' Visualises the undirected weighted association network. The default
#' renderer uses \pkg{igraph}'s base-graphics plot. If \code{engine = "ggraph"}
#' is requested and \pkg{ggraph} is installed, a \pkg{ggplot2}-based plot is
#' returned as a \code{ggplot} object.
#'
#' @param x A \code{catgraph} object.
#' @param engine Character. Either \code{"igraph"} (default) or \code{"ggraph"}.
#' @param layout Character. Layout algorithm. For \code{engine = "igraph"},
#'   one of \code{"fr"} (Fruchterman-Reingold, default), \code{"kk"}
#'   (Kamada-Kawai), \code{"circle"}, \code{"grid"}, \code{"graphopt"},
#'   \code{"nicely"}, or \code{"random"}. For \code{engine = "ggraph"},
#'   any layout string accepted by \code{ggraph::ggraph}.
#' @param edge_width_range Numeric vector of length 2. Minimum and maximum
#'   line widths mapped to edge weights. Default \code{c(0.5, 4)}.
#' @param vertex_size Numeric. Vertex size for \code{engine = "igraph"}.
#'   Default \code{30}.
#' @param vertex_color Character. Vertex fill colour.
#' @param edge_color Character. Edge colour.
#' @param label_size Numeric. Label character expansion. Default \code{0.8}.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param ... Additional arguments passed to the underlying renderer.
#'
#' @return For \code{engine = "igraph"}: invisibly returns \code{x} (called
#'   for its side effect).
#'   For \code{engine = "ggraph"}: a \code{ggplot} object.
#'
#' @details
#' In v0.3.0 and earlier, the \code{layout} argument was silently ignored by
#' the igraph branch (Fruchterman-Reingold was always used). From 0.4.0
#' onwards, \code{layout} is respected for both engines.
#'
#' @examples
#' df <- expand_table(Titanic)
#' cg <- catgraph(df)
#' plot(cg)
#' plot(cg, layout = "kk")
#' plot(cg, layout = "circle")
#'
#' @importFrom igraph plot.igraph E V
#' @export
plot.catgraph <- function(x,
                          engine           = c("igraph", "ggraph"),
                          layout           = "fr",
                          edge_width_range = c(0.5, 4),
                          vertex_size      = 30,
                          vertex_color     = "#AFA9EC",
                          edge_color       = "#888780",
                          label_size       = 0.8,
                          title            = NULL,
                          ...) {
  engine <- match.arg(engine)

  g <- x$graph
  if (igraph::ecount(g) == 0L) {
    warning("catgraph has no edges to plot.")
  }

  w <- igraph::E(g)$weight
  w[is.na(w)] <- 0

  # Rescale weights to edge widths (guard against single-edge / constant w)
  if (length(w) == 0L) {
    scaled_w <- numeric(0)
  } else {
    w_range <- range(w, na.rm = TRUE)
    if (diff(w_range) < .Machine$double.eps) {
      scaled_w <- rep(mean(edge_width_range), length(w))
    } else {
      scaled_w <- edge_width_range[1] +
        (w - w_range[1]) / diff(w_range) * diff(edge_width_range)
    }
  }

  if (engine == "igraph") {
    # v0.4.0 FIX: honour the layout argument instead of hard-coding FR.
    lay <- .catgraph_igraph_layout(g, layout)

    igraph::plot.igraph(
      g,
      layout             = lay,
      vertex.size        = vertex_size,
      vertex.color       = vertex_color,
      vertex.label.cex   = label_size,
      vertex.label.color = "#2C2C2A",
      edge.width         = scaled_w,
      edge.color         = edge_color,
      main               = title,
      ...
    )
    return(invisible(x))
  }

  # ------------------------------------------------------------- ggraph branch
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Packages 'ggraph' and 'ggplot2' are needed for engine = 'ggraph'. ",
      "Install them with: install.packages(c('ggraph', 'ggplot2'))"
    )
  }

  igraph::E(g)$scaled_width <- scaled_w

  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = .data$scaled_width,
                   alpha = .data$weight),
      colour = edge_color,
      show.legend = FALSE
    ) +
    ggraph::scale_edge_width(range = edge_width_range) +
    ggraph::geom_node_point(size = vertex_size / 5, color = vertex_color) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$name),
      repel = TRUE, size = label_size * 4
    ) +
    ggplot2::theme_void()

  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}


# Internal helper: map layout string -> igraph layout matrix
.catgraph_igraph_layout <- function(g, layout) {
  layout <- tolower(as.character(layout)[[1L]])
  fn <- switch(
    layout,
    fr       = igraph::layout_with_fr,
    kk       = igraph::layout_with_kk,
    circle   = igraph::layout_in_circle,
    grid     = igraph::layout_on_grid,
    graphopt = igraph::layout_with_graphopt,
    nicely   = igraph::layout_nicely,
    random   = igraph::layout_randomly,
    NULL
  )
  if (is.null(fn)) {
    stop("Unknown `layout`: '", layout,
         "'. Valid options: fr, kk, circle, grid, graphopt, nicely, random.")
  }
  fn(g)
}


# -------------------------------------------------- S3: as_igraph ----------

#' Extract the underlying igraph object from a catgraph
#'
#' @param x A \code{catgraph} object.
#' @param ... Ignored.
#' @return An \pkg{igraph} undirected graph.
#' @export
as_igraph <- function(x, ...) UseMethod("as_igraph")

#' @export
as_igraph.catgraph <- function(x, ...) x$graph


# -------------------------------------------------- as_ggraph --------------

#' Coerce a catgraph to a ggraph-compatible tbl_graph
#'
#' @param x A \code{catgraph} object.
#' @param ... Ignored.
#' @return A \code{tbl_graph} object (requires \pkg{tidygraph}).
#' @export
as_ggraph <- function(x, ...) UseMethod("as_ggraph")

#' @export
as_ggraph.catgraph <- function(x, ...) {
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    stop(
      "Package 'tidygraph' is required for as_ggraph(). ",
      "Install with: install.packages('tidygraph')"
    )
  }
  tidygraph::as_tbl_graph(x$graph)
}
