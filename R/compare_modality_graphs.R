#' Compare multiple modality networks on one panel
#'
#' Visualises two or more \code{catmodgraph} objects side-by-side with a
#' shared node layout, for exploring how category-level marginal association
#' structure differs across populations, sites, or time points.
#'
#' @param x A named list of \code{catmodgraph} objects. The list names
#'   become panel titles. A minimum of two graphs is required.
#' @param restrict Character. How to handle modalities that appear in
#'   only some of the input graphs. One of:
#'   \describe{
#'     \item{\code{"common"} (default)}{Restrict every panel to
#'       modalities present in all input graphs. Safest for
#'       comparison.}
#'     \item{\code{"union"}}{Use the union of modalities across all
#'       inputs. Modalities absent from a given graph appear as
#'       isolated vertices in that panel.}
#'   }
#' @param pruning Character. How to filter edges before plotting. One of
#'   \code{"individual"} (default), \code{"none"}. Pooled pruning is
#'   not supported for modality networks because it requires pooling
#'   the underlying row data and re-estimating, which may change the
#'   modality set non-trivially. Use \code{"individual"} with a shared
#'   \code{min_weight} to ensure comparable per-panel thresholds.
#' @param min_weight Numeric. Minimum phi threshold for
#'   \code{"individual"} pruning. Default \code{0.10}.
#' @param max_p Numeric. Maximum p-value for \code{"individual"} pruning.
#'   Default \code{0.05}.
#' @param signed Logical. If \code{TRUE}, edges are coloured by sign of
#'   the stored \code{std_resid} attribute (green = attraction, red =
#'   repulsion). Default \code{FALSE}.
#' @param layout_fn Function. An \pkg{igraph} layout function applied to
#'   the union graph. Default \code{igraph::layout_with_fr}.
#' @param vertex_size Numeric. Vertex size. Default \code{14}.
#' @param edge_scale Numeric. Multiplier for edge widths. Default
#'   \code{6}.
#' @param ... Further arguments passed to \code{\link[igraph]{plot.igraph}}.
#'
#' @return Invisibly returns the reference union graph.
#'
#' @details
#' Modality graphs are inherently harder to compare than variable-level
#' graphs because modality sets may differ across groups (e.g., a
#' response category endorsed in population A but not population B).
#' The \code{restrict = "common"} default avoids this by reducing all
#' panels to a shared vocabulary; \code{"union"} preserves all nodes
#' but panels will differ in vertex presence.
#'
#' @section Formal testing:
#' For inferential comparison of modality networks, see
#' \code{\link{test_modality_graph_equality}} (omnibus) and
#' \code{\link{test_modality_edge_differences}} (edge-wise post-hoc).
#'
#' @examples
#' # Build two modality graphs from subsets of HairEyeColor
#' df <- expand_table(HairEyeColor)
#' mg_f <- build_modality_graph(df[df$Sex == "Female", c("Hair", "Eye")])
#' mg_m <- build_modality_graph(df[df$Sex == "Male",   c("Hair", "Eye")])
#'
#' compare_modality_graphs(list(Female = mg_f, Male = mg_m))
#'
#' @seealso \code{\link{build_modality_graph}}, \code{\link{compare_catgraphs}}
#' @importFrom graphics par legend
#' @importFrom grDevices hcl.colors adjustcolor
#' @importFrom igraph layout_with_fr V E vcount ecount as_edgelist
#'   induced_subgraph add_vertices add_edges simplify edge_attr_names
#' @importFrom stats setNames
#' @importFrom scales rescale
#' @export
compare_modality_graphs <- function(x,
                                    restrict      = c("common", "union"),
                                    pruning       = c("individual", "none"),
                                    min_weight    = 0.10,
                                    max_p         = 0.05,
                                    signed        = FALSE,
                                    layout_fn     = igraph::layout_with_fr,
                                    vertex_size   = 14,
                                    edge_scale    = 6,
                                    ...) {
  
  .validate_compare_input(x, class_name = "catmodgraph")
  
  restrict <- match.arg(restrict)
  pruning  <- match.arg(pruning)
  
  if (!is.logical(signed) || length(signed) != 1L || is.na(signed)) {
    stop("`signed` must be TRUE or FALSE.", call. = FALSE)
  }
  
  # Apply per-panel pruning if requested
  if (pruning == "individual") {
    x <- lapply(x, function(mg) {
      prune_modality_edges(mg, min_weight = min_weight, max_p = max_p)
    })
  }
  
  # Resolve vertex set across panels
  vertex_sets <- lapply(x, function(mg) igraph::V(mg$graph)$name)
  if (restrict == "common") {
    common <- Reduce(intersect, vertex_sets)
    if (length(common) < 2L) {
      stop("Fewer than 2 modalities are common across all inputs; ",
           "try restrict = 'union' or relax pruning.", call. = FALSE)
    }
    panel_graphs <- lapply(x, function(mg) {
      keep <- igraph::V(mg$graph)$name %in% common
      igraph::induced_subgraph(mg$graph, vids = which(keep))
    })
  } else {
    all_v <- Reduce(union, vertex_sets)
    panel_graphs <- lapply(x, function(mg) {
      missing_v <- setdiff(all_v, igraph::V(mg$graph)$name)
      g <- mg$graph
      if (length(missing_v) > 0L) {
        g <- igraph::add_vertices(g, length(missing_v), name = missing_v)
      }
      g
    })
  }
  
  # Shared layout on union graph, rescaled to [-1, 1] for consistent panels
  union_g <- .union_graph(panel_graphs)
  lay     <- .shared_layout(union_g, layout_fn, rescale_to = c(-1, 1))
  
  # Colour modality nodes by variable
  mod_vars   <- sub("=.*", "", igraph::V(union_g)$name)
  var_levels <- sort(unique(mod_vars))
  pal        <- grDevices::hcl.colors(length(var_levels), palette = "Dark 3")
  col_map    <- stats::setNames(pal, var_levels)
  
  # ------- Plot multi-panel layout -----------------------------------------
  n <- length(panel_graphs)
  
  # Tighter margins when > 3 panels; also scale label size
  tight <- n > 3L
  op <- graphics::par(
    mfrow = c(1, n),
    mar   = if (tight) c(1, 0.5, 2.5, 0.5) else c(2, 2, 3, 2),
    oma   = c(0, 0, 0, 0)
  )
  on.exit(graphics::par(op), add = TRUE)
  
  label_cex <- if (tight) 0.55 else 0.65
  
  for (i in seq_len(n)) {
    g <- panel_graphs[[i]]
    v_mod <- sub("=.*", "", igraph::V(g)$name)
    v_col <- col_map[v_mod]
    
    if (signed && "std_resid" %in% igraph::edge_attr_names(g) &&
        igraph::ecount(g) > 0L) {
      std_res <- igraph::E(g)$std_resid
      edge_cols <- vapply(
        seq_along(std_res),
        function(j) {
          r <- std_res[j]
          if (is.na(r)) return(grDevices::adjustcolor("grey60", alpha.f = 0.4))
          base_col <- if (r > 0) "#1D9E75" else "#D85A30"
          a <- pmax(pmin(abs(r) / 10, 1), 0.15)
          grDevices::adjustcolor(base_col, alpha.f = a)
        },
        FUN.VALUE = character(1)
      )
    } else {
      edge_cols <- grDevices::adjustcolor("grey40", alpha.f = 0.6)
    }
    
    edge_w <- if (igraph::ecount(g) > 0L) {
      pmax(igraph::E(g)$weight, 0.01) * edge_scale
    } else numeric(0)
    
    plot(
      g,
      layout              = .align_layout(g, lay),
      vertex.color        = v_col,
      vertex.size         = vertex_size,
      vertex.label        = igraph::V(g)$name,
      vertex.label.cex    = label_cex,
      vertex.label.color  = "black",
      vertex.label.dist   = if (tight) 0.8 else 0,
      vertex.label.degree = -pi/2,   # labels below vertices
      edge.width          = edge_w,
      edge.color          = edge_cols,
      main                = sprintf("%s (%d edges)",
                                    names(panel_graphs)[i],
                                    igraph::ecount(g)),
      ...
    )
    
    if (i == 1L) {
      legend("topleft", legend = var_levels, col = pal,
             pch = 19, pt.cex = 1.1, bty = "n", cex = 0.65,
             title = "variable")
    }
  }
  
  invisible(union_g)
}