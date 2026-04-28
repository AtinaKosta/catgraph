#' Construct a bipartite respondent-modality graph
#'
#' Builds a two-mode (bipartite) graph where one partition of vertices
#' represents respondents (rows of the data) and the other represents
#' modalities (factor levels). An edge connects respondent \eqn{i} to
#' modality \eqn{V{=}l} if row \eqn{i} has value \eqn{l} on variable \eqn{V}.
#'
#' Bipartite (two-mode) graphs preserve the full respondent-to-modality
#' incidence structure. They are the unprojected counterpart of the
#' modality co-association graph returned by
#' \code{\link{build_modality_graph}}: projecting a \code{catbipartite}
#' onto its modality partition yields a respondents-in-common modality
#' network, which is a \strong{count-based co-occurrence network}
#' distinct from the chi-square / phi-filtered network produced by
#' \code{build_modality_graph}.
#'
#' The bipartite graph is the correct tool when the scientific object of
#' interest is the raw incidence relationship between units and categories
#' — for example, affiliation networks (persons to events), ecological
#' surveys (sites to species), or survey data seen as respondents
#' endorsing response levels.
#'
#' @param data A data frame of categorical variables.
#' @param remove_na Logical. If \code{TRUE} (default), rows with any
#'   missing value are dropped before graph construction.
#' @param row_prefix Character. Prefix used to name respondent-side
#'   vertices. Default \code{"r"}, so respondents are labelled
#'   \code{"r1", "r2", ...}.
#'
#' @return An object of class \code{"catbipartite"}, a list with:
#'   \describe{
#'     \item{\code{graph}}{An \code{igraph} undirected graph with vertex
#'       attribute \code{type} (FALSE = respondent, TRUE = modality).}
#'     \item{\code{modalities}}{Data frame describing modality-side
#'       vertices, with columns \code{node}, \code{variable},
#'       \code{modality}.}
#'     \item{\code{n_rows}}{Number of respondent-side vertices.}
#'     \item{\code{n_modalities}}{Number of modality-side vertices.}
#'     \item{\code{data}}{The processed categorical data used to build
#'       the graph (after NA removal).}
#'   }
#'
#' @details
#' The resulting graph satisfies the \pkg{igraph} convention for bipartite
#' graphs: \code{igraph::bipartite_projection(g)} can be called directly
#' to obtain either the respondent-side or the modality-side projection.
#'
#' \strong{Scope.} A \code{catbipartite} is a raw incidence graph; it
#' applies no statistical filter. Edges reflect endorsement, not
#' association. For statistically-filtered category co-association,
#' see \code{\link{build_modality_graph}}.
#'
#' @examples
#' df <- expand_table(Titanic)
#' bg <- bipartite_modality_graph(df)
#' bg
#' summary(bg)
#'
#' # Project onto modalities: edges weighted by shared respondents
#' proj <- igraph::bipartite_projection(bg$graph, which = "true")
#'
#' @seealso \code{\link{build_modality_graph}} for the chi-square /
#'   phi-filtered modality co-association graph.
#' @importFrom igraph graph_from_data_frame V
#' @export
bipartite_modality_graph <- function(data,
                                     remove_na  = TRUE,
                                     row_prefix = "r") {
  
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  
  if (ncol(data) < 2L) {
    stop("`data` must contain at least two variables.", call. = FALSE)
  }
  
  if (!is.logical(remove_na) || length(remove_na) != 1L || is.na(remove_na)) {
    stop("`remove_na` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.character(row_prefix) || length(row_prefix) != 1L ||
      is.na(row_prefix) || nchar(row_prefix) == 0L) {
    stop("`row_prefix` must be a single non-empty character string.",
         call. = FALSE)
  }
  
  dat <- data
  if (remove_na) {
    dat <- stats::na.omit(dat)
  }
  if (nrow(dat) == 0L) {
    stop("No complete rows remain after NA removal.", call. = FALSE)
  }
  
  dat[] <- lapply(dat, function(x) if (is.factor(x)) x else factor(x))
  
  var_names <- names(dat)
  n_rows    <- nrow(dat)
  
  modalities_list <- lapply(var_names, function(v) {
    lv <- levels(dat[[v]])
    data.frame(
      node     = paste0(v, "=", lv),
      variable = v,
      modality = lv,
      stringsAsFactors = FALSE
    )
  })
  modalities <- do.call(rbind, modalities_list)
  
  row_nodes <- paste0(row_prefix, seq_len(n_rows))
  
  edges_per_row <- lapply(seq_len(n_rows), function(i) {
    data.frame(
      from = row_nodes[i],
      to   = vapply(var_names, function(v) {
        paste0(v, "=", as.character(dat[[v]][i]))
      }, character(1L)),
      stringsAsFactors = FALSE
    )
  })
  edge_df <- do.call(rbind, edges_per_row)
  
  verts <- data.frame(
    name = c(row_nodes, modalities$node),
    type = c(rep(FALSE, n_rows), rep(TRUE, nrow(modalities))),
    stringsAsFactors = FALSE
  )
  
  g <- igraph::graph_from_data_frame(
    d        = edge_df,
    directed = FALSE,
    vertices = verts
  )
  
  out <- list(
    graph        = g,
    modalities   = modalities,
    n_rows       = n_rows,
    n_modalities = nrow(modalities),
    data         = dat
  )
  class(out) <- "catbipartite"
  out
}


#' @export
print.catbipartite <- function(x, ...) {
  cat("catbipartite object (respondent-modality incidence graph)\n")
  cat("  Respondent nodes :", x$n_rows, "\n")
  cat("  Modality nodes   :", x$n_modalities, "\n")
  cat("  Edges            :", igraph::ecount(x$graph), "\n")
  cat("  Variables        :",
      paste(unique(x$modalities$variable), collapse = ", "), "\n")
  invisible(x)
}


#' @export
summary.catbipartite <- function(object, ...) {
  list(
    n_rows                  = object$n_rows,
    n_modalities            = object$n_modalities,
    n_edges                 = igraph::ecount(object$graph),
    variables               = unique(object$modalities$variable),
    modalities_per_variable = table(object$modalities$variable)
  )
}


#' Plot a bipartite respondent-modality graph
#'
#' Visualises a \code{catbipartite} object with respondents on one side
#' and modalities on the other. Respondent nodes are plotted small and
#' unlabelled; modality nodes are plotted large and coloured by
#' originating variable.
#'
#' @param x A \code{catbipartite} object.
#' @param show_respondents Logical. If \code{TRUE} (default),
#'   respondent-side vertices are drawn as small dots. If \code{FALSE},
#'   only the modality partition is plotted (useful when the number of
#'   respondents would make the plot unreadable).
#' @param max_respondents Integer or \code{NULL}. If the number of
#'   respondents exceeds this value and \code{show_respondents = TRUE},
#'   a random sample of \code{max_respondents} is drawn. Set to
#'   \code{NULL} to disable. Default \code{500}.
#' @param vertex_size_mod Numeric. Vertex size for modality nodes.
#'   Default \code{18}.
#' @param vertex_size_row Numeric. Vertex size for respondent nodes.
#'   Default \code{2}.
#' @param ... Further arguments passed to \code{plot.igraph()}.
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' df <- expand_table(Titanic)
#' bg <- bipartite_modality_graph(df)
#' plot(bg)
#' plot(bg, show_respondents = FALSE)
#'
#' @importFrom graphics plot legend
#' @importFrom grDevices hcl.colors adjustcolor
#' @importFrom igraph V E vcount ecount layout_as_bipartite induced_subgraph
#' @importFrom graphics plot legend
#' @importFrom grDevices hcl.colors adjustcolor
#' @importFrom igraph V E vcount ecount layout_as_bipartite induced_subgraph
#' @importFrom stats setNames
#' @export
plot.catbipartite <- function(x,
                              show_respondents = TRUE,
                              max_respondents  = 500L,
                              vertex_size_mod  = 18,
                              vertex_size_row  = 2,
                              ...) {
  
  if (!inherits(x, "catbipartite")) {
    stop("`x` must be a catbipartite object.", call. = FALSE)
  }
  
  if (!is.logical(show_respondents) || length(show_respondents) != 1L ||
      is.na(show_respondents)) {
    stop("`show_respondents` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.null(max_respondents) &&
      (!is.numeric(max_respondents) || length(max_respondents) != 1L ||
       is.na(max_respondents) || max_respondents < 1)) {
    stop("`max_respondents` must be NULL or a single positive integer.",
         call. = FALSE)
  }
  
  g <- x$graph
  vtype <- igraph::V(g)$type
  
  if (!show_respondents) {
    keep <- which(vtype)
    g <- igraph::induced_subgraph(g, vids = keep)
  } else if (!is.null(max_respondents) &&
             x$n_rows > as.integer(max_respondents)) {
    row_vids <- which(!vtype)
    sampled  <- sample(row_vids, as.integer(max_respondents))
    keep     <- sort(c(sampled, which(vtype)))
    g        <- igraph::induced_subgraph(g, vids = keep)
  }
  
  vtype_local <- igraph::V(g)$type
  
  mod_vars_full <- sub("=.*", "", igraph::V(g)$name)
  var_levels    <- sort(unique(mod_vars_full[vtype_local]))
  pal           <- grDevices::hcl.colors(length(var_levels), palette = "Dark 3")
  mod_col_map   <- setNames(pal, var_levels)
  
  vertex_cols <- ifelse(
    vtype_local,
    mod_col_map[mod_vars_full],
    grDevices::adjustcolor("grey70", alpha.f = 0.4)
  )
  
  vertex_sizes  <- ifelse(vtype_local, vertex_size_mod, vertex_size_row)
  vertex_labels <- ifelse(vtype_local, igraph::V(g)$name, NA)
  
  lay <- igraph::layout_as_bipartite(g, types = vtype_local)
  
  plot(
    g,
    layout             = lay,
    vertex.color       = vertex_cols,
    vertex.size        = vertex_sizes,
    vertex.label       = vertex_labels,
    vertex.label.cex   = 0.7,
    vertex.label.color = "black",
    edge.width         = 0.3,
    edge.color         = grDevices::adjustcolor("grey50", alpha.f = 0.2),
    ...
  )
  
  legend(
    "topright",
    legend = var_levels,
    col    = pal,
    pch    = 19,
    pt.cex = 1.2,
    bty    = "n",
    cex    = 0.75,
    title  = "variable"
  )
  
  invisible(x)
}