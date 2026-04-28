#' Extract pairwise association weights as a matrix or tidy data frame
#'
#' Returns the pairwise effect size values from a \code{catgraph} object as
#' a symmetric matrix or a long-format tidy data frame suitable for further
#' analysis or export.
#'
#' @param x A \code{catgraph} object.
#' @param format Character. Output format: \code{"matrix"} (default) returns
#'   a symmetric numeric matrix with \code{NA} on the diagonal;
#'   \code{"tidy"} returns a data frame with one row per pair.
#' @param include_p Logical. When \code{format = "tidy"}, whether to include
#'   the p-value column. Default \code{TRUE}.
#' @param include_n Logical. When \code{format = "tidy"}, whether to include
#'   the pairwise observation count. Default \code{TRUE}.
#'
#' @return
#' \describe{
#'   \item{\code{format = "matrix"}}{A symmetric numeric matrix of dimension
#'     p x p (p = number of variables). Diagonal is \code{NA}. Row and column
#'     names are the variable names.}
#'   \item{\code{format = "tidy"}}{A data frame with columns \code{var1},
#'     \code{var2}, \code{effect_size}, \code{metric}, \code{type},
#'     and optionally \code{p_value} and \code{n}.}
#' }
#'
#' @examples
#' df <- as.data.frame(Titanic)
#' df_exp <- df[rep(seq_len(nrow(df)), df$Freq), -5]
#' cg <- catgraph(df_exp)
#' assoc_matrix(cg)
#' assoc_matrix(cg, format = "tidy")
#'
#' @seealso \code{\link{catgraph}}, \code{\link{prune_edges}}
#' @importFrom igraph as_edgelist E vcount V
#' @export
assoc_matrix <- function(x,
                         format    = c("matrix", "tidy"),
                         include_p = TRUE,
                         include_n = TRUE) {
  
  if (!inherits(x, "catgraph")) {
    stop("`x` must be a catgraph object.", call. = FALSE)
  }
  
  if (is.null(x$graph) || !igraph::is_igraph(x$graph)) {
    stop("`x$graph` must be a valid igraph object.", call. = FALSE)
  }
  
  format <- match.arg(format)
  g <- x$graph
  vars <- igraph::V(g)$name
  
  if (length(vars) == 0L) {
    stop("The graph has no vertices.", call. = FALSE)
  }
  
  if (igraph::ecount(g) == 0L) {
    if (format == "tidy") {
      out <- data.frame(
        var1 = character(0),
        var2 = character(0),
        effect_size = numeric(0),
        metric = character(0),
        type = character(0),
        stringsAsFactors = FALSE
      )
      if (include_p) out$p_value <- numeric(0)
      if (include_n) out$n <- numeric(0)
      return(out)
    }
    
    p <- length(vars)
    return(matrix(
      NA_real_,
      nrow = p,
      ncol = p,
      dimnames = list(vars, vars)
    ))
  }
  
  required_attrs <- c("weight", "metric", "type")
  edge_attrs <- igraph::edge_attr_names(g)
  missing_required <- setdiff(required_attrs, edge_attrs)
  
  if (length(missing_required) > 0L) {
    stop(
      sprintf(
        "The graph is missing required edge attribute(s): %s.",
        paste(missing_required, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  has_p <- "p_value" %in% edge_attrs
  has_n <- "n" %in% edge_attrs
  
  el <- igraph::as_edgelist(g, names = TRUE)
  
  if (format == "tidy") {
    tbl <- data.frame(
      var1 = el[, 1],
      var2 = el[, 2],
      effect_size = igraph::E(g)$weight,
      metric = igraph::E(g)$metric,
      type = igraph::E(g)$type,
      stringsAsFactors = FALSE
    )
    
    if (include_p) tbl$p_value <- if (has_p) igraph::E(g)$p_value else NA_real_
    if (include_n) tbl$n <- if (has_n) igraph::E(g)$n else NA_real_
    
    ord <- order(tbl$effect_size, decreasing = TRUE, na.last = TRUE)
    return(tbl[ord, , drop = FALSE])
  }
  
  p <- length(vars)
  mat <- matrix(
    NA_real_,
    nrow = p,
    ncol = p,
    dimnames = list(vars, vars)
  )
  
  w <- igraph::E(g)$weight
  for (i in seq_len(nrow(el))) {
    mat[el[i, 1], el[i, 2]] <- w[i]
    mat[el[i, 2], el[i, 1]] <- w[i]
  }
  
  mat
}