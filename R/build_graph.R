#' Build the underlying igraph association network
#'
#' Computes pairwise effect sizes (phi or Cramer's V) for all pairs of
#' categorical columns in a data frame and returns the underlying
#' \pkg{igraph} object. This is the lower-level computational engine used by
#' \code{\link{catgraph}}. Most users should call \code{\link{catgraph}}
#' unless they specifically need direct access to the raw \pkg{igraph}
#' representation.
#'
#' Computes effect sizes (phi or Cramer's V) for all pairs of categorical
#' columns in a data frame and returns an \pkg{igraph} object whose edge
#' weights correspond to those effect sizes. This is the main computational
#' engine used by the top-level \code{\link{catgraph}} constructor.
#'
#' \strong{Scope.} The returned graph represents \emph{pairwise marginal}
#' association strength. It is not a conditional-independence graphical model:
#' an edge between \code{A} and \code{B} does not imply that the variables
#' remain dependent after controlling for the other variables in the data.
#' See the package vignette section "Scope and interpretation" for details.
#'
#' @param data A data frame or tibble. All columns are treated as categorical.
#'   Non-factor, non-character, non-logical columns are coerced to character
#'   with a message. Columns with only one unique observed value (after
#'   pairwise deletion) are dropped with a warning.
#' @param method Character. Association metric used to weight edges. One of:
#'   \code{"cramers_v"} (default, classical phi / Cramer's V),
#'   \code{"cramers_v_corrected"} (bias-corrected via Bergsma 2013),
#'   \code{"nmi"} (Normalised Mutual Information),
#'   \code{"ami"} (Adjusted Mutual Information, corrects NMI for chance), or
#'   \code{"bayesian_cramers_v"} (Dirichlet-smoothed Cramér's V).
#' @param alpha Numeric. Dirichlet prior concentration for
#'   \code{method = "bayesian_cramers_v"}. Default \code{0.5} (Jeffreys
#'   prior). Ignored for all other methods.
#' @param corrected Logical. Deprecated shortcut: if \code{TRUE}, overrides
#'   \code{method} to \code{"cramers_v_corrected"}. Kept for backward
#'   compatibility. Default \code{FALSE}.
#' @param correct Logical. Yates' continuity correction for chi-square.
#'   Default \code{FALSE}.
#' @param simulate_p Logical. Use Monte Carlo simulation for p-values.
#'   Default \code{FALSE}.
#' @param B Integer. Number of Monte Carlo resamples. Default \code{2000L}.
#'
#' @return An \pkg{igraph} undirected graph. Pairs with a true zero effect
#'   size (no association whatsoever) are represented as \emph{absent edges}
#'   rather than near-zero edges, so the graph is sparse rather than
#'   structurally complete. The attribute \code{"processed_data"} on the
#'   returned graph holds the data frame actually used for estimation (after
#'   coercion and constant-column removal), which downstream functions such
#'   as \code{\link{catgraph_ci}} use when resampling.
#'   The graph attribute \code{"pair_results"} stores the full pairwise
#'   results table before zero-weight edges are omitted.
#'   
#' For ordinary package use, prefer \code{\link{catgraph}}, which wraps this
#' graph together with processed data, metadata, and S3 methods.
#'
#'   Vertex and edge attributes:
#' \describe{
#'   \item{Vertices}{One per column in \code{data}, with vertex attribute
#'     \code{name} set to the column name. Isolated vertices are preserved
#'     even if all their pairs have zero effect size.}
#'   \item{Edge attribute \code{weight}}{The phi or Cramer's V value.}
#'   \item{Edge attribute \code{metric}}{\code{"phi"} or \code{"cramers_v"}.}
#'   \item{Edge attribute \code{corrected}}{Whether bias correction was applied.}
#'   \item{Edge attribute \code{p_value}}{Chi-square p-value.}
#'   \item{Edge attribute \code{statistic}}{Chi-square statistic.}
#'   \item{Edge attribute \code{df}}{Degrees of freedom.}
#'   \item{Edge attribute \code{n}}{Pairwise-complete observation count for
#'     that pair. Values can differ across edges when missingness is present
#'     (pairwise deletion).}
#'   \item{Edge attribute \code{type}}{\code{"2x2"} or \code{"RxC"}.}
#'   \item{Edge attribute \code{estimable}}{Logical indicating whether the
#'     pairwise effect size was estimable before zero-weight omission.}
#' }
#'
#' @details
#' All variable pairs with non-zero effect size are included by default. Use
#' \code{\link{prune_edges}} to remove edges below a weight or adjusted-p
#' threshold after construction.
#'
#' \strong{Note on zero-weight pairs.} In earlier versions of the package
#' (<= 0.3.0), zero-weight pairs were stored as edges with weight
#' \code{.Machine$double.eps} to guarantee a fully connected graph. This
#' made the graph structurally complete and silently inflated density-based
#' measures. From 0.4.0 onwards, zero-weight pairs are absent edges; a dense
#' similarity matrix suitable for heatmaps is available separately via
#' \code{\link{assoc_similarity}} or \code{\link{assoc_matrix}}.
#'
#' @references
#' Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
#'   \emph{Journal of the Korean Statistical Society}, 42(3), 323--328.
#'   \doi{10.1016/j.jkss.2012.10.002}
#'
#' Good, I. J. (1965). \emph{The Estimation of Probabilities: An Essay on
#'   Modern Bayesian Methods}. MIT Press.
#'
#' Cover, T. M., & Thomas, J. A. (2006). \emph{Elements of Information
#'   Theory} (2nd ed.). Wiley. \doi{10.1002/047174882X}
#'
#' Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic
#'   measures for clusterings comparison: Variants, properties,
#'   normalisation and correction for chance.
#'   \emph{Journal of Machine Learning Research}, 11, 2837--2854.
#'   \url{https://jmlr.org/papers/v11/vinh10a.html}
#'
#' Csardi, G., & Nepusz, T. (2006). The igraph software package for complex
#'   network research. \emph{InterJournal, Complex Systems}, 1695.
#'   \url{https://igraph.org}
#'
#' @examples
#' data(HairEyeColor)
#' df <- expand_table(HairEyeColor)
#' g  <- build_graph(df[, c("Hair", "Eye")])
#' igraph::E(g)$weight
#'
#' @seealso \code{\link{catgraph}}, \code{\link{prune_edges}},
#'   \code{\link{effect_size}}, \code{\link{assoc_similarity}}
#' @importFrom igraph graph_from_adjacency_matrix set_edge_attr E V
#' @importFrom utils combn head
#' @export
build_graph <- function(data,
                        method     = "cramers_v",
                        corrected  = FALSE,
                        correct    = FALSE,
                        simulate_p = FALSE,
                        B          = 2000L,
                        alpha      = 0.5) {
  # --- resolve method ------------------------------------------------------
  # The old `corrected` flag is kept for backward compatibility: if the
  # caller passes corrected = TRUE but leaves method at its default, we
  # honour the old behaviour transparently.
  valid_methods <- c("cramers_v", "cramers_v_corrected", "nmi", "ami",
                     "bayesian_cramers_v")
  if (!method %in% valid_methods) {
    stop(
      "`method` must be one of: ",
      paste(valid_methods, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (isTRUE(corrected) && method == "cramers_v") {
    method <- "cramers_v_corrected"
  }
  use_corrected <- method == "cramers_v_corrected"
  use_nmi       <- method %in% c("nmi", "ami")
  use_adjusted  <- method == "ami"
  use_bayesian  <- method == "bayesian_cramers_v"
  
  if (use_bayesian && (!is.numeric(alpha) || length(alpha) != 1L ||
                       is.na(alpha) || alpha <= 0)) {
    stop("`alpha` must be a single positive number.", call. = FALSE)
  }
  
  non_cat <- vapply(data, function(col) {
    !is.factor(col) && !is.character(col) && !is.logical(col)
  }, logical(1))
  
  if (any(non_cat)) {
    message(
      "Coercing the following columns to character: ",
      paste(names(data)[non_cat], collapse = ", ")
    )
    data[non_cat] <- lapply(data[non_cat], as.character)
  }
  
  n_unique <- vapply(data, function(col) length(unique(col[!is.na(col)])),
                     integer(1))
  constant <- n_unique < 2L
  
  if (any(constant)) {
    warning(
      paste0(
        "The following columns have fewer than 2 unique values and will be dropped: ",
        paste(names(data)[constant], collapse = ", ")
      ),
      call. = FALSE
    )
    data <- data[, !constant, drop = FALSE]
  }
  
  if (ncol(data) < 2L) {
    stop("After removing constant columns, fewer than 2 columns remain.", call. = FALSE)
  }
  
  vars <- names(data)
  
  if (is.null(vars) || anyNA(vars) || any(vars == "")) {
    stop("`data` must have non-empty column names.", call. = FALSE)
  }
  
  if (anyDuplicated(vars)) {
    stop("`data` must have unique column names.", call. = FALSE)
  }
  
  p <- length(vars)
  
  pairs   <- combn(p, 2L, simplify = FALSE)
  n_pairs <- length(pairs)
  
  from_v    <- character(n_pairs)
  to_v      <- character(n_pairs)
  weights   <- numeric(n_pairs)
  metrics   <- character(n_pairs)
  corr_v    <- logical(n_pairs)
  pvals     <- numeric(n_pairs)
  stats     <- numeric(n_pairs)
  dfs       <- numeric(n_pairs)
  ns        <- integer(n_pairs)
  types     <- character(n_pairs)
  estimable <- logical(n_pairs)
  
  for (i in seq_along(pairs)) {
    idx_a <- pairs[[i]][1L]
    idx_b <- pairs[[i]][2L]
    
    from_v[i] <- vars[idx_a]
    to_v[i]   <- vars[idx_b]
    
    if (use_nmi) {
      es <- nmi_assoc(
        data[[idx_a]], data[[idx_b]],
        adjusted = use_adjusted,
        x_name   = vars[idx_a],
        y_name   = vars[idx_b]
      )
    } else if (use_bayesian) {
      es <- bayesian_cramers_v(
        data[[idx_a]], data[[idx_b]],
        alpha  = alpha,
        x_name = vars[idx_a],
        y_name = vars[idx_b]
      )
    } else {
      es <- effect_size(
        data[[idx_a]], data[[idx_b]],
        corrected  = use_corrected,
        correct    = correct,
        simulate_p = simulate_p,
        B          = B,
        x_name     = vars[idx_a],
        y_name     = vars[idx_b]
      )
    }
    
    estimable[i] <- !is.na(es$effect_size)
    weights[i]   <- if (estimable[i]) es$effect_size else 0
    metrics[i]   <- if (is.na(es$metric)) NA_character_ else es$metric
    corr_v[i]    <- use_corrected
    pvals[i]     <- if (is.na(es$p_value)) NA_real_ else es$p_value
    stats[i]     <- if (is.na(es$statistic)) NA_real_ else es$statistic
    dfs[i]       <- if (is.na(es$df)) NA_real_ else es$df
    ns[i]        <- if (is.na(es$n)) NA_integer_ else es$n
    types[i]     <- if (is.na(es$type)) NA_character_ else es$type
  }
  
  pair_results <- data.frame(
    var1 = from_v,
    var2 = to_v,
    weight = weights,
    metric = metrics,
    corrected = corr_v,
    p_value = pvals,
    statistic = stats,
    df = dfs,
    n = ns,
    type = types,
    estimable = estimable,
    stringsAsFactors = FALSE
  )
  
  adj <- matrix(0, nrow = p, ncol = p, dimnames = list(vars, vars))
  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1L]
    b <- pairs[[i]][2L]
    adj[a, b] <- weights[i]
    adj[b, a] <- weights[i]
  }
  
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode     = "undirected",
    weighted = TRUE,
    diag     = FALSE
  )
  
  edge_list <- igraph::as_edgelist(g)
  
  if (nrow(edge_list) > 0L) {
    edge_order <- match(
      paste(edge_list[, 1], edge_list[, 2], sep = "\r"),
      paste(from_v, to_v, sep = "\r")
    )
    
    na_idx <- is.na(edge_order)
    if (any(na_idx)) {
      edge_order[na_idx] <- match(
        paste(edge_list[na_idx, 1], edge_list[na_idx, 2], sep = "\r"),
        paste(to_v, from_v, sep = "\r")
      )
    }
    
    if (any(is.na(edge_order))) {
      stop(
        "Failed to align computed pairwise results with graph edges.",
        call. = FALSE
      )
    }
    
    igraph::E(g)$metric    <- metrics[edge_order]
    igraph::E(g)$corrected <- corr_v[edge_order]
    igraph::E(g)$p_value   <- pvals[edge_order]
    igraph::E(g)$statistic <- stats[edge_order]
    igraph::E(g)$df        <- dfs[edge_order]
    igraph::E(g)$n         <- ns[edge_order]
    igraph::E(g)$type      <- types[edge_order]
    igraph::E(g)$estimable <- estimable[edge_order]
  }
  
  attr(g, "processed_data") <- data
  attr(g, "pair_results") <- pair_results
  g
}