#' Dense pairwise similarity matrix of categorical variables
#'
#' Computes the full \eqn{p \times p} matrix of pairwise effect sizes for
#' categorical variables, including pairs with zero association.
#'
#' This function is the \strong{correct input for heatmap-style visualisation}
#' and any analysis that requires a dense similarity matrix. The
#' \pkg{igraph} object returned by \code{\link{catgraph}} is the correct
#' input for \emph{topology} (centrality, clustering, density, community
#' detection): it represents zero-association pairs as absent edges and
#' therefore would give misleading heatmaps.
#'
#' @param data A data frame of categorical variables (same requirements as
#'   \code{\link{catgraph}}).
#' @param corrected Logical. If \code{TRUE}, use the bias-corrected
#'   estimator (Bergsma, 2013). Default \code{FALSE}.
#' @param correct Logical. Yates' continuity correction for chi-square.
#'   Default \code{FALSE}.
#' @param simulate_p Logical. Monte Carlo simulation for p-values (affects
#'   only the p-value matrix, not the effect-size matrix). Default
#'   \code{FALSE}.
#' @param B Integer. Monte Carlo resamples. Default \code{2000L}.
#' @param what Character. What to return: \code{"effect_size"} (default),
#'   \code{"p_value"}, \code{"n"} (pairwise-complete observation count), or
#'   \code{"all"} (a list of matrices).
#'
#' @return A symmetric numeric matrix (or a list of matrices when
#'   \code{what = "all"}). Diagonal is \code{NA}. Row and column names are
#'   the variable names.
#'
#' @details
#' This function duplicates the computation done by
#' \code{\link{build_graph}} but does not collapse the result into a graph,
#' so all pairs are represented. In v0.3.0 and earlier, the same output was
#' extracted from the graph via \code{assoc_matrix()}, but because the
#' graph forced zero-weight pairs to \code{.Machine$double.eps}, the
#' resulting matrix silently conflated "zero association" with "near-zero
#' association". From 0.4.0 onwards, use \code{assoc_similarity()} when you
#' want the \emph{full} dense matrix and \code{assoc_matrix()} (the graph
#' extractor) when you want the matrix of \emph{actual} edges.
#'
#' @examples
#' df <- expand_table(Titanic)
#' S <- assoc_similarity(df)
#' round(S, 3)
#'
#' # All three components at once
#' out <- assoc_similarity(df, what = "all")
#' str(out, max.level = 1)
#'
#' @seealso \code{\link{build_graph}}, \code{\link{assoc_matrix}},
#'   \code{\link{plot_heatmap}}
#' @importFrom utils combn
#' @export
assoc_similarity <- function(data,
                             corrected  = FALSE,
                             correct    = FALSE,
                             simulate_p = FALSE,
                             B          = 2000L,
                             what       = c("effect_size", "p_value",
                                            "n", "all")) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.")
  }
  if (ncol(data) < 2L) {
    stop("`data` must have at least 2 columns.")
  }
  what <- match.arg(what)

  # Mirror build_graph()'s coercion and constant-column handling so the two
  # stay in sync. Silent here (build_graph() already messages).
  non_cat <- vapply(data, function(col) {
    !is.factor(col) && !is.character(col) && !is.logical(col)
  }, logical(1))
  if (any(non_cat)) {
    data[non_cat] <- lapply(data[non_cat], as.character)
  }
  n_unique <- vapply(data, function(col) length(unique(col[!is.na(col)])),
                     integer(1))
  data <- data[, n_unique >= 2L, drop = FALSE]
  if (ncol(data) < 2L) {
    stop("After removing constant columns, fewer than 2 columns remain.")
  }

  vars <- names(data)
  p    <- length(vars)
  pairs <- combn(p, 2L, simplify = FALSE)

  S  <- matrix(NA_real_, p, p, dimnames = list(vars, vars))
  P  <- matrix(NA_real_, p, p, dimnames = list(vars, vars))
  N  <- matrix(NA_real_, p, p, dimnames = list(vars, vars))

  for (pp in pairs) {
    a <- pp[1L]; b <- pp[2L]
    es <- suppressWarnings(effect_size(
      data[[a]], data[[b]],
      corrected = corrected, correct = correct,
      simulate_p = simulate_p, B = B
    ))
    S[a, b] <- S[b, a] <- if (is.na(es$effect_size)) NA_real_ else es$effect_size
    P[a, b] <- P[b, a] <- if (is.na(es$p_value))     NA_real_ else es$p_value
    N[a, b] <- N[b, a] <- if (is.na(es$n))           NA_real_ else es$n
  }

  switch(
    what,
    effect_size = S,
    p_value     = P,
    n           = N,
    all         = list(effect_size = S, p_value = P, n = N)
  )
}
