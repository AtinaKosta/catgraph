#' Prune edges from a catgraph by effect size or adjusted p-value
#'
#' Removes edges whose effect size or (adjusted) p-value does not meet a
#' specified threshold, returning a new \code{catgraph} object with a
#' sparser graph. Multiple-testing adjustment is applied across all edges
#' by default.
#'
#' @param x A \code{catgraph} object.
#' @param min_weight Numeric, non-negative. Edges with effect size strictly
#'   below this value are removed. Since phi and Cramer's V lie in [0, 1],
#'   meaningful thresholds are in that range; values >= 1 remove all edges.
#'   Effect-size pruning is the \strong{primary} filter and is always a
#'   safer choice than p-value pruning, because the package is centred on
#'   effect sizes. Default \code{0} (no filtering).
#' @param max_p Numeric in [0, 1]. Edges with \emph{adjusted} p-value
#'   strictly above this value are removed. Adjustment method is controlled
#'   by \code{p_adjust}. Default \code{1} (no filtering).
#' @param p_adjust Character. Multiple-testing correction applied across all
#'   edges (\code{choose(p, 2)} simultaneous tests). One of:
#'   \describe{
#'     \item{\code{"BH"}}{Benjamini-Hochberg false discovery rate (default).
#'       Recommended for exploratory work.}
#'     \item{\code{"holm"}}{Holm-Bonferroni step-down; strong family-wise
#'       error rate control.}
#'     \item{\code{"bonferroni"}}{Bonferroni; conservative FWER control.}
#'     \item{\code{"none"}}{Raw p-values (unadjusted). Not recommended when
#'       many variables are analysed; retained for reproducing pre-0.4.0
#'       behaviour.}
#'   }
#' @param remove_isolates Logical. If \code{TRUE}, vertices with degree 0
#'   after pruning are also removed. Default \code{FALSE}.
#'
#' @return A \code{catgraph} object with the filtered graph. The graph gains
#'   two new edge attributes: \code{p_value_adj} (adjusted p-values) and
#'   \code{p_adjust_method} (the method string).
#'
#' @details
#' Pruning uses the edge attribute \code{weight} (the active effect size)
#' and the \emph{adjusted} \code{p_value}. Both thresholds apply
#' simultaneously; an edge is retained only when \strong{both} conditions
#' are met.
#'
#' Edges with \code{NA} weights or p-values (from degenerate pairs) are
#' always removed.
#'
#' \strong{Chained calls and multiplicity scoping.} Multiple-testing
#' correction is applied across the edges present in the graph at the
#' time of the call. When \code{prune_edges()} is called on a graph that
#' has \emph{already} been pruned with a non-\code{"none"} p-value
#' adjustment, the second call re-adjusts on the surviving subset, not
#' on the original \code{choose(p, 2)} tests. This is anti-conservative
#' and the function emits a warning in this case. To change the
#' adjustment method mid-analysis, rebuild the \code{catgraph} with
#' \code{\link{catgraph}()} and prune once; do not chain two adjusted
#' prunes. A single \code{prune_edges()} call that specifies both
#' \code{min_weight} and \code{max_p} is always safe because the BH /
#' Holm denominators are computed before any filtering.
#'
#' Conventional Cohen (1988) thresholds for phi and Cramer's V: small
#' \eqn{\approx 0.1}, medium \eqn{\approx 0.3}, large \eqn{\geq 0.5}.
#'
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#'   rate: a practical and powerful approach to multiple testing.
#'   \emph{JRSS-B}, 57(1), 289--300.
#'   \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#'   \emph{Scandinavian Journal of Statistics}, 6(2), 65--70.
#'
#' Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral
#'   Sciences} (2nd ed.). Lawrence Erlbaum Associates.
#'
#' @examples
#' df <- expand_table(Titanic)
#' cg <- catgraph(df)
#'
#' # Default: BH-adjusted p-values, effect-size floor of 0.1
#' cg_pruned <- prune_edges(cg, min_weight = 0.1, max_p = 0.05)
#' cg_pruned
#'
#' # Stricter: Holm adjustment
#' prune_edges(cg, min_weight = 0.1, max_p = 0.05, p_adjust = "holm")
#'
#' # Pre-0.4.0 behaviour (raw p-values)
#' prune_edges(cg, min_weight = 0.1, max_p = 0.05, p_adjust = "none")
#'
#' @seealso \code{\link{catgraph}}, \code{\link{detect_clusters}}
#' @importFrom igraph delete_edges delete_vertices E degree ecount vcount
#' @importFrom stats p.adjust na.omit
#' @export
prune_edges <- function(x,
                        min_weight      = 0,
                        max_p           = 1,
                        p_adjust        = c("BH", "holm", "bonferroni", "none"),
                        remove_isolates = FALSE) {
  
  if (!inherits(x, "catgraph")) {
    stop("`x` must be a catgraph object.", call. = FALSE)
  }
  
  p_adjust <- match.arg(p_adjust)
  
  if (!is.numeric(min_weight) || length(min_weight) != 1L || is.na(min_weight) ||
      min_weight < 0) {
    stop("`min_weight` must be a single non-negative number.", call. = FALSE)
  }
  
  if (!is.numeric(max_p) || length(max_p) != 1L || is.na(max_p) ||
      max_p < 0 || max_p > 1) {
    stop("`max_p` must be a single number in [0, 1].", call. = FALSE)
  }
  
  if (!is.logical(remove_isolates) || length(remove_isolates) != 1L || is.na(remove_isolates)) {
    stop("`remove_isolates` must be TRUE or FALSE.", call. = FALSE)
  }
  
  g <- x$graph
  
  if (igraph::ecount(g) == 0L) {
    if (p_adjust != "none" &&
        "p_value_adj" %in% igraph::edge_attr_names(g)) {
      warning(
        paste0(
          "prune_edges() called on a catgraph with no edges. ",
          "A prior multiple-testing adjustment is present but cannot be re-applied ",
          "(no edges to adjust). Returning the empty graph unchanged."
        ),
        call. = FALSE
      )
    }
    x$n_vars  <- igraph::vcount(g)
    x$n_pairs <- igraph::ecount(g)
    return(x)
  }
  
  w <- igraph::E(g)$weight
  p <- igraph::E(g)$p_value
  
  prior_adj <- "p_value_adj" %in% igraph::edge_attr_names(g)
  prior_method <- if (prior_adj) {
    unique(stats::na.omit(igraph::E(g)$p_adjust_method))
  } else {
    character(0)
  }
  
  if (prior_adj && p_adjust != "none") {
    warning(
      paste0(
        "This catgraph object has already been adjusted for multiple testing (method: ",
        if (length(prior_method) == 1L) prior_method else "unknown/mixed",
        "). Re-adjusting on the surviving subset uses m = ", length(p),
        " rather than the original choose(p, 2) = ", choose(x$n_vars, 2L),
        ", which is anti-conservative. If you want to change the adjustment method, ",
        "rebuild the catgraph with catgraph() and prune once. Proceeding with the ",
        "requested adjustment; this is correct only if you genuinely intend a nested-family adjustment."
      ),
      call. = FALSE
    )
  }
  
  if (p_adjust == "none") {
    p_adj <- p
  } else {
    p_adj <- stats::p.adjust(p, method = p_adjust)
  }
  
  igraph::E(g)$p_value_adj     <- p_adj
  igraph::E(g)$p_adjust_method <- p_adjust
  
  drop <- is.na(w) | is.na(p_adj) | (w < min_weight) | (p_adj > max_p)
  
  if (!is.null(x$pair_results)) {
    pr <- x$pair_results
    el <- igraph::as_edgelist(g)
    
    edge_tbl <- data.frame(
      var1 = el[, 1],
      var2 = el[, 2],
      p_value_adj = p_adj,
      p_adjust_method = p_adjust,
      retained_after_prune = !drop,
      stringsAsFactors = FALSE
    )
    
    key_pr_1 <- paste(pr$var1, pr$var2, sep = "\r")
    key_pr_2 <- paste(pr$var2, pr$var1, sep = "\r")
    key_ed   <- paste(edge_tbl$var1, edge_tbl$var2, sep = "\r")
    
    idx <- match(key_pr_1, key_ed)
    rev_idx <- is.na(idx)
    if (any(rev_idx)) idx[rev_idx] <- match(key_pr_2[rev_idx], key_ed)
    
    pr$p_value_adj <- NA_real_
    pr$p_adjust_method <- NA_character_
    pr$retained_after_prune <- FALSE
    
    matched <- !is.na(idx)
    pr$p_value_adj[matched] <- edge_tbl$p_value_adj[idx[matched]]
    pr$p_adjust_method[matched] <- edge_tbl$p_adjust_method[idx[matched]]
    pr$retained_after_prune[matched] <- edge_tbl$retained_after_prune[idx[matched]]
    
    x$pair_results <- pr
  }
  
  if (any(drop)) {
    g <- igraph::delete_edges(g, which(drop))
  }
  
  if (remove_isolates) {
    iso <- which(igraph::degree(g) == 0L)
    if (length(iso) > 0L) {
      g <- igraph::delete_vertices(g, iso)
    }
  }
  
  x$graph   <- g
  x$n_vars  <- igraph::vcount(g)
  x$n_pairs <- igraph::ecount(g)
  x
}