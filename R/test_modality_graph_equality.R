#' Permutation test for equality of modality-graph structure
#'
#' Performs an omnibus permutation test of whether two samples show the same
#' modality-level association structure.
#'
#' @importFrom stats sd
#'
#' @param x,y Two \code{catmodgraph} objects constructed from datasets
#'   with the same variable set (same column names). The underlying row
#'   data must be stored in each object's \code{$data} component, which
#'   is standard output from \code{\link{build_modality_graph}}.
#' @param n_perm Integer. Number of permutations. Default \code{500}.
#'   For publication-grade p-values use \code{2000} or more.
#' @param statistic Character. Test statistic, one of \code{"frobenius"}
#'   (default), \code{"jaccard"}, or \code{"max"}. See Details.
#' @param test_type Character. How to compute the test statistic. One of
#'   \code{"unfiltered"} (default) or \code{"pipeline"}. See Details.
#' @param strata Optional vector of length equal to the combined sample
#'   size. If supplied, permutations are conducted \emph{within} levels
#'   of \code{strata}, preserving the joint distribution of the
#'   stratification variable under the null. Use to remove the
#'   confounding of study membership with a known nuisance variable
#'   (e.g., fieldwork year, geographic region).
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical. If \code{TRUE} (default), prints progress.
#'
#' @return An object of class \code{catmodtest} with components:
#'   \describe{
#'     \item{\code{statistic}}{Character name of the test statistic used.}
#'     \item{\code{observed}}{Numeric, the observed value of the test
#'       statistic on the two input graphs.}
#'     \item{\code{null_distribution}}{Numeric vector of length
#'       \code{n_perm} with the test statistic computed under each
#'       label permutation.}
#'     \item{\code{p_value}}{Numeric empirical p-value,
#'       \code{(sum(null >= obs) + 1) / (n_perm + 1)}.}
#'     \item{\code{n_perm}}{Integer number of permutations.}
#'     \item{\code{n_x}, \code{n_y}}{Sample sizes of the two inputs.}
#'     \item{\code{test_type}}{Character, which pipeline mode was used.}
#'     \item{\code{strata_used}}{Logical indicator of whether
#'       stratified permutation was applied.}
#'   }
#'
#' @details
#' The test evaluates whether the observed difference between two
#' modality-level association matrices is larger than expected under random
#' reassignment of sample labels. Rejection supports a difference in the
#' overall marginal association structure. The test is omnibus: it does not
#' identify which specific modality pairs drive the difference. Use
#' \code{\link{test_modality_edge_differences}} for edge-wise follow-up.
#'
#' \strong{Test statistics.} Let \eqn{A} and \eqn{B} denote the weighted
#' adjacency matrices of the two graphs on a common node set.
#' \itemize{
#'   \item \code{"frobenius"} (default) uses
#'     \eqn{\|A - B\|_F^2 = \sum_{i,j} (A_{ij} - B_{ij})^2},
#'     sensitive to all edge-weight changes and dominated by
#'     high-weight edges.
#'   \item \code{"jaccard"} uses \eqn{1 - |E_A \cap E_B|/|E_A \cup E_B|},
#'     the complement of edge-set agreement. Topological, ignores
#'     weight magnitudes.
#'   \item \code{"max"} uses \eqn{\max_{i,j} |A_{ij} - B_{ij}|},
#'     sensitive to any single sharp edge-weight change.
#' }
#'
#' \strong{Pipeline modes.}
#' \itemize{
#'   \item \code{"unfiltered"} (default): the test statistic is
#'     computed on the \emph{full unpruned} phi matrix from each
#'     sample (i.e., every pair of modalities contributes its raw
#'     phi). This tests the joint distribution cleanly, with no
#'     confounding between "edges differ" and "different edges
#'     survived pruning."
#'   \item \code{"pipeline"}: each permutation re-runs
#'     \code{\link{build_modality_graph}} plus
#'     \code{\link{prune_modality_edges}} and compares the pruned
#'     adjacency matrices. Matches the actual analysis a user ran
#'     but the resulting null mixes edge-weight and edge-set
#'     changes. Slower.
#' }
#'
#' \strong{Stratification.} When \code{strata} is supplied,
#' permutations rearrange study labels \emph{within} strata, so the
#' joint distribution of the stratification variable is preserved
#' under the null.
#' This evaluates sample-label differences conditional on the supplied
#' strata and should be interpreted as a stratified permutation analysis,
#' not as causal adjustment.
#' \strong{Assumptions.} Respondents are assumed i.i.d. within each
#' input. For clustered or repeated-measures data the test is
#' anticonservative.
#'
#' @references
#' Anderson, M. J. (2001). A new method for non-parametric multivariate
#'   analysis of variance. \emph{Austral Ecology}, 26(1), 32-46.
#'
#' van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J.,
#'   Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2022).
#'   Comparing network structures on three aspects: A permutation test.
#'   \emph{Psychological Methods}. \doi{10.1037/met0000476}
#'
#' @examples
#' # Split survey_health by sex and test whether the joint structure
#' # differs. Using a small n_perm for the example; in practice use 2000+.
#' data(survey_health)
#' df_f <- subset(survey_health, sex == "female")[, -1]
#' df_m <- subset(survey_health, sex == "male")[, -1]
#'
#' mg_f <- build_modality_graph(df_f)
#' mg_m <- build_modality_graph(df_m)
#'
#' test_result <- test_modality_graph_equality(
#'   mg_f, mg_m, n_perm = 200, seed = 1, verbose = FALSE
#' )
#' print(test_result)
#'
#' @seealso \code{\link{test_modality_edge_differences}} for post-hoc
#'   edge-wise testing; \code{\link{compare_modality_graphs}} for
#'   visual comparison.
#' @importFrom stats quantile
#' @importFrom igraph as_adjacency_matrix
#' @export
test_modality_graph_equality <- function(x, y,
                                         n_perm    = 500L,
                                         statistic = c("frobenius",
                                                       "jaccard", "max"),
                                         test_type = c("unfiltered",
                                                       "pipeline"),
                                         min_weight = 0.10,
                                         max_p      = 0.05,
                                         strata    = NULL,
                                         seed      = NULL,
                                         verbose   = TRUE) {
  
  # ---- Input validation ----
  if (!inherits(x, "catmodgraph"))
    stop("`x` must be a catmodgraph object.", call. = FALSE)
  if (!inherits(y, "catmodgraph"))
    stop("`y` must be a catmodgraph object.", call. = FALSE)
  if (is.null(x$data) || is.null(y$data))
    stop("Both `x` and `y` must carry their source data in `$data` ",
         "(default output of build_modality_graph()).", call. = FALSE)
  
  vars_x <- sort(colnames(x$data))
  vars_y <- sort(colnames(y$data))
  if (!identical(vars_x, vars_y))
    stop("`x` and `y` must be built from datasets with the same ",
         "variables (columns).", call. = FALSE)
  
  statistic <- match.arg(statistic)
  test_type <- match.arg(test_type)

  #' @param min_weight Numeric. Minimum edge-weight threshold used only when
  #'   \code{test_type = "pipeline"}. Default \code{0.10}.
  #' @param max_p Numeric. Maximum edge p-value used only when
  #'   \code{test_type = "pipeline"}. Default \code{0.05}.
  
  if (!is.numeric(n_perm) || length(n_perm) != 1L || n_perm < 10)
    stop("`n_perm` must be a single integer >= 10.", call. = FALSE)
  n_perm <- as.integer(n_perm)
  
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  
  # ---- Set up combined data ----
  df_x <- x$data[, vars_x, drop = FALSE]
  df_y <- y$data[, vars_x, drop = FALSE]
  combined <- rbind(df_x, df_y)
  labels   <- c(rep("x", nrow(df_x)), rep("y", nrow(df_y)))
  
  if (!is.null(strata)) {
    if (length(strata) != length(labels))
      stop("`strata` must have length ", length(labels),
           " (combined sample size).", call. = FALSE)
    strata <- as.factor(strata)
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---- Compute observed statistic ----
  if (verbose) cat("Computing observed test statistic...\n")
  
  observed <- .catmod_compute_stat(
    df_a = df_x, df_b = df_y,
    statistic = statistic, test_type = test_type,
    prune_args = list(min_weight = min_weight, max_p = max_p)
  )
  
  # ---- Permutation loop ----
  if (verbose)
    cat("Running", n_perm, "permutations (",
        if (is.null(strata)) "unstratified" else "stratified",
        ")...\n")
  
  null_vals <- numeric(n_perm)
  pb <- if (verbose) utils::txtProgressBar(0, n_perm, style = 3) else NULL
  
  for (b in seq_len(n_perm)) {
    perm_labels <- if (is.null(strata)) {
      sample(labels)
    } else {
      .catmod_strata_permute(labels, strata)
    }
    
    df_a_perm <- combined[perm_labels == "x", , drop = FALSE]
    df_b_perm <- combined[perm_labels == "y", , drop = FALSE]
    
    null_vals[b] <- tryCatch(
      .catmod_compute_stat(
        df_a = df_a_perm, df_b = df_b_perm,
        statistic = statistic, test_type = test_type,
        prune_args = list(min_weight = 0.10, max_p = 0.05)
      ),
      error = function(e) NA_real_
    )
    
    if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
  }
  if (!is.null(pb)) close(pb)
  
  # ---- Tally p-value ----
  n_valid <- sum(!is.na(null_vals))
  if (n_valid < 10)
    warning("Fewer than 10 valid permutations; p-value unreliable.",
            call. = FALSE)
  
  p_value <- (sum(null_vals >= observed, na.rm = TRUE) + 1) / (n_valid + 1)
  
  out <- structure(
    list(
      statistic         = statistic,
      observed          = observed,
      null_distribution = null_vals,
      p_value           = p_value,
      n_perm            = n_perm,
      n_valid           = n_valid,
      n_x               = nrow(df_x),
      n_y               = nrow(df_y),
      test_type         = test_type,
      strata_used       = !is.null(strata)
    ),
    class = "catmodtest"
  )
  
  out
}


# ----------------------------------------------------------- helpers ----

#' Internal: compute the test statistic from two data frames
#'
#' @keywords internal
#' @noRd
.catmod_compute_stat <- function(df_a, df_b, statistic, test_type,
                                 prune_args) {
  if (test_type == "unfiltered") {
    mat_a <- .catmod_phi_matrix(df_a)
    mat_b <- .catmod_phi_matrix(df_b)
  } else {
    g_a <- build_modality_graph(df_a)
    g_b <- build_modality_graph(df_b)
    g_a <- prune_modality_edges(g_a,
                                min_weight = prune_args$min_weight,
                                max_p      = prune_args$max_p)
    g_b <- prune_modality_edges(g_b,
                                min_weight = prune_args$min_weight,
                                max_p      = prune_args$max_p)
    mat_a <- igraph::as_adjacency_matrix(g_a$graph, attr = "weight",
                                         sparse = FALSE)
    mat_b <- igraph::as_adjacency_matrix(g_b$graph, attr = "weight",
                                         sparse = FALSE)
  }
  
  # Align on common modality set
  common <- intersect(rownames(mat_a), rownames(mat_b))
  if (length(common) < 2) return(NA_real_)
  
  mat_a <- mat_a[common, common, drop = FALSE]
  mat_b <- mat_b[common, common, drop = FALSE]
  mat_a[is.na(mat_a)] <- 0
  mat_b[is.na(mat_b)] <- 0
  
  diff <- mat_a - mat_b
  
  switch(statistic,
         frobenius = sum(diff^2),
         jaccard   = {
           e_a <- which(mat_a > 0, arr.ind = FALSE)
           e_b <- which(mat_b > 0, arr.ind = FALSE)
           u <- length(union(e_a, e_b))
           if (u == 0) 0 else 1 - length(intersect(e_a, e_b)) / u
         },
         max       = max(abs(diff))
  )
}

#' Internal: compute the full phi matrix (unpruned)
#'
#' @keywords internal
#' @noRd
.catmod_phi_matrix <- function(df) {
  # Build modality graph without pruning; extract weighted adjacency
  g <- suppressWarnings(build_modality_graph(df))
  igraph::as_adjacency_matrix(g$graph, attr = "weight", sparse = FALSE)
}

#' Internal: stratified permutation
#'
#' @keywords internal
#' @noRd
.catmod_strata_permute <- function(labels, strata) {
  new_labels <- labels
  for (lv in levels(strata)) {
    idx <- which(strata == lv)
    new_labels[idx] <- sample(labels[idx])
  }
  new_labels
}


# ----------------------------------------------------------- methods ----

#' @export
print.catmodtest <- function(x, ...) {
  cat("Permutation test of modality-graph equality\n")
  cat("  Statistic       :", x$statistic, "\n")
  cat("  Test type       :", x$test_type, "\n")
  cat("  Stratified      :", x$strata_used, "\n")
  cat("  Sample sizes    : n_x =", x$n_x, ", n_y =", x$n_y, "\n")
  cat("  Permutations    :", x$n_valid, "valid /", x$n_perm, "\n")
  cat("  Observed stat   :", signif(x$observed, 4), "\n")
  cat("  Null mean / sd  :",
      signif(mean(x$null_distribution, na.rm = TRUE), 4), "/",
      signif(sd(x$null_distribution, na.rm = TRUE), 4), "\n")
  cat("  Empirical p     :", signif(x$p_value, 4), "\n\n")
  
  decision <- if (x$p_value < 0.05) {
    "REJECT null at alpha = 0.05 (graphs differ significantly)"
  } else {
    "FAIL TO REJECT null at alpha = 0.05"
  }
  cat(" ", decision, "\n")
  
  invisible(x)
}

#' @export
summary.catmodtest <- function(object, ...) {
  q <- stats::quantile(object$null_distribution,
                       probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  list(
    statistic   = object$statistic,
    test_type   = object$test_type,
    observed    = object$observed,
    null_q025   = unname(q[1]),
    null_median = unname(q[2]),
    null_q975   = unname(q[3]),
    p_value     = object$p_value,
    decision    = if (object$p_value < 0.05) "reject" else "fail_to_reject"
  )
}

#' @export
plot.catmodtest <- function(x, bins = 40, ...) {
  graphics::hist(
    x$null_distribution,
    breaks = bins,
    col    = "#AFA9EC",
    border = "white",
    main   = sprintf("Permutation null: %s (p = %.3f)",
                     x$statistic, x$p_value),
    xlab   = sprintf("%s statistic under H_0", x$statistic),
    ...
  )
  graphics::abline(v = x$observed, col = "#D85A30", lwd = 3)
  graphics::legend("topright",
                   legend = c("null", "observed"),
                   fill   = c("#AFA9EC", NA),
                   border = c("white", NA),
                   lty    = c(NA, 1),
                   lwd    = c(NA, 3),
                   col    = c(NA, "#D85A30"),
                   bty    = "n")
  invisible(x)
}
