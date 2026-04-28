#' Edge-wise post-hoc test for modality-network differences
#'
#' After a significant global test with
#' \code{\link{test_modality_graph_equality}}, this function identifies
#' which specific edges contribute most strongly to the observed difference. For each edge, an
#' empirical two-sided p-value is computed from label-permutation
#' distributions of the absolute edge-weight difference, and
#' Benjamini-Hochberg FDR correction is applied across all tested
#' edges.
#'
#' @param x,y Two \code{catmodgraph} objects on the same variable set,
#'   each carrying their \code{$data} component.
#' @param n_perm Integer, number of permutations. Default \code{500}.
#' @param edges Character, which edges to test. One of \code{"all"}
#'   (default, every upper-triangle pair of modalities) or
#'   \code{"union"} (only edges present in at least one of the two
#'   input graphs). Using \code{"union"} reduces the multiple-testing
#'   burden.
#' @param p_adjust Character, adjustment method for multiple testing,
#'   passed to \code{\link[stats]{p.adjust}}. Default \code{"BH"}.
#' @param strata Optional stratification vector of length equal to the
#'   combined sample size. See \code{\link{test_modality_graph_equality}}.
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Logical. If \code{TRUE} (default), prints progress.
#'
#' @return An object of class \code{catmodedgetest} with components:
#'   \describe{
#'     \item{\code{edge_table}}{Data frame with one row per tested edge:
#'       columns \code{from}, \code{to}, \code{weight_x}, \code{weight_y},
#'       \code{obs_diff}, \code{p_empirical}, \code{p_adjusted}.}
#'     \item{\code{n_perm}}{Integer number of permutations.}
#'     \item{\code{n_x}, \code{n_y}}{Sample sizes.}
#'     \item{\code{edges}}{Character, which edge-subset criterion was used.}
#'     \item{\code{p_adjust_method}}{Adjustment method.}
#'     \item{\code{strata_used}}{Logical.}
#'   }
#'
#' @details
#' For each edge, the permutation null is that the observed edge-weight
#' difference is no larger than expected under random reassignment of
#' sample labels. The test statistic is the observed edge-weight
#' difference; the null distribution is obtained by recomputing that
#' difference across \code{n_perm} random label reassignments of the
#' combined data. Edges not present in the pruned graphs are treated
#' as having weight 0.
#'
#' The function is designed to be called \emph{after} a significant
#' global test. Calling it without a significant omnibus inflates the
#' family-wise Type I error beyond the nominal FDR level, because the
#' omnibus test acts as a gatekeeper under the closed-testing
#' principle.
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
#' head(edge_test$edge_table)
#' }
#'
#' @seealso \code{\link{test_modality_graph_equality}} for the global test.
#' @importFrom stats p.adjust
#' @importFrom igraph as_adjacency_matrix
#' @export
test_modality_edge_differences <- function(x, y,
                                           n_perm   = 500L,
                                           edges    = c("all", "union"),
                                           p_adjust = "BH",
                                           strata   = NULL,
                                           seed     = NULL,
                                           verbose  = TRUE) {
  
  # ---- Input validation ----
  if (!inherits(x, "catmodgraph"))
    stop("`x` must be a catmodgraph object.", call. = FALSE)
  if (!inherits(y, "catmodgraph"))
    stop("`y` must be a catmodgraph object.", call. = FALSE)
  if (is.null(x$data) || is.null(y$data))
    stop("Both `x` and `y` must carry their source data in `$data`.",
         call. = FALSE)
  
  vars_x <- sort(colnames(x$data))
  vars_y <- sort(colnames(y$data))
  if (!identical(vars_x, vars_y))
    stop("`x` and `y` must be built from datasets with the same variables.",
         call. = FALSE)
  
  edges <- match.arg(edges)
  
  if (!is.numeric(n_perm) || length(n_perm) != 1L || n_perm < 10)
    stop("`n_perm` must be a single integer >= 10.", call. = FALSE)
  n_perm <- as.integer(n_perm)
  
  if (!is.character(p_adjust) || length(p_adjust) != 1L)
    stop("`p_adjust` must be a single string.", call. = FALSE)
  
  # ---- Observed phi matrices ----
  if (verbose) cat("Computing observed edge differences...\n")
  
  mat_x <- .catmod_phi_matrix(x$data[, vars_x, drop = FALSE])
  mat_y <- .catmod_phi_matrix(y$data[, vars_x, drop = FALSE])
  
  common <- intersect(rownames(mat_x), rownames(mat_y))
  if (length(common) < 2)
    stop("Fewer than 2 modalities are common to both graphs.",
         call. = FALSE)
  mat_x <- mat_x[common, common, drop = FALSE]
  mat_y <- mat_y[common, common, drop = FALSE]
  mat_x[is.na(mat_x)] <- 0
  mat_y[is.na(mat_y)] <- 0
  
  obs_diff <- mat_x - mat_y
  
  # ---- Build edge table (upper-triangle only) ----
  upper_idx <- which(upper.tri(obs_diff), arr.ind = TRUE)
  edge_table <- data.frame(
    from      = rownames(obs_diff)[upper_idx[, 1]],
    to        = colnames(obs_diff)[upper_idx[, 2]],
    weight_x  = mat_x[upper_idx],
    weight_y  = mat_y[upper_idx],
    obs_diff  = obs_diff[upper_idx],
    stringsAsFactors = FALSE
  )
  
  # Restrict to union of edges if requested
  if (edges == "union") {
    keep <- (edge_table$weight_x != 0) | (edge_table$weight_y != 0)
    edge_table <- edge_table[keep, , drop = FALSE]
    if (nrow(edge_table) == 0L)
      stop("No edges present in either graph; nothing to test.",
           call. = FALSE)
  }
  
  n_edges <- nrow(edge_table)
  if (verbose) cat("Testing", n_edges, "edges...\n")
  
  # ---- Combined data and permutation ----
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
  
  perm_diffs <- matrix(0, nrow = n_edges, ncol = n_perm)
  edge_keys  <- paste(edge_table$from, edge_table$to, sep = "___")
  
  pb <- if (verbose) utils::txtProgressBar(0, n_perm, style = 3) else NULL
  
  for (b in seq_len(n_perm)) {
    perm_labels <- if (is.null(strata)) {
      sample(labels)
    } else {
      .catmod_strata_permute(labels, strata)
    }
    
    df_a_perm <- combined[perm_labels == "x", , drop = FALSE]
    df_b_perm <- combined[perm_labels == "y", , drop = FALSE]
    
    m_a <- tryCatch(.catmod_phi_matrix(df_a_perm), error = function(e) NULL)
    m_b <- tryCatch(.catmod_phi_matrix(df_b_perm), error = function(e) NULL)
    
    if (is.null(m_a) || is.null(m_b)) {
      # leave column at 0
      if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
      next
    }
    
    cp <- intersect(rownames(m_a), rownames(m_b))
    if (length(cp) < 2L) {
      if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
      next
    }
    m_a <- m_a[cp, cp, drop = FALSE]; m_a[is.na(m_a)] <- 0
    m_b <- m_b[cp, cp, drop = FALSE]; m_b[is.na(m_b)] <- 0
    perm_diff_mat <- m_a - m_b
    
    # Pull out the observed edges that exist in this permutation
    ui <- which(upper.tri(perm_diff_mat), arr.ind = TRUE)
    perm_keys <- paste(rownames(perm_diff_mat)[ui[, 1]],
                       colnames(perm_diff_mat)[ui[, 2]], sep = "___")
    perm_vals <- perm_diff_mat[ui]
    
    # Lookup
    ix <- match(edge_keys, perm_keys)
    valid <- !is.na(ix)
    perm_diffs[valid, b] <- perm_vals[ix[valid]]
    # remaining rows stay at 0 (edge didn't exist in this permutation)
    
    if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
  }
  if (!is.null(pb)) close(pb)
  
  # ---- Empirical two-sided p-values + BH adjustment ----
  p_emp <- vapply(seq_len(n_edges), function(i) {
    (sum(abs(perm_diffs[i, ]) >= abs(edge_table$obs_diff[i])) + 1) /
      (n_perm + 1)
  }, numeric(1L))
  
  p_adj <- stats::p.adjust(p_emp, method = p_adjust)
  
  edge_table$p_empirical <- p_emp
  edge_table$p_adjusted  <- p_adj
  edge_table$abs_diff    <- abs(edge_table$obs_diff)
  
  edge_table <- edge_table[order(edge_table$p_adjusted,
                                 -edge_table$abs_diff), , drop = FALSE]
  
  out <- structure(
    list(
      edge_table      = edge_table,
      n_perm          = n_perm,
      n_x             = nrow(df_x),
      n_y             = nrow(df_y),
      edges           = edges,
      p_adjust_method = p_adjust,
      strata_used     = !is.null(strata)
    ),
    class = "catmodedgetest"
  )
  
  out
}


# ----------------------------------------------------------- methods ----

#' @export
print.catmodedgetest <- function(x, n = 10, ...) {
  cat("Edge-wise permutation test of modality-network differences\n")
  cat("  Edges tested    :", nrow(x$edge_table), "\n")
  cat("  Edge selection  :", x$edges, "\n")
  cat("  Sample sizes    : n_x =", x$n_x, ", n_y =", x$n_y, "\n")
  cat("  Permutations    :", x$n_perm, "\n")
  cat("  p-adjust method :", x$p_adjust_method, "\n")
  cat("  Stratified      :", x$strata_used, "\n")
  cat("  Edges sig. at 5%%:",
      sum(x$edge_table$p_adjusted < 0.05), "\n\n")
  
  cat("  Top", n, "edges by adjusted p-value:\n")
  print(utils::head(x$edge_table, n), digits = 3, row.names = FALSE)
  invisible(x)
}

#' @export
summary.catmodedgetest <- function(object, ...) {
  n_sig <- sum(object$edge_table$p_adjusted < 0.05)
  list(
    n_edges_tested = nrow(object$edge_table),
    n_significant  = n_sig,
    fdr_level      = 0.05,
    strongest      = utils::head(object$edge_table, 10)
  )
}

#' @export
plot.catmodedgetest <- function(x, n = 20, ...) {
  et <- utils::head(x$edge_table, n)
  et$edge_label <- paste(et$from, "\u2194", et$to)
  et <- et[order(et$abs_diff), , drop = FALSE]
  
  graphics::par(mar = c(5, 14, 3, 2))
  on.exit(graphics::par(mar = c(5, 4, 4, 2)), add = TRUE)
  
  cols <- ifelse(et$p_adjusted < 0.05, "#534AB7", "#D3D1C7")
  graphics::barplot(
    et$obs_diff,
    horiz      = TRUE,
    names.arg  = et$edge_label,
    las        = 1,
    col        = cols,
    border     = NA,
    cex.names  = 0.75,
    xlab       = "Observed edge-weight difference (x - y)",
    main       = sprintf("Top %d edges by adjusted p-value", nrow(et)),
    ...
  )
  graphics::abline(v = 0, col = "black", lwd = 1)
  graphics::legend("topright",
                   legend = c("significant (adj. p < 0.05)",
                              "not significant"),
                   fill   = c("#534AB7", "#D3D1C7"),
                   border = NA, bty = "n", cex = 0.8)
  invisible(x)
}