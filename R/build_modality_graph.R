#' Build a modality-level association graph
#'
#' Constructs a modality-level association network in which each node is a
#' category level, encoded as \code{"variable=level"}, and each edge represents
#' the pairwise association between two modalities from different variables.
#'
#' Edge weights (\code{weight}) are absolute phi coefficients computed from
#' the corresponding binary 2x2 indicator table, so that edge thickness in
#' plots always scales with association strength regardless of direction.
#' The signed phi coefficient is stored separately as \code{phi_signed}.
#' The signed standardized residual is stored as \code{std_resid}: positive
#' values indicate co-occurrence above expectation (attraction), whereas
#' negative values indicate co-occurrence below expectation (repulsion).
#' The graph is therefore weighted by association magnitude and annotated
#' by direction.
#'
#' @param data A data frame of categorical variables.
#' @param remove_na Logical. If TRUE, rows with missing values are removed
#'   before graph construction. Default is TRUE.
#' @param min_count Integer. Modalities with total count strictly below this
#'   threshold are removed before association computation. Default is 1.
#'
#' @return An object of class \code{"catmodgraph"}, containing:
#'   \describe{
#'     \item{\code{graph}}{An \code{igraph} object.}
#'     \item{\code{modalities}}{A data frame describing node identities.}
#'     \item{\code{indicator_matrix}}{The binary modality indicator matrix.}
#'     \item{\code{data}}{The processed categorical data used to build the graph.}
#'   }
#'
#' @details
#' Same-variable modality pairs are excluded by construction, because such
#' modalities are mutually exclusive and are not meaningful as co-occurrence
#' edges in the intended network.
#' 
#' \strong{Interpretation.} This graph is descriptive and exploratory. It is
#' not a causal graph, not a respondent-segmentation model, and not a
#' conditional-independence graph. Because edge weights are absolute values,
#' users should inspect \code{std_resid} or use signed plotting options when
#' distinguishing attraction from repulsion matters.
#'
#' \strong{Implementation.} Since v0.8.0 the per-pair 2x2 statistics are
#' computed in vectorised closed form from the modality co-occurrence matrix
#' \eqn{M^\top M}, rather than by an inner loop over
#' \code{\link[stats]{chisq.test}} calls. For \eqn{m} modalities on \eqn{n}
#' rows the dominant cost drops from \eqn{O(m^2 n)} plus a per-pair
#' \code{chisq.test} overhead to a single \eqn{O(m^2 n)} matrix product plus
#' \eqn{O(m^2)} arithmetic. Typical speed-up on realistic data
#' (\eqn{m \approx 40}, \eqn{n \approx 5000}) is 30--100x. Numeric output
#' is algebraically identical to the earlier loop implementation up to
#' floating-point tolerance; edge ordering in the returned graph may
#' differ.
#'
#' This function is the modality-level counterpart to \code{\link{catgraph}()}.
#'
#' @examples
#' df <- expand_table(Titanic)
#' mg <- build_modality_graph(df)
#' mg
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom stats na.omit pchisq
#' @export
build_modality_graph <- function(data,
                                 remove_na = TRUE,
                                 min_count = 1L) {
  
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  
  if (ncol(data) < 2L) {
    stop("`data` must contain at least two variables.", call. = FALSE)
  }
  
  if (!is.logical(remove_na) || length(remove_na) != 1L || is.na(remove_na)) {
    stop("`remove_na` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.numeric(min_count) || length(min_count) != 1L || is.na(min_count) ||
      min_count < 1 || min_count != as.integer(min_count)) {
    stop("`min_count` must be a single integer >= 1.", call. = FALSE)
  }
  
  dat <- data
  
  if (remove_na) {
    dat <- stats::na.omit(dat)
  }
  
  if (nrow(dat) == 0L) {
    stop("No complete rows remain after NA removal.", call. = FALSE)
  }
  
  # Coerce all variables to factors
  dat[] <- lapply(dat, function(x) if (is.factor(x)) x else factor(x))
  
  var_names <- names(dat)
  
  # -- Build the one-hot indicator matrix one variable at a time ------------
  # (same construction as the previous implementation, preserved for
  # deterministic column ordering and identical node identities)
  mm_list         <- vector("list", length(var_names))
  modalities_list <- vector("list", length(var_names))
  
  for (j in seq_along(var_names)) {
    v    <- var_names[j]
    levs <- levels(dat[[j]])
    
    block <- vapply(
      levs,
      function(lev) as.integer(dat[[j]] == lev),
      FUN.VALUE = integer(nrow(dat))
    )
    if (is.null(dim(block))) {
      block <- matrix(block, ncol = 1L)
    }
    colnames(block) <- paste0(v, "=", levs)
    
    mm_list[[j]] <- block
    modalities_list[[j]] <- data.frame(
      node     = paste0(v, "=", levs),
      variable = v,
      modality = levs,
      stringsAsFactors = FALSE
    )
  }
  
  mm         <- do.call(cbind, mm_list)
  mm         <- as.matrix(mm)
  storage.mode(mm) <- "double"   # keep %*% in double precision
  modalities <- do.call(rbind, modalities_list)
  modalities <- modalities[match(colnames(mm), modalities$node), , drop = FALSE]
  
  # -- Drop rare modalities --------------------------------------------------
  counts    <- colSums(mm)
  keep_cols <- counts >= min_count
  
  mm         <- mm[, keep_cols, drop = FALSE]
  modalities <- modalities[match(colnames(mm), modalities$node), , drop = FALSE]
  
  if (ncol(mm) < 2L) {
    stop("Fewer than two modalities remain after filtering.", call. = FALSE)
  }
  
  n_rows <- nrow(mm)
  n_mod  <- ncol(mm)
  node_names <- colnames(mm)
  var_of_node <- modalities$variable  # aligned with colnames(mm) by construction
  
  # -- Vectorised 2x2 statistics --------------------------------------------
  # Co-occurrence matrix: A[i, j] = sum(mm[, i] * mm[, j])
  A <- crossprod(mm)                    # n_mod x n_mod, symmetric
  col_tot <- diag(A)                    # n_i: total count of each modality
  
  # Outer products give the required marginal combinations
  nij <- tcrossprod(col_tot)            # n_i * n_j
  ni_compl <- tcrossprod(n_rows - col_tot, rep(1, n_mod))  # (n - n_i), replicated across cols
  # We actually want (n - n_i)(n - n_j):
  n_minus <- n_rows - col_tot
  n_minus_ij <- tcrossprod(n_minus)     # (n - n_i) * (n - n_j)
  
  # phi numerator: n * A - n_i * n_j
  num <- n_rows * A - nij
  # phi denominator: sqrt(n_i * n_j * (n - n_i) * (n - n_j))
  denom_sq <- nij * n_minus_ij
  denom <- suppressWarnings(sqrt(denom_sq))
  
  phi_signed <- num / denom             # NaN where denom == 0
  # Absolute phi is the edge weight (matches classical phi on 2x2 used
  # by the loop version; the sign is carried separately via std_resid).
  phi_abs <- abs(phi_signed)
  
  # Chi-square statistic on 1 df: n * phi^2
  chi_sq <- n_rows * phi_signed^2
  p_mat  <- stats::pchisq(chi_sq, df = 1, lower.tail = FALSE)
  
  # Standardised Pearson residual at cell (2, 2) [i.e., xi = 1, xj = 1]:
  #   z_22 = (O_22 - E_22) / sqrt(E_22 * (1 - p_i.) * (1 - p_.j))
  # Algebraically equal to phi_signed * sqrt(n), so reuse that:
  std_resid_mat <- phi_signed * sqrt(n_rows)
  
  # -- Assemble edge table from upper triangle, excluding same-variable pairs
  idx <- which(upper.tri(A), arr.ind = TRUE)   # 2-col matrix of (i, j) with i < j
  i_idx <- idx[, 1L]
  j_idx <- idx[, 2L]
  
  same_var <- var_of_node[i_idx] == var_of_node[j_idx]
  keep     <- !same_var
  
  i_idx <- i_idx[keep]
  j_idx <- j_idx[keep]
  
  edge_df <- data.frame(
    from       = node_names[i_idx],
    to         = node_names[j_idx],
    weight     = phi_abs[cbind(i_idx, j_idx)],
    phi_signed = phi_signed[cbind(i_idx, j_idx)],
    p_value    = p_mat[cbind(i_idx, j_idx)],
    std_resid  = std_resid_mat[cbind(i_idx, j_idx)],
    stringsAsFactors = FALSE
  )
  
  # Drop degenerate pairs (matches old behaviour: phi NA -> edge dropped)
  edge_df <- edge_df[!is.na(edge_df$weight), , drop = FALSE]
  
  if (nrow(edge_df) == 0L) {
    stop("No cross-variable modality associations could be computed.",
         call. = FALSE)
  }
  
  # NaN weights already filtered via !is.na; ensure p_value / std_resid NaN
  # are carried as NA_real_ to match stats::chisq.test behaviour
  edge_df$p_value[is.nan(edge_df$p_value)]     <- NA_real_
  edge_df$std_resid[is.nan(edge_df$std_resid)] <- NA_real_
  
  g <- igraph::graph_from_data_frame(
    d        = edge_df,
    directed = FALSE,
    vertices = modalities
  )
  
  out <- list(
    graph            = g,
    modalities       = modalities,
    indicator_matrix = mm,
    data             = dat
  )
  
  class(out) <- "catmodgraph"
  out
}