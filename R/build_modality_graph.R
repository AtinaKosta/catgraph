#' Build a modality-level association graph
#'
#' Constructs a modality-level association network in which each node is a
#' category level, encoded as \code{"variable=level"}, and each edge represents
#' the pairwise association between two modalities from different variables.
#'
#' Edge weights (\code{weight}) reflect association strength between modality
#' pairs according to the chosen \code{method}. Since every modality pair
#' reduces to a binary 2x2 indicator table, all four metric families are
#' well-defined: absolute phi (frequentist), NMI on the 2x2 table
#' (information-theoretic), Dirichlet-smoothed phi (Bayesian), or
#' AMI (chance-corrected information-theoretic). The signed phi coefficient is always stored as
#' \code{phi_signed} regardless of \code{method}, and the signed standardised
#' residual is always stored as \code{std_resid}: positive values indicate
#' co-occurrence above expectation (attraction), negative values indicate
#' co-occurrence below expectation (repulsion). Edge weights therefore always
#' reflect association magnitude, and direction is always accessible via
#' \code{std_resid}.
#'
#' @param data A data frame of categorical variables.
#' @param method Character. Association metric for edge weights. One of
#'   \code{"cramers_v"} (default, absolute phi on each 2x2 table),
#'   \code{"cramers_v_corrected"} (bias-corrected phi, Bergsma 2013),
#'   \code{"nmi"} (Normalised Mutual Information on each 2x2 table),
#'   \code{"ami"} (Adjusted Mutual Information), or
#'   \code{"bayesian_cramers_v"} (Dirichlet-smoothed phi). Must match
#'   the \code{method} used in the corresponding \code{\link{catgraph}}
#'   call for consistency.
#' @param alpha Numeric. Dirichlet prior concentration for
#'   \code{method = "bayesian_cramers_v"}. Default \code{0.5}
#'   (Jeffreys prior). Ignored for all other methods.
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
#'     \item{\code{method}}{Character string recording which association
#'       metric was used.}
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
                                 method    = "cramers_v",
                                 alpha     = 0.5,
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
  
  # --- resolve method -------------------------------------------------------
  valid_methods <- c("cramers_v", "cramers_v_corrected", "nmi", "ami",
                     "bayesian_cramers_v")
  if (!method %in% valid_methods) {
    stop(
      "`method` must be one of: ",
      paste(valid_methods, collapse = ", "), ".",
      call. = FALSE
    )
  }
  use_corrected <- method == "cramers_v_corrected"
  use_nmi       <- method %in% c("nmi", "ami")
  use_adjusted  <- method == "ami"
  use_bayesian  <- method == "bayesian_cramers_v"
  
  if (use_bayesian && (!is.numeric(alpha) || length(alpha) != 1L ||
                       is.na(alpha) || alpha <= 0)) {
    stop("`alpha` must be a single positive number.", call. = FALSE)
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
  
  # -- Vectorised 2x2 statistics (always computed — phi and std_resid are
  #    stored regardless of method, as they carry the direction signal) ------
  # Co-occurrence matrix: A[i, j] = sum(mm[, i] * mm[, j])
  A <- crossprod(mm)                    # n_mod x n_mod, symmetric
  col_tot <- diag(A)                    # n_i: total count of each modality
  
  # Outer products give the required marginal combinations
  n_minus    <- n_rows - col_tot
  nij        <- tcrossprod(col_tot)     # n_i * n_j
  n_minus_ij <- tcrossprod(n_minus)     # (n - n_i) * (n - n_j)
  
  # phi numerator: n * A - n_i * n_j
  num      <- n_rows * A - nij
  # phi denominator: sqrt(n_i * n_j * (n - n_i) * (n - n_j))
  denom_sq <- nij * n_minus_ij
  denom    <- suppressWarnings(sqrt(denom_sq))
  
  phi_signed <- num / denom             # NaN where denom == 0
  phi_abs    <- abs(phi_signed)
  
  # Chi-square statistic on 1 df: n * phi^2
  chi_sq <- n_rows * phi_signed^2
  p_mat  <- stats::pchisq(chi_sq, df = 1, lower.tail = FALSE)
  
  # Standardised Pearson residual: phi_signed * sqrt(n)
  std_resid_mat <- phi_signed * sqrt(n_rows)
  
  # -- Compute weight matrix depending on method ----------------------------
  if (use_bayesian) {
    # Bayesian path: Dirichlet-smoothed phi on each 2x2 table.
    #
    # For modality pair (i, j) the four smoothed cell counts are:
    #   n11_s = A[i,j] + alpha          (both present)
    #   n10_s = col_tot[i] - A[i,j] + alpha   (i present, j absent)
    #   n01_s = col_tot[j] - A[i,j] + alpha   (i absent,  j present)
    #   n00_s = n - col_tot[i] - col_tot[j] + A[i,j] + alpha
    #
    # Smoothed total: n_s = n + 4 * alpha  (4 cells in every 2x2 table)
    # Smoothed row/col marginals:
    #   row1_s = n11_s + n10_s = col_tot[i] + 2*alpha
    #   row0_s = n01_s + n00_s = (n - col_tot[i]) + 2*alpha
    #   col1_s = n11_s + n01_s = col_tot[j] + 2*alpha
    #   col0_s = n10_s + n00_s = (n - col_tot[j]) + 2*alpha
    #
    # phi_smooth = (n_s * n11_s - row1_s * col1_s) /
    #              sqrt(row1_s * col1_s * row0_s * col0_s)
    #
    # All operations vectorised over the full n_mod x n_mod matrix.
    
    n_s    <- n_rows + 4 * alpha
    n11_s  <- A + alpha
    # row marginals of smoothed table
    row1_s <- outer(col_tot + 2 * alpha,   rep(1,     n_mod))
    row0_s <- outer(n_rows - col_tot + 2 * alpha, rep(1, n_mod))
    # col marginals of smoothed table
    col1_s <- outer(rep(1, n_mod), col_tot + 2 * alpha)
    col0_s <- outer(rep(1, n_mod), n_rows - col_tot + 2 * alpha)
    
    num_s    <- n_s * n11_s - row1_s * col1_s
    denom_s  <- suppressWarnings(sqrt(row1_s * col1_s * row0_s * col0_s))
    
    phi_smooth <- num_s / denom_s         # NaN where denom = 0
    weight_mat <- matrix(
      pmin(1, pmax(0, as.vector(abs(phi_smooth)))),
      nrow = n_mod, ncol = n_mod
    )
    diag(weight_mat) <- 0
    
  } else if (!use_nmi) {
    # Frequentist path: absolute phi (classical or bias-corrected)
    if (!use_corrected) {
      weight_mat <- phi_abs
    } else {
      # Bergsma (2013) bias correction applied cell-wise to the 2x2 phi^2
      phi2_corr <- pmax(0, phi_signed^2 - 1 / (n_rows - 1))
      weight_mat <- sqrt(phi2_corr)
    }
  } else {
    # Information-theoretic path: NMI (or AMI) on each 2x2 table
    # All quantities derivable from A (co-occurrence) and col_tot (marginals).
    # For the 2x2 table formed by modalities i and j:
    #   cell (1,1) = A[i,j]         cell (1,0) = col_tot[i] - A[i,j]
    #   cell (0,1) = col_tot[j] - A[i,j]   cell (0,0) = n - col_tot[i] - col_tot[j] + A[i,j]
    #
    # We compute NMI vectorised over all pairs using matrix arithmetic.
    
    # Safe entropy helper operating on a matrix of counts (one row = one cell)
    # h(p) = -p * log(p), with 0*log(0) = 0
    .safe_ent <- function(x, total) {
      p <- x / total
      out <- -p * log(p)
      out[x <= 0] <- 0
      out
    }
    
    # Four 2x2 cell count matrices (each n_mod x n_mod)
    n11 <- A
    n10 <- outer(col_tot, rep(1, n_mod)) - A
    n01 <- outer(rep(1, n_mod), col_tot) - A
    n00 <- n_rows - n10 - n01 - n11
    
    # Joint entropy H(X,Y)
    hxy_mat <- .safe_ent(n11, n_rows) +
      .safe_ent(n10, n_rows) +
      .safe_ent(n01, n_rows) +
      .safe_ent(n00, n_rows)
    
    # Marginal entropies
    hx_vec <- .safe_ent(col_tot, n_rows) +
      .safe_ent(n_rows - col_tot, n_rows)
    hx_mat <- outer(hx_vec, rep(1, n_mod))
    hy_mat <- outer(rep(1, n_mod), hx_vec)
    
    # Mutual information
    mi_mat <- matrix(
      pmax(0, as.vector(hx_mat + hy_mat - hxy_mat)),
      nrow = n_mod, ncol = n_mod
    )
    
    denom_nmi <- sqrt(hx_mat * hy_mat)
    
    if (!use_adjusted) {
      nmi_mat <- mi_mat / denom_nmi
      nmi_mat[denom_nmi < .Machine$double.eps] <- 0
      weight_mat <- matrix(
        pmin(1, pmax(0, as.vector(nmi_mat))),
        nrow = n_mod, ncol = n_mod
      )
      diag(weight_mat) <- 0
      
    } else {
      ami_mat <- matrix(0, n_mod, n_mod)
      idx_pairs <- which(upper.tri(A), arr.ind = TRUE)
      for (k in seq_len(nrow(idx_pairs))) {
        ii <- idx_pairs[k, 1L]
        jj <- idx_pairs[k, 2L]
        ai <- col_tot[ii]
        bj <- col_tot[jj]
        nij_min <- max(1L, ai + bj - n_rows)
        nij_max <- min(ai, bj)
        emi <- 0
        if (nij_min <= nij_max) {
          for (nij_val in seq(nij_min, nij_max)) {
            log_p <- lchoose(ai, nij_val) +
              lchoose(n_rows - ai, bj - nij_val) -
              lchoose(n_rows, bj)
            term  <- exp(log_p) * (nij_val / n_rows) *
              (log(nij_val) - log(ai) - log(bj) + log(n_rows))
            emi   <- emi + term
          }
        }
        emi <- max(0, emi)
        mi_ij   <- mi_mat[ii, jj]
        dn      <- denom_nmi[ii, jj]
        ami_val <- if (abs(dn - emi) < .Machine$double.eps) 0 else
          (mi_ij - emi) / (dn - emi)
        ami_val <- min(1, max(0, ami_val))
        ami_mat[ii, jj] <- ami_val
        ami_mat[jj, ii] <- ami_val
      }
      weight_mat <- matrix(as.vector(ami_mat), nrow = n_mod, ncol = n_mod)
      diag(weight_mat) <- 0
    }
  }
  
  # Zero same-variable pairs in weight_mat to prevent any leakage into edges.
  # Build a logical mask: TRUE wherever both nodes belong to the same variable.
  same_var_mask <- outer(var_of_node, var_of_node, "==")
  weight_mat[same_var_mask] <- 0
  # Reconstruct matrix dimensions — logical indexing can strip the matrix class
  weight_mat <- matrix(as.vector(weight_mat), nrow = n_mod, ncol = n_mod)
  diag(weight_mat) <- 0
  
  # -- Assemble edge table from upper triangle, excluding same-variable pairs
  idx     <- which(upper.tri(A), arr.ind = TRUE)
  i_idx   <- idx[, 1L]
  j_idx   <- idx[, 2L]
  
  same_var <- var_of_node[i_idx] == var_of_node[j_idx]
  keep     <- !same_var
  
  i_idx <- i_idx[keep]
  j_idx <- j_idx[keep]
  
  edge_df <- data.frame(
    from       = node_names[i_idx],
    to         = node_names[j_idx],
    weight     = weight_mat[cbind(i_idx, j_idx)],
    phi_signed = phi_signed[cbind(i_idx, j_idx)],
    p_value    = p_mat[cbind(i_idx, j_idx)],
    std_resid  = std_resid_mat[cbind(i_idx, j_idx)],
    stringsAsFactors = FALSE
  )
  
  # Drop degenerate pairs
  edge_df <- edge_df[!is.na(edge_df$weight), , drop = FALSE]
  
  if (nrow(edge_df) == 0L) {
    stop("No cross-variable modality associations could be computed.",
         call. = FALSE)
  }
  
  edge_df$p_value[is.nan(edge_df$p_value)]     <- NA_real_
  edge_df$std_resid[is.nan(edge_df$std_resid)] <- NA_real_
  
  g <- igraph::graph_from_data_frame(
    d        = edge_df,
    directed = FALSE,
    vertices = modalities
  )
  # Collapse any duplicate edges that can arise when weight_mat is fully
  # populated (NMI / AMI paths) — keep the maximum weight of any duplicates
  if (igraph::any_multiple(g)) {
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                          edge.attr.comb = list(weight    = "max",
                                                phi_signed = "first",
                                                p_value   = "first",
                                                std_resid  = "first"))
  }
  
  out <- list(
    graph            = g,
    modalities       = modalities,
    indicator_matrix = mm,
    data             = dat,
    method           = method,
    alpha            = if (use_bayesian) alpha else NA_real_
  )
  
  class(out) <- "catmodgraph"
  out
}