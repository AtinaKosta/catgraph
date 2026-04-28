#' Joint categorical distribution diagnostic across groups
#'
#' Runs a descriptive cross-group diagnostic for categorical data: marginal
#' comparisons per variable, plus modality-level omnibus and edge-wise
#' association-structure comparisons. Intended as the user-facing entry point for cross-group
#' diagnostics in meta-analysis, multi-site studies, and repeated
#' cross-sectional surveys.
#'
#' @param data A data frame of categorical variables plus a grouping column.
#' @param group Character. Name of the grouping column in \code{data}.
#'   Must be categorical with at least 2 levels.
#' @param variables Optional character vector. Names of columns to include
#'   in the joint analysis. Default \code{NULL} uses all columns in
#'   \code{data} other than \code{group}.
#' @param n_perm Integer. Permutations for the omnibus test, passed to
#'   \code{\link{test_modality_graph_equality}}. Default \code{500L}.
#' @param n_perm_edge Integer. Permutations for the edge-wise post-hoc
#'   test. Default equal to \code{n_perm}.
#' @param alpha Numeric. Significance level reported for all three
#'   testing families. Default \code{0.05}.
#' @param strata Optional stratification vector passed through to the
#'   omnibus and edge-wise tests. Length must equal \code{nrow(data)}.
#' @param run_edgewise Logical. If \code{TRUE} (default), runs
#'   \code{\link{test_modality_edge_differences}} for each pair whose
#'   omnibus test rejects at \code{alpha}. Set to \code{FALSE} to skip
#'   the post-hoc step (faster, useful when you only want the balance
#'   summary).
#' @param seed Optional integer seed.
#' @param verbose Logical. If \code{TRUE} (default), prints progress.
#'
#' @return An object of class \code{jointbalance} with components:
#'   \describe{
#'     \item{\code{group}}{Name of the grouping variable.}
#'     \item{\code{group_levels}}{Character vector of group levels.}
#'     \item{\code{variables}}{Character vector of analysed variables.}
#'     \item{\code{marginal}}{Data frame, one row per variable:
#'       \code{variable}, \code{chisq}, \code{df}, \code{p_value},
#'       \code{p_adjusted}, \code{cramers_v}. Adjustment is BH across
#'       variables.}
#'     \item{\code{pairwise_omnibus}}{Data frame, one row per group
#'       pair: \code{group_x}, \code{group_y}, \code{n_x}, \code{n_y},
#'       \code{statistic}, \code{observed}, \code{p_value},
#'       \code{p_bonferroni}. Adjustment is Bonferroni across pairs.}
#'     \item{\code{pairwise_edgewise}}{Named list of
#'       \code{catmodedgetest} objects, one per pair that rejected
#'       the omnibus (empty if none did or \code{run_edgewise = FALSE}).
#'       Names are \code{"\{group_x\}_vs_\{group_y\}"}.}
#'     \item{\code{modality_graphs}}{Named list of \code{catmodgraph}
#'       objects, one per group level.}
#'     \item{\code{alpha}}{The alpha used.}
#'     \item{\code{call}}{The matched call.}
#'   }
#'
#' @details
#' \strong{What this function reports, in plain terms.}
#' \enumerate{
#'   \item Marginal balance (the classical "Table 1" check): is each
#'     variable distributed the same across groups, one variable at a
#'     time?
#'   \item Joint balance (the modality network test): is the full
#'     joint categorical structure the same across groups?
#'   \item Where the joint structure differs (the edge-wise post-hoc):
#'     which modality-pair associations account for any joint
#'     disagreement?
#' }
#' Marginal similarity is \emph{not sufficient} to establish similarity of
#' joint categorical association structure. Two groups with identical marginals can still differ in
#' joint structure --- for example, if variable A and variable B are
#' positively associated in group 1 and negatively associated in
#' group 2, both groups will show the same one-way marginals for A
#' and for B, but the joint distribution will differ. This is exactly
#' the scenario the modality-layer test is designed to catch, and it
#' is the reason a joint-structure diagnostic can reveal discrepancies
#' that marginal-only Table 1 summaries may miss.
#' 
#' \strong{Interpretation.} The diagnostic compares observed categorical
#' distributions and modality-level association patterns across groups. It
#' does not establish exchangeability, causal comparability, or absence of
#' residual confounding. Results should be interpreted as evidence of
#' distributional and associational discrepancies, not as proof that groups
#' are or are not analytically interchangeable.
#'
#' \strong{Testing families and multiplicity.} Three separate
#' multiple-testing corrections are applied:
#' \itemize{
#'   \item Marginal tests: Benjamini-Hochberg FDR across the \eqn{p}
#'     variables. Exploratory; the marginal panel of the plot is
#'     meant to triage, not to make definitive claims.
#'   \item Omnibus tests across group pairs: Bonferroni across the
#'     \eqn{\binom{k}{2}} pairs. Conservative by design, because
#'     these pairwise tests are logically nested (rejecting
#'     "wave 1 vs wave 2" and "wave 1 vs wave 3" both implicate
#'     wave 1).
#'   \item Edge-wise post-hoc within each pair: Benjamini-Hochberg
#'     across edges, as done by
#'     \code{\link{test_modality_edge_differences}}. This is a
#'     conditional family, valid under the closed-testing principle
#'     only if the omnibus rejected.
#' }
#' These are \strong{three separate testing families}, not one big
#' correction. The print method makes this distinction explicit.
#'
#' \strong{When not to use this.} \code{joint_balance()} is a
#' descriptive diagnostic for fully categorical data. It is not a
#' causal adjustment tool and does not replace propensity-score or
#' entropy-balancing workflows. For continuous covariates use the
#' \pkg{cobalt} or \pkg{tableone} packages; this function silently
#' coerces numeric columns to factors if they appear in
#' \code{variables}, which is almost never what you want.
#'
#' @references
#' The marginal/joint distinction is a long-standing one in
#' log-linear modelling; see Agresti (2013), \emph{Categorical Data
#' Analysis} (3rd ed.), Wiley, chapter 9, for homogeneity-of-
#' association tests built from three-way log-linear models. The
#' permutation approach used here is an alternative that avoids the
#' large-sample assumptions of the likelihood-ratio test.
#'
#' @examples
#' \donttest{
#' data(survey_health)
#' jb <- joint_balance(
#'   survey_health, group = "sex",
#'   n_perm = 200, n_perm_edge = 200, seed = 1, verbose = FALSE
#' )
#' jb
#' summary(jb)
#' plot(jb)
#' }
#'
#' @seealso \code{\link{test_modality_graph_equality}},
#'   \code{\link{test_modality_edge_differences}},
#'   \code{\link{plot_modality_difference}},
#'   \code{\link{compare_modality_graphs}}
#' @importFrom stats chisq.test p.adjust
#' @importFrom utils combn
#' @export
joint_balance <- function(data,
                          group,
                          variables    = NULL,
                          n_perm       = 500L,
                          n_perm_edge  = NULL,
                          alpha        = 0.05,
                          strata       = NULL,
                          run_edgewise = TRUE,
                          seed         = NULL,
                          verbose      = TRUE) {
  
  mc <- match.call()
  
  # ---- Input validation -----------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1L ||
      !group %in% names(data)) {
    stop("`group` must be the name of a column in `data`.", call. = FALSE)
  }
  if (!is.null(variables)) {
    if (!is.character(variables) || length(variables) < 2L) {
      stop("`variables` must be a character vector of length >= 2.",
           call. = FALSE)
    }
    missing_vars <- setdiff(variables, names(data))
    if (length(missing_vars) > 0L) {
      stop("Variables not found in `data`: ",
           paste(missing_vars, collapse = ", "), ".", call. = FALSE)
    }
    if (group %in% variables) {
      stop("`group` must not appear in `variables`.", call. = FALSE)
    }
  } else {
    variables <- setdiff(names(data), group)
    if (length(variables) < 2L) {
      stop("`data` must contain at least 2 non-group variables.",
           call. = FALSE)
    }
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.logical(run_edgewise) || length(run_edgewise) != 1L ||
      is.na(run_edgewise)) {
    stop("`run_edgewise` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(n_perm_edge)) n_perm_edge <- n_perm
  
  # ---- Prepare data --------------------------------------------------------
  grp_col <- data[[group]]
  if (!is.factor(grp_col)) grp_col <- factor(grp_col)
  grp_col <- droplevels(grp_col)
  grp_levels <- levels(grp_col)
  
  if (length(grp_levels) < 2L) {
    stop("Grouping variable `", group,
         "` must have at least 2 levels after dropping NAs.", call. = FALSE)
  }
  
  dat <- data[, c(group, variables), drop = FALSE]
  keep_row <- !is.na(grp_col)
  dat      <- dat[keep_row, , drop = FALSE]
  grp_col  <- grp_col[keep_row]
  if (!is.null(strata)) {
    if (length(strata) != nrow(data)) {
      stop("`strata` must have length nrow(data) = ", nrow(data), ".",
           call. = FALSE)
    }
    strata <- strata[keep_row]
  }
  
  # Coerce variables to factors; warn on numeric coercion
  for (v in variables) {
    if (is.numeric(dat[[v]])) {
      warning("Variable `", v, "` is numeric; coercing to factor. ",
              "joint_balance() is a categorical diagnostic. ",
              "If `", v, "` is genuinely continuous, check balance with ",
              "tableone::CreateTableOne() or cobalt::bal.tab() instead.",
              call. = FALSE)
    }
    if (!is.factor(dat[[v]])) dat[[v]] <- factor(dat[[v]])
  }
  
  # ---- Step 1: marginal chi-square per variable ----------------------------
  if (verbose) cat("Running marginal tests (", length(variables),
                   " variables)...\n", sep = "")
  
  marg_rows <- lapply(variables, function(v) {
    tab <- table(dat[[v]], grp_col)
    # Guard against empty rows/cols after subsetting
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
    if (any(dim(tab) < 2L)) {
      return(data.frame(variable = v, chisq = NA_real_, df = NA_integer_,
                        p_value = NA_real_, cramers_v = NA_real_,
                        stringsAsFactors = FALSE))
    }
    chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    n   <- sum(tab)
    k   <- min(dim(tab))
    cv  <- sqrt(as.numeric(chi$statistic) / (n * (k - 1)))
    data.frame(
      variable  = v,
      chisq     = unname(chi$statistic),
      df        = unname(chi$parameter),
      p_value   = unname(chi$p.value),
      cramers_v = cv,
      stringsAsFactors = FALSE
    )
  })
  marginal <- do.call(rbind, marg_rows)
  marginal$p_adjusted <- stats::p.adjust(marginal$p_value, method = "BH")
  
  # ---- Step 2: build one modality graph per group level --------------------
  if (verbose) cat("Building modality graphs per group level...\n")
  
  modality_graphs <- lapply(grp_levels, function(lv) {
    sub_df <- dat[grp_col == lv, variables, drop = FALSE]
    if (nrow(sub_df) < 2L) {
      stop("Group level '", lv, "' has < 2 rows; cannot build modality graph.",
           call. = FALSE)
    }
    suppressWarnings(build_modality_graph(sub_df))
  })
  names(modality_graphs) <- grp_levels
  
  # ---- Step 3: pairwise omnibus tests --------------------------------------
  pair_idx <- utils::combn(length(grp_levels), 2L, simplify = FALSE)
  n_pairs  <- length(pair_idx)
  
  if (verbose) cat("Running ", n_pairs, " pairwise omnibus test(s)...\n",
                   sep = "")
  
  omni_rows <- vector("list", n_pairs)
  for (p in seq_len(n_pairs)) {
    i <- pair_idx[[p]][1L]; j <- pair_idx[[p]][2L]
    gx <- grp_levels[i];    gy <- grp_levels[j]
    
    if (verbose) cat("  [", p, "/", n_pairs, "] ", gx, " vs ", gy, "\n",
                     sep = "")
    
    mask_pair <- grp_col %in% c(gx, gy)
    strata_p  <- if (!is.null(strata)) strata[mask_pair] else NULL
    
    tst <- test_modality_graph_equality(
      modality_graphs[[gx]], modality_graphs[[gy]],
      n_perm    = n_perm,
      strata    = strata_p,
      seed      = seed,
      verbose   = FALSE
    )
    
    omni_rows[[p]] <- data.frame(
      group_x   = gx,
      group_y   = gy,
      n_x       = tst$n_x,
      n_y       = tst$n_y,
      statistic = tst$statistic,
      observed  = tst$observed,
      p_value   = tst$p_value,
      stringsAsFactors = FALSE
    )
  }
  pairwise_omnibus <- do.call(rbind, omni_rows)
  pairwise_omnibus$p_bonferroni <- pmin(pairwise_omnibus$p_value * n_pairs, 1)
  
  # ---- Step 4: edge-wise post-hoc for rejecting pairs ---------------------
  pairwise_edgewise <- list()
  if (run_edgewise) {
    reject_idx <- which(pairwise_omnibus$p_bonferroni < alpha)
    if (length(reject_idx) > 0L && verbose) {
      cat("Running edge-wise post-hoc on ", length(reject_idx),
          " rejecting pair(s)...\n", sep = "")
    }
    for (r in reject_idx) {
      gx <- pairwise_omnibus$group_x[r]
      gy <- pairwise_omnibus$group_y[r]
      mask_pair <- grp_col %in% c(gx, gy)
      strata_p  <- if (!is.null(strata)) strata[mask_pair] else NULL
      
      et <- test_modality_edge_differences(
        modality_graphs[[gx]], modality_graphs[[gy]],
        n_perm   = n_perm_edge,
        edges    = "union",
        p_adjust = "BH",
        strata   = strata_p,
        seed     = seed,
        verbose  = FALSE
      )
      pairwise_edgewise[[paste0(gx, "_vs_", gy)]] <- et
    }
  }
  
  # ---- Assemble -----------------------------------------------------------
  structure(
    list(
      group             = group,
      group_levels      = grp_levels,
      variables         = variables,
      marginal          = marginal,
      pairwise_omnibus  = pairwise_omnibus,
      pairwise_edgewise = pairwise_edgewise,
      modality_graphs   = modality_graphs,
      alpha             = alpha,
      call              = mc
    ),
    class = "jointbalance"
  )
}