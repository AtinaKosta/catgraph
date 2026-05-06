#' Compute phi or Cramer's V effect size for a pair of categorical variables
#'
#' Dispatches to the phi coefficient for 2x2 tables and to Cramer's V for
#' larger tables. Optionally applies the bias correction of Bergsma (2013).
#' All calculations are based on the chi-square statistic and sample size
#' returned by \code{\link{compute_assoc}}.
#'
#' @param x A factor, character, or logical vector.
#' @param y A factor, character, or logical vector of the same length as
#'   \code{x}.
#' @param corrected Logical. If \code{TRUE}, returns the bias-corrected
#'   version of Cramer's V (Bergsma, 2013) for n x m tables, or the
#'   bias-corrected phi for 2x2 tables. If \code{FALSE} (default), returns
#'   the classical estimators.
#' @param correct Logical. Yates' continuity correction for chi-square.
#'   Passed to \code{\link{compute_assoc}}. Default \code{FALSE}.
#' @param simulate_p Logical. Monte Carlo p-value simulation. Passed to
#'   \code{\link{compute_assoc}}. Default \code{FALSE}.
#' @param B Integer. Monte Carlo resamples. Default \code{2000L}.
#' @param x_name,y_name Optional character. Variable names used in warning
#'   messages to identify which pair triggered a sparsity warning. If
#'   \code{NULL} (default), the names are deparsed from the call. Primarily
#'   used internally by \code{\link{build_graph}}; users do not normally
#'   supply these.
#'
#' @return A named list with:
#' \describe{
#'   \item{effect_size}{The phi or Cramer's V value (numeric, in [0, 1] for
#'     classical; may be 0 when bias-corrected formula yields a negative
#'     value, which is clamped to 0).}
#'   \item{metric}{Character: \code{"phi"} or \code{"cramers_v"}.}
#'   \item{corrected}{Logical indicating whether bias correction was applied.}
#'   \item{type}{Observed contingency-table type: \code{"2x2"},
#'     \code{"RxC"}, or \code{"degenerate"}.}
#'   \item{statistic}{The chi-square statistic.}
#'   \item{p_value}{The p-value from the chi-square test.}
#'   \item{df}{Degrees of freedom (NA if Monte Carlo simulation was used).}
#'   \item{n}{Number of pairwise-complete observations.}
#' }
#'
#' @details
#' \strong{Phi coefficient (2x2 tables):}
#'
#' \deqn{\phi = \sqrt{\chi^2 / n}}
#'
#' This is equivalent to the Pearson correlation coefficient computed on two
#' binary variables coded 0/1 (Agresti, 2002, p. 60).
#'
#' \strong{Bias-corrected phi (Bergsma, 2013):}
#'
#' \deqn{\tilde{\phi} = \max\!\left(0,\; \phi^2 - \frac{1}{n-1}\right)^{1/2}}
#'
#' \strong{Classical Cramer's V (n x m tables):}
#'
#' \deqn{V = \sqrt{\frac{\chi^2 / n}{\min(r-1,\, c-1)}}}
#'
#' where \eqn{r} and \eqn{c} are the number of rows and columns in the
#' contingency table (Cramer, 1946).
#'
#' \strong{Bias-corrected Cramer's V (Bergsma, 2013):}
#'
#' Let \eqn{\phi^2 = \chi^2/n}. The bias-corrected estimator is computed as
#'
#' \deqn{\phi^2_{corr} = \max\!\left(0,\; \phi^2 -
#'   \frac{(r-1)(c-1)}{n-1}\right)}
#'
#' with corrected effective dimensions
#'
#' \deqn{r_{corr} = r - \frac{(r-1)^2}{n-1}, \qquad
#'       c_{corr} = c - \frac{(c-1)^2}{n-1}}
#'
#' and
#'
#' \deqn{\tilde{V} =
#'   \sqrt{\frac{\phi^2_{corr}}{\min(r_{corr}-1,\; c_{corr}-1)}}}
#'
#' whenever the denominator is positive. If the corrected denominator is
#' non-positive, the function returns 0 with a warning.
#' 
#' \strong{Comparability across table dimensions.} Both phi and Cramer's V
#' are normalised to \eqn{[0, 1]}, but this does \emph{not} make them
#' directly comparable across tables of different dimensions. Under
#' independence, the sampling distribution of V depends on table dimension
#' and sample size: for the same true association strength, V from a
#' \eqn{5 \times 5} table is not exchangeable with V from a \eqn{2 \times 2}
#' table. A V of 0.25 observed on a \eqn{5 \times 5} table and a V of 0.25
#' observed on a \eqn{2 \times 2} table represent comparable values on the
#' normalised scale but may correspond to different sampling-distribution
#' quantiles under independence. In a \code{catgraph} object, variables
#' with many more categories than their neighbours may therefore score
#' higher on Cramer's V at any given level of true dependence, which can
#' inflate their apparent centrality. Bias correction via
#' \code{corrected = TRUE} (Bergsma, 2013) partially mitigates this. See
#' the package vignette, "Methodological caveats", item 2 (mixed metrics).
#'
#' @references
#' Agresti, A. (2002). \emph{Categorical Data Analysis} (2nd ed.).
#'   John Wiley & Sons. \doi{10.1002/0471249688}
#'
#' Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
#'   \emph{Journal of the Korean Statistical Society}, 42(3), 323--328.
#'   \doi{10.1016/j.jkss.2012.10.002}
#'
#' Cramer, H. (1946). \emph{Mathematical Methods of Statistics}.
#'   Princeton University Press.
#'
#' @examples
#' set.seed(1)
#' x <- sample(c("A", "B"), 80, replace = TRUE)
#' y <- sample(c("yes", "no"), 80, replace = TRUE)
#' effect_size(x, y)
#' effect_size(x, y, corrected = TRUE)
#'
#' # Multinomial example
#' z <- sample(c("low", "mid", "high"), 80, replace = TRUE)
#' effect_size(x, z)
#' effect_size(x, z, corrected = TRUE)
#'
#' @seealso \code{\link{compute_assoc}}, \code{\link{detect_type}},
#'   \code{\link{catgraph}}
#' @export
effect_size <- function(x, y,
                        corrected  = FALSE,
                        correct    = FALSE,
                        simulate_p = FALSE,
                        B          = 2000L,
                        x_name     = NULL,
                        y_name     = NULL) {
  
  assoc <- compute_assoc(x, y,
                         correct    = correct,
                         simulate_p = simulate_p,
                         B          = B,
                         x_name     = x_name,
                         y_name     = y_name)
  
  chi2 <- assoc$statistic
  n    <- assoc$n
  type <- assoc$type
  tab  <- assoc$table
  
  if (isTRUE(type == "degenerate") || is.na(chi2) || n < 2L) {
    return(list(
      effect_size = NA_real_,
      metric      = NA_character_,
      corrected   = corrected,
      type        = type,
      statistic   = chi2,
      p_value     = assoc$p_value,
      df          = assoc$df,
      n           = n
    ))
  }
  
  r <- nrow(tab)
  cc <- ncol(tab)
  
  if (type == "2x2") {
    if (!corrected) {
      es <- sqrt(chi2 / n)
      metric <- "phi"
    } else {
      es <- sqrt(max(0, (chi2 / n) - 1 / (n - 1)))
      metric <- "phi"
    }
  } else {
    k <- min(r, cc)
    
    if (!corrected) {
      es <- sqrt((chi2 / n) / (k - 1))
      metric <- "cramers_v"
    } else {
      phi2 <- chi2 / n
      phi2_corr <- max(0, phi2 - ((r - 1) * (cc - 1)) / (n - 1))
      r_corr <- r - ((r - 1)^2) / (n - 1)
      cc_corr <- cc - ((cc - 1)^2) / (n - 1)
      denom <- min(r_corr - 1, cc_corr - 1)
      
      if (denom <= 0 || !is.finite(denom)) {
        warning(
          "Denominator for bias-corrected Cramer's V is <= 0 or non-finite (n may be too small relative to table dimensions). Returning 0.",
          call. = FALSE
        )
        es <- 0
      } else {
        es <- sqrt(phi2_corr / denom)
      }
      metric <- "cramers_v"
    }
  }
  
  es <- max(0, min(1, es))
  
  list(
    effect_size = es,
    metric      = metric,
    corrected   = corrected,
    type        = type,
    statistic   = chi2,
    p_value     = assoc$p_value,
    df          = assoc$df,
    n           = n
  )
}

#' Normalised Mutual Information for a pair of categorical variables
#'
#' Computes the symmetric Normalised Mutual Information (NMI) between two
#' categorical variables as an association measure suitable for use as an
#' edge weight in an undirected \code{catgraph}.
#'
#' NMI is defined as:
#'
#' \deqn{NMI(X, Y) = \frac{I(X; Y)}{\sqrt{H(X) \cdot H(Y)}}}
#'
#' where \eqn{I(X;Y) = H(X) + H(Y) - H(X,Y)} is the mutual information and
#' \eqn{H(\cdot)} denotes Shannon entropy (in nats, using natural logarithm).
#' The geometric-mean normalisation ensures the result lies in \eqn{[0, 1]}
#' and is symmetric: \eqn{NMI(X, Y) = NMI(Y, X)}.
#'
#' \strong{Small-sample behaviour.} Raw mutual information is upward-biased
#' when contingency tables are sparse. \code{nmi_assoc()} offers two
#' corrections via the \code{adjusted} argument:
#' \itemize{
#'   \item \code{adjusted = FALSE} (default): returns plain NMI computed
#'     from observed frequencies. Suitable when tables are well-populated.
#'   \item \code{adjusted = TRUE}: returns Adjusted Mutual Information (AMI),
#'     which subtracts the expected MI under random permutation and
#'     rescales, so that the expected value under independence is 0.
#'     Recommended for sparse tables or variables with many categories.
#' }
#'
#' \strong{Relationship to Cramer's V.} Both NMI and Cramer's V are
#' symmetric, bounded in \eqn{[0, 1]}, and equal 0 under independence.
#' They measure different things: Cramer's V quantifies departure from
#' independence relative to a chi-square null; NMI quantifies the fraction
#' of uncertainty in one variable explained by the other. The two will
#' generally agree on strong vs. weak associations but can differ on
#' variables with unequal numbers of categories or highly skewed marginals.
#'
#' @param x A factor, character, or logical vector.
#' @param y A factor, character, or logical vector of the same length as
#'   \code{x}.
#' @param adjusted Logical. If \code{TRUE}, returns Adjusted Mutual
#'   Information (AMI), which corrects for chance. Default \code{FALSE}.
#' @param x_name,y_name Optional character. Variable names used in warning
#'   messages when contingency tables are sparse. Defaults to the deparsed
#'   expressions of \code{x} and \code{y}.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{effect_size}}{NMI or AMI value, numeric in \eqn{[0, 1]}
#'     (AMI may be slightly negative due to numerical noise; clamped to 0).}
#'   \item{\code{metric}}{Character: \code{"nmi"} or \code{"ami"}.}
#'   \item{\code{adjusted}}{Logical indicating whether AMI was computed.}
#'   \item{\code{type}}{Contingency-table type: \code{"2x2"} or
#'     \code{"RxC"} or \code{"degenerate"}.}
#'   \item{\code{statistic}}{The chi-square statistic (from the underlying
#'     \code{\link{compute_assoc}} call), for reference.}
#'   \item{\code{p_value}}{Chi-square p-value, for reference.}
#'   \item{\code{df}}{Degrees of freedom.}
#'   \item{\code{n}}{Number of pairwise-complete observations.}
#' }
#'
#' @details
#' Entropy is computed using natural logarithms (nats). The choice of base
#' does not affect NMI because it cancels in the ratio. Zero-count cells
#' contribute 0 to entropy sums (the standard convention: 0 * log(0) = 0).
#'
#' The AMI formula follows Vinh, Epps & Bailey (2010), adapted for
#' two-variable contingency tables rather than clustering partitions.
#' Expected MI is computed analytically from the marginal counts using
#' the hypergeometric model.
#'
#' @references
#' Cover, T. M., & Thomas, J. A. (2006). \emph{Elements of Information
#'   Theory} (2nd ed.). Wiley. \doi{10.1002/047174882X}
#'
#' Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic
#'   measures for clusterings comparison: Variants, properties,
#'   normalisation and correction for chance.
#'   \emph{Journal of Machine Learning Research}, 11, 2837--2854.
#'   \url{https://jmlr.org/papers/v11/vinh10a.html}
#'
#' @examples
#' set.seed(1)
#' x <- sample(c("A", "B", "C"), 120, replace = TRUE)
#' y <- sample(c("yes", "no"),    120, replace = TRUE)
#' nmi_assoc(x, y)
#' nmi_assoc(x, y, adjusted = TRUE)
#'
#' @seealso \code{\link{effect_size}}, \code{\link{compute_assoc}}
#' @export
nmi_assoc <- function(x, y, adjusted = FALSE,
                      x_name = NULL, y_name = NULL) {
  
  # --- reuse compute_assoc for input validation, NA removal, and chi-square
  assoc <- compute_assoc(x, y, x_name = x_name, y_name = y_name)
  
  n    <- assoc$n
  tab  <- assoc$table
  type <- assoc$type
  
  # propagate degenerate / NA cases
  if (isTRUE(type == "degenerate") || is.null(tab) || n < 2L) {
    return(list(
      effect_size = NA_real_,
      metric      = if (adjusted) "ami" else "nmi",
      adjusted    = adjusted,
      type        = type,
      statistic   = assoc$statistic,
      p_value     = assoc$p_value,
      df          = assoc$df,
      n           = n
    ))
  }
  
  # --- entropy helper: H = -sum(p * log(p)), 0*log(0) treated as 0
  .H <- function(counts) {
    p <- counts / sum(counts)
    p <- p[p > 0]
    -sum(p * log(p))
  }
  
  # marginal and joint distributions
  p_tab <- tab / n                        # joint proportions
  p_x   <- rowSums(p_tab)                 # marginal of x
  p_y   <- colSums(p_tab)                 # marginal of y
  
  hx  <- .H(rowSums(tab))
  hy  <- .H(colSums(tab))
  hxy <- .H(as.vector(tab))
  
  mi  <- hx + hy - hxy                   # mutual information (nats)
  mi  <- max(0, mi)                       # guard against tiny negatives
  
  denom_nmi <- sqrt(hx * hy)
  
  # degenerate: one variable is constant after NA removal (entropy = 0)
  if (denom_nmi < .Machine$double.eps) {
    es <- 0
    metric <- if (adjusted) "ami" else "nmi"
    return(list(
      effect_size = es,
      metric      = metric,
      adjusted    = adjusted,
      type        = type,
      statistic   = assoc$statistic,
      p_value     = assoc$p_value,
      df          = assoc$df,
      n           = n
    ))
  }
  
  if (!adjusted) {
    # --- plain NMI
    es     <- mi / denom_nmi
    es     <- min(1, max(0, es))
    metric <- "nmi"
    
  } else {
    # --- AMI: subtract expected MI under the hypergeometric model
    # (Vinh et al. 2010, equation 2)
    r_counts <- rowSums(tab)   # marginal counts for x
    c_counts <- colSums(tab)   # marginal counts for y
    R <- length(r_counts)
    C <- length(c_counts)
    
    # Expected MI: sum over all (i,j) cells
    emi <- 0
    for (i in seq_len(R)) {
      for (j in seq_len(C)) {
        ai <- r_counts[i]
        bj <- c_counts[j]
        # n_ij ranges from max(0, ai+bj-n) to min(ai, bj)
        nij_min <- max(1L, ai + bj - n)   # must be >= 1 for log to be defined
        nij_max <- min(ai, bj)
        if (nij_min > nij_max) next
        for (nij in seq(nij_min, nij_max)) {
          # log-prob of this cell count under hypergeometric
          log_p <- lchoose(ai, nij) + lchoose(n - ai, bj - nij) - lchoose(n, bj)
          term  <- exp(log_p) * (nij / n) * (log(nij) - log(ai) - log(bj) + log(n))
          emi   <- emi + term
        }
      }
    }
    emi <- max(0, emi)
    
    # normaliser for AMI uses the same sqrt(H(X)*H(Y)) denominator as NMI
    # (one of several valid choices; keeps AMI comparable with NMI)
    ami_num   <- mi - emi
    ami_denom <- denom_nmi - emi
    if (abs(ami_denom) < .Machine$double.eps) {
      es <- 0
    } else {
      es <- ami_num / ami_denom
    }
    es     <- min(1, max(0, es))
    metric <- "ami"
  }
  
  list(
    effect_size = es,
    metric      = metric,
    adjusted    = adjusted,
    type        = type,
    statistic   = assoc$statistic,
    p_value     = assoc$p_value,
    df          = assoc$df,
    n           = n
  )
}
  
#' Bayesian Cramér's V for a pair of categorical variables
#' Computes a Bayesian estimate of Cramér's V by applying a symmetric
#' Dirichlet prior to the contingency table cell counts before computing
#' the association measure. This shrinks edge weights toward zero for
#' sparse tables, producing more stable estimates than classical or
#' bias-corrected Cramér's V when expected cell frequencies are small.
#'
#' @param x A factor, character, or logical vector.
#' @param y A factor, character, or logical vector of the same length as
#'   \code{x}.
#' @param alpha Numeric. Dirichlet prior concentration parameter added to
#'   each cell count before computing the association. Must be > 0.
#'   Default \code{0.5} (Jeffreys prior). Use \code{alpha = 1} for the
#'   Laplace (uniform) prior.
#' @param x_name,y_name Optional character. Variable names used in warning
#'   messages when contingency tables are sparse. Defaults to the deparsed
#'   expressions of \code{x} and \code{y}.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{effect_size}}{Bayesian Cramér's V, numeric in
#'     \eqn{[0, 1]}.}
#'   \item{\code{metric}}{Character: \code{"bayesian_cramers_v"}.}
#'   \item{\code{alpha}}{The prior concentration used.}
#'   \item{\code{type}}{Contingency-table type: \code{"2x2"},
#'     \code{"RxC"}, or \code{"degenerate"}.}
#'   \item{\code{statistic}}{The chi-square statistic computed on the
#'     \strong{smoothed} table, for reference.}
#'   \item{\code{p_value}}{Chi-square p-value from the
#'     \strong{original} (unsmoothed) table. Smoothing is applied only
#'     to the effect-size estimate, not to the test.}
#'   \item{\code{df}}{Degrees of freedom.}
#'   \item{\code{n}}{Number of pairwise-complete observations
#'     (unsmoothed).}
#' }
#'
#' @details
#' \strong{Dirichlet smoothing.} Under a symmetric Dirichlet(\eqn{\alpha})
#' prior on the \eqn{r \times c} cell probability vector, the posterior
#' mean of each cell probability is:
#'
#' \deqn{\hat{p}_{ij} = \frac{n_{ij} + \alpha}{n + \alpha \cdot r \cdot c}}
#'
#' Cramér's V is then computed from the chi-square statistic derived from
#' these smoothed proportions rather than from the raw counts. This is
#' equivalent to computing Cramér's V on a pseudo-count table
#' \eqn{\tilde{n}_{ij} = n_{ij} + \alpha} with effective sample size
#' \eqn{\tilde{n} = n + \alpha \cdot r \cdot c}.
#'
#' \strong{Jeffreys prior (\eqn{\alpha = 0.5}).} This is the standard
#' non-informative choice for categorical data. It corresponds to adding
#' half a pseudocount to each cell, which stabilises the chi-square
#' statistic for sparse tables without materially distorting the estimate
#' when tables are well-populated.
#'
#' \strong{Relationship to classical Cramér's V.} As \eqn{n \to \infty},
#' the smoothed estimate converges to the classical estimator because the
#' pseudocounts \eqn{\alpha} become negligible relative to \eqn{n}.
#' For large samples (such as the Titanic dataset with n = 2201) the
#' difference is therefore very small. The practical advantage of the
#' Bayesian estimator appears on small samples or sparse contingency
#' tables.
#'
#' \strong{p-value.} The p-value is taken from the unsmoothed chi-square
#' test (via \code{\link{compute_assoc}}). This is intentional: smoothing
#' inflates the effective sample size and would otherwise produce
#' anti-conservative p-values. Users who want a fully Bayesian decision
#' criterion should use posterior credible intervals (not yet implemented)
#' rather than the p-value.
#'
#' @references
#' Good, I. J. (1965). \emph{The Estimation of Probabilities: An Essay on
#'   Modern Bayesian Methods}. MIT Press.
#'
#' Agresti, A. (2002). \emph{Categorical Data Analysis} (2nd ed.).
#'   Wiley. \doi{10.1002/0471249688}
#'
#'
#' @examples
#' set.seed(1)
#' x <- sample(c("A", "B", "C"), 120, replace = TRUE)
#' y <- sample(c("yes", "no"),    120, replace = TRUE)
#'
#' # Classical
#' effect_size(x, y)$effect_size
#'
#' # Bayesian (Jeffreys prior)
#' bayesian_cramers_v(x, y)$effect_size
#'
#' # Bayesian (Laplace prior)
#' bayesian_cramers_v(x, y, alpha = 1)$effect_size
#'
#' # On a sparse table: Bayesian estimate is more stable
#' x_sparse <- sample(c("A","B","C","D"), 20, replace = TRUE)
#' y_sparse <- sample(c("P","Q","R","S"), 20, replace = TRUE)
#' effect_size(x_sparse, y_sparse)$effect_size
#' bayesian_cramers_v(x_sparse, y_sparse)$effect_size
#'
#' @seealso \code{\link{effect_size}}, \code{\link{nmi_assoc}},
#'   \code{\link{compute_assoc}}
#' @export
bayesian_cramers_v <- function(x, y, alpha = 0.5,
                               x_name = NULL, y_name = NULL) {

  # --- validate alpha -------------------------------------------------------
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha <= 0) {
    stop("`alpha` must be a single positive number.", call. = FALSE)
  }

  # --- reuse compute_assoc for input validation, NA removal, table ----------
  assoc <- compute_assoc(x, y, x_name = x_name, y_name = y_name)

  n    <- assoc$n
  tab  <- assoc$table
  type <- assoc$type

  # propagate degenerate / NA cases
  if (isTRUE(type == "degenerate") || is.null(tab) || n < 2L) {
    return(list(
      effect_size = NA_real_,
      metric      = "bayesian_cramers_v",
      alpha       = alpha,
      type        = type,
      statistic   = assoc$statistic,
      p_value     = assoc$p_value,
      df          = assoc$df,
      n           = n
    ))
  }

  r  <- nrow(tab)
  cc <- ncol(tab)

  # --- Dirichlet smoothing --------------------------------------------------
  # Pseudo-count table: add alpha to every cell
  tab_smooth <- tab + alpha

  # Effective sample size after smoothing
  n_smooth <- n + alpha * r * cc

  # Smoothed expected counts under independence:
  # E_ij = (row_i_total * col_j_total) / n_smooth
  row_tots <- rowSums(tab_smooth)
  col_tots <- colSums(tab_smooth)
  expected_smooth <- outer(row_tots, col_tots) / n_smooth

  # Chi-square on smoothed table
  chi2_smooth <- sum((tab_smooth - expected_smooth)^2 / expected_smooth)

  # Cramér's V from smoothed chi-square
  k <- min(r, cc)
  phi2_smooth <- chi2_smooth / n_smooth

  if (k <= 1) {
    # degenerate after smoothing — should not happen but guard anyway
    warning(
      "Smoothed table has effective min dimension <= 1; returning 0.",
      call. = FALSE
    )
    es <- 0
  } else {
    es <- sqrt(phi2_smooth / (k - 1))
    es <- min(1, max(0, es))
  }

  list(
    effect_size = es,
    metric      = "bayesian_cramers_v",
    alpha       = alpha,
    type        = type,
    statistic   = chi2_smooth,        # chi-square on smoothed table
    p_value     = assoc$p_value,      # p-value from unsmoothed test
    df          = assoc$df,
    n           = n                   # original (unsmoothed) n
  )
}