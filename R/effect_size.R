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
