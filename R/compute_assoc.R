#' Compute chi-square association between two categorical variables
#'
#' Performs a Pearson chi-square test of independence on the pairwise-complete
#' observations of two categorical variables. Returns the test statistic,
#' p-value, degrees of freedom, and the contingency table used.
#'
#' @param x A factor, character, or logical vector.
#' @param y A factor, character, or logical vector of the same length as
#'   \code{x}.
#' @param correct Logical. Whether to apply Yates' continuity correction for
#'   2x2 tables. Default \code{FALSE} to keep results comparable with effect
#'   size formulas that do not incorporate the correction
#'   (Agresti, 2002, p. 77).
#' @param simulate_p Logical. If \code{TRUE}, the p-value is estimated by
#'   Monte Carlo simulation. Recommended when expected cell frequencies are
#'   small. Default \code{FALSE}.
#' @param B Integer. Monte Carlo resamples when \code{simulate_p = TRUE}.
#'   Default \code{2000L}.
#' @param x_name,y_name Optional character. Variable names used in warning
#'   messages to identify which pair triggered a sparsity warning. If
#'   \code{NULL} (default), the names are deparsed from the call. Primarily
#'   used internally by \code{\link{build_graph}}; users do not normally
#'   supply these.
#'   
#' @return A named list with:
#' \describe{
#'   \item{statistic}{The chi-square test statistic.}
#'   \item{p_value}{The p-value of the test.}
#'   \item{df}{Degrees of freedom (NA under simulate_p = TRUE).}
#'   \item{n}{Number of pairwise-complete observations.}
#'   \item{table}{The contingency table.}
#'   \item{type}{Table type from \code{\link{detect_type}}.}
#' }
#'
#' @details
#' Rows with \code{NA} in either \code{x} or \code{y} are removed before
#' tabulation (pairwise deletion). Two warning conditions are checked:
#'
#' \enumerate{
#'   \item \strong{Any expected count < 5} (Cochran, 1954): a minor sparsity
#'         warning. Consider \code{simulate_p = TRUE}.
#'   \item \strong{Severe sparsity} - either more than 20\% of cells have
#'         expected count < 5, or in a table of minimum dimension \eqn{\geq 3}
#'         the observations-to-cells ratio is below 5. In this regime the
#'         chi-square approximation is unreliable \emph{and} Cramer's V becomes
#'         numerically unstable, especially in its classical (uncorrected)
#'         form. Consider collapsing categories, enabling bias correction
#'         (\code{corrected = TRUE} in \code{\link{effect_size}}), or using
#'         \code{simulate_p = TRUE}.
#' }
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
#'   \doi{10.1002/0471249688}
#'
#' Cochran, W. G. (1954). Some methods for strengthening the common
#'   chi-square tests. \emph{Biometrics}, 10(4), 417--451.
#'   \doi{10.2307/3001616}
#'
#' @examples
#' set.seed(42)
#' x <- sample(c("A", "B"), 100, replace = TRUE)
#' y <- sample(c("yes", "no"), 100, replace = TRUE)
#' result <- compute_assoc(x, y)
#' result$statistic
#' result$p_value
#'
#' @seealso \code{\link{effect_size}}, \code{\link{detect_type}}
#' @importFrom stats chisq.test
#' @export
compute_assoc <- function(x, y,
                          correct    = FALSE,
                          simulate_p = FALSE,
                          B          = 2000L,
                          x_name     = NULL,
                          y_name     = NULL) {
  if (length(x) != length(y)) {
    stop("`x` and `y` must have the same length.", call. = FALSE)
  }
  
  allowed_x <- is.factor(x) || is.character(x) || is.logical(x)
  allowed_y <- is.factor(y) || is.character(y) || is.logical(y)
  
  if (!allowed_x) {
    stop("`x` must be a factor, character, or logical vector.", call. = FALSE)
  }
  
  if (!allowed_y) {
    stop("`y` must be a factor, character, or logical vector.", call. = FALSE)
  }
  
  if (is.null(x_name)) {
    x_expr <- substitute(x)
    x_name <- if (is.symbol(x_expr)) deparse(x_expr) else "x"
  }
  
  if (is.null(y_name)) {
    y_expr <- substitute(y)
    y_name <- if (is.symbol(y_expr)) deparse(y_expr) else "y"
  }
  
  complete_idx <- !is.na(x) & !is.na(y)
  x_obs <- x[complete_idx]
  y_obs <- y[complete_idx]
  n_obs <- sum(complete_idx)
  
  if (n_obs < 2L) {
    warning("Fewer than 2 complete observations; returning NA results.", call. = FALSE)
    return(list(
      statistic = NA_real_,
      p_value   = NA_real_,
      df        = NA_real_,
      n         = n_obs,
      table     = NULL,
      type      = NA_character_
    ))
  }
  
  tab <- table(x = x_obs, y = y_obs, dnn = c("x", "y"))
  
  if (nrow(tab) < 2L || ncol(tab) < 2L) {
    warning(
      sprintf(
        "Contingency table for pair (%s, %s) has fewer than 2 levels in at least one variable after removing missing values; returning NA results.",
        x_name, y_name
      ),
      call. = FALSE
    )
    return(list(
      statistic = NA_real_,
      p_value   = NA_real_,
      df        = NA_real_,
      n         = n_obs,
      table     = tab,
      type      = "degenerate"
    ))
  }
  
  tab_type <- if (nrow(tab) == 2L && ncol(tab) == 2L) "2x2" else "RxC"
  
  expected   <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
  n_cells    <- nrow(tab) * ncol(tab)
  prop_small <- mean(expected < 5)
  any_small  <- any(expected < 5)
  min_dim    <- min(nrow(tab), ncol(tab))
  
  if (any_small) {
    warning(
      paste0(
        "At least one expected cell frequency is < 5 for pair (",
        x_name, ", ", y_name,
        "). Consider setting simulate_p = TRUE."
      ),
      call. = FALSE
    )
  }
  
  severe <- prop_small > 0.20 || (min_dim >= 3L && n_obs / n_cells < 5)
  
  if (severe) {
    warning(
      sprintf(
        paste0(
          "Sparse contingency table for pair (%s, %s): %dx%d = %d cells, ",
          "%d obs, %.0f%% cells with E < 5. Cramer's V and chi-square p-values ",
          "may be unstable; consider collapsing categories, enabling bias correction, ",
          "or simulate_p = TRUE."
        ),
        x_name, y_name, nrow(tab), ncol(tab), n_cells, n_obs, 100 * prop_small
      ),
      call. = FALSE
    )
  }
  
  test <- suppressWarnings(
    chisq.test(
      tab,
      correct = correct,
      simulate.p.value = simulate_p,
      B = B
    )
  )
  
  list(
    statistic = unname(test$statistic),
    p_value   = test$p.value,
    df        = if (simulate_p) NA_real_ else unname(test$parameter),
    n         = n_obs,
    table     = tab,
    type      = tab_type
  )
}