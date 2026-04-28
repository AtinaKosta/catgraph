#' Bootstrap confidence intervals for phi or Cramer's V
#'
#' Estimates a confidence interval for a pairwise effect size (phi or
#' Cramer's V) using a non-parametric percentile bootstrap. Rows are
#' resampled with replacement from the pairwise-complete observations,
#' and the chosen effect size estimator is recomputed on each resample.
#' The resulting empirical distribution is used to derive the interval.
#'
#' @param x A factor, character, or logical vector.
#' @param y A factor, character, or logical vector of the same length as
#'   \code{x}.
#' @param R Integer. Number of bootstrap resamples. Default \code{1000L}.
#'   Values below 500 are accepted with a warning.
#' @param conf Numeric in (0, 1). Confidence level. Default \code{0.95}.
#' @param type Character. Bootstrap interval type: \code{"percentile"}
#'   (default) or \code{"bca"} (bias-corrected and accelerated, Efron 1987).
#'   The BCa interval generally has better coverage but requires the
#'   jackknife influence values and is slower.
#' @param corrected Logical. Whether to use the bias-corrected estimator
#'   (Bergsma, 2013). Default \code{FALSE}.
#' @param correct Logical. Yates' continuity correction. Default \code{FALSE}.
#' @param seed Integer or \code{NULL}. Optional random seed for
#'   reproducibility. Default \code{NULL} (no seed set).
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{estimate}}{The point estimate of phi or Cramer's V on the
#'     original data.}
#'   \item{\code{ci_lower}}{Lower confidence bound.}
#'   \item{\code{ci_upper}}{Upper confidence bound.}
#'   \item{\code{conf}}{The requested confidence level.}
#'   \item{\code{type}}{The interval type used.}
#'   \item{\code{R}}{Number of resamples actually used (may be lower than
#'     requested if degenerate resamples were removed).}
#'   \item{\code{metric}}{Character: \code{"phi"} or \code{"cramers_v"}.}
#'   \item{\code{corrected}}{Logical.}
#'   \item{\code{n}}{Number of pairwise-complete observations.}
#'   \item{\code{boot_dist}}{Numeric vector of length \code{R} containing
#'     the full bootstrap distribution (useful for plotting).}
#' }
#'
#' @details
#' \strong{Percentile interval} (Efron & Tibshirani, 1993, Ch. 13):
#'
#' The \eqn{\alpha/2} and \eqn{1-\alpha/2} quantiles of the bootstrap
#' distribution \eqn{\hat{\theta}^*_1, \ldots, \hat{\theta}^*_R} are used
#' directly as the lower and upper bounds.
#'
#' \strong{BCa interval} (Efron, 1987):
#'
#' Adjusts the quantiles used for the interval by a bias-correction term
#' \eqn{\hat{z}_0} (estimated from the proportion of bootstrap resamples
#' below the point estimate) and an acceleration term \eqn{\hat{a}}
#' (estimated from the jackknife influence values). This gives second-order
#' accurate intervals without requiring a transformation. The formulas
#' follow DiCiccio & Efron (1996).
#'
#' \strong{Degenerate resamples}: Bootstrap resamples that yield a
#' single-level variable (all observations in one category) produce
#' \code{NA} effect sizes and are discarded. A warning is issued if more
#' than 5\% of resamples are discarded.
#'
#' \strong{Effect size on a boundary}: Both phi and Cramer's V are bounded
#' at 0. Bootstrap distributions for weak associations are therefore
#' right-skewed and pile up at 0. The BCa correction partially addresses
#' this; for near-zero associations the percentile lower bound may be 0,
#' which is correct.
#'
#' @references
#' DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
#'   \emph{Statistical Science}, 11(3), 189--212.
#'   \doi{10.1214/ss/1032280214}
#'
#' Efron, B. (1987). Better bootstrap confidence intervals.
#'   \emph{Journal of the American Statistical Association}, 82(397), 171--185.
#'   \doi{10.1080/01621459.1987.10478410}
#'
#' Efron, B., & Tibshirani, R. J. (1993).
#'   \emph{An Introduction to the Bootstrap}. Chapman & Hall.
#'   \doi{10.1007/978-1-4899-4541-9}
#'
#' @examples
#' set.seed(42)
#' x <- sample(c("A", "B"), 120, replace = TRUE,
#'             prob = c(0.4, 0.6))
#' y <- ifelse(x == "A",
#'             sample(c("yes", "no"), 120, replace = TRUE, prob = c(0.7, 0.3)),
#'             sample(c("yes", "no"), 120, replace = TRUE, prob = c(0.3, 0.7)))
#'
#' # Small R for a fast example; use R >= 1000 in real work.
#' ci <- bootstrap_ci(x, y, R = 200, seed = 1)
#' ci$estimate
#' ci$ci_lower
#' ci$ci_upper
#'
#' \donttest{
#' # BCa interval (slower: adds a leave-one-out jackknife)
#' ci_bca <- bootstrap_ci(x, y, R = 500, type = "bca", seed = 1)
#' }
#'
#' @seealso \code{\link{effect_size}}, \code{\link{catgraph_ci}}
#' @importFrom stats quantile qnorm pnorm
#' @export
bootstrap_ci <- function(x, y,
                         R         = 1000L,
                         conf      = 0.95,
                         type      = c("percentile", "bca"),
                         corrected = FALSE,
                         correct   = FALSE,
                         seed      = NULL) {

  type <- match.arg(type)

  if (R < 500L) {
    warning("R < 500 may yield unstable confidence intervals.")
  }
  if (conf <= 0 || conf >= 1) {
    stop("`conf` must be strictly between 0 and 1.")
  }
  if (!is.null(seed)) set.seed(seed)

  stopifnot(length(x) == length(y))

  # Pairwise deletion on original data
  ok  <- !is.na(x) & !is.na(y)
  xc  <- x[ok]
  yc  <- y[ok]
  n   <- sum(ok)

  if (n < 4L) {
    stop("Fewer than 4 pairwise-complete observations; cannot bootstrap.")
  }

  # Point estimate on original data
  orig <- effect_size(xc, yc, corrected = corrected, correct = correct)
  est  <- orig$effect_size

  # --------------------------------------------------- bootstrap loop
  boot_vals <- numeric(R)
  n_degen   <- 0L

  for (i in seq_len(R)) {
    idx  <- sample.int(n, n, replace = TRUE)
    xb   <- xc[idx]
    yb   <- yc[idx]

    # Skip degenerate resamples (single-level variable)
    if (length(unique(xb)) < 2L || length(unique(yb)) < 2L) {
      boot_vals[i] <- NA_real_
      n_degen      <- n_degen + 1L
      next
    }

    res <- tryCatch(
      effect_size(xb, yb, corrected = corrected, correct = correct),
      error   = function(e) list(effect_size = NA_real_),
      warning = function(w) {
        suppressWarnings(
          effect_size(xb, yb, corrected = corrected, correct = correct)
        )
      }
    )
    boot_vals[i] <- if (is.na(res$effect_size)) NA_real_ else res$effect_size
  }

  # Remove NAs from degenerate resamples
  valid <- boot_vals[!is.na(boot_vals)]
  R_eff <- length(valid)

  if (n_degen / R > 0.05) {
    warning(sprintf(
      "%.1f%% of bootstrap resamples were degenerate (single-level variable) ",
      100 * n_degen / R
    ), "and were discarded. Interval may be unreliable.")
  }

  if (R_eff < 10L) {
    stop("Fewer than 10 valid bootstrap resamples; cannot compute interval.")
  }

  alpha <- 1 - conf

  # ------------------------------------------------ interval calculation
  if (type == "percentile") {
    ci_lo <- stats::quantile(valid, probs = alpha / 2,     names = FALSE)
    ci_hi <- stats::quantile(valid, probs = 1 - alpha / 2, names = FALSE)

  } else {
    # --------------------------------------------------------- BCa interval
    # Bias-correction z0: proportion of boot values < original estimate
    z0 <- stats::qnorm(mean(valid < est))

    # Acceleration a: jackknife influence values
    jack_vals <- numeric(n)
    for (j in seq_len(n)) {
      xj  <- xc[-j]
      yj  <- yc[-j]
      if (length(unique(xj)) < 2L || length(unique(yj)) < 2L) {
        jack_vals[j] <- 0
        next
      }
      rj <- tryCatch(
        effect_size(xj, yj, corrected = corrected, correct = correct),
        error   = function(e) list(effect_size = est),
        warning = function(w) suppressWarnings(
          effect_size(xj, yj, corrected = corrected, correct = correct)
        )
      )
      jack_vals[j] <- if (is.na(rj$effect_size)) est else rj$effect_size
    }

    jack_mean <- mean(jack_vals)
    num  <- sum((jack_mean - jack_vals)^3)
    denom <- 6 * (sum((jack_mean - jack_vals)^2))^(3/2)
    a <- if (abs(denom) < .Machine$double.eps) 0 else num / denom

    # Adjusted quantiles
    z_alpha_lo <- stats::qnorm(alpha / 2)
    z_alpha_hi <- stats::qnorm(1 - alpha / 2)

    p_lo <- stats::pnorm(z0 + (z0 + z_alpha_lo) / (1 - a * (z0 + z_alpha_lo)))
    p_hi <- stats::pnorm(z0 + (z0 + z_alpha_hi) / (1 - a * (z0 + z_alpha_hi)))

    # Clamp probabilities to valid range
    p_lo <- max(0, min(1, p_lo))
    p_hi <- max(0, min(1, p_hi))

    ci_lo <- stats::quantile(valid, probs = p_lo, names = FALSE)
    ci_hi <- stats::quantile(valid, probs = p_hi, names = FALSE)
  }

  list(
    estimate  = est,
    ci_lower  = ci_lo,
    ci_upper  = ci_hi,
    conf      = conf,
    type      = type,
    R         = R_eff,
    metric    = orig$metric,
    corrected = corrected,
    n         = n,
    boot_dist = valid
  )
}


#' Add bootstrap confidence intervals to all edges of a catgraph
#'
#' Calls \code{\link{bootstrap_ci}} for every edge in a \code{catgraph}
#' object and stores the lower and upper bounds as additional edge
#' attributes (\code{ci_lower}, \code{ci_upper}, \code{ci_conf},
#' \code{ci_type}).
#'
#' @param x A \code{catgraph} object.
#' @param R Integer. Number of bootstrap resamples per pair. Default
#'   \code{1000L}.
#' @param conf Numeric. Confidence level. Default \code{0.95}.
#' @param type Character. \code{"percentile"} or \code{"bca"}.
#'   Default \code{"percentile"}.
#' @param seed Integer or \code{NULL}. Base seed; each pair uses
#'   \code{seed + i} to ensure per-pair reproducibility without global
#'   state. Default \code{NULL}.
#' @param verbose Logical. Print a progress counter. Default \code{TRUE}.
#'
#' @return The input \code{catgraph} object with three new edge attributes:
#'   \code{ci_lower}, \code{ci_upper}, \code{ci_conf}, \code{ci_type}.
#'
#' @details
#' For large graphs (many variable pairs) this function can be slow because
#' it runs \code{R} resamples per edge. Consider lowering \code{R} or
#' running on a pruned graph (via \code{\link{prune_edges}}) to reduce
#' computation.
#'
#' @examples
#' df <- as.data.frame(Titanic)
#' df_exp <- df[rep(seq_len(nrow(df)), df$Freq), -5]
#' cg <- catgraph(df_exp)
#' \donttest{
#' cg <- catgraph_ci(cg, R = 500, seed = 1)
#' igraph::E(cg$graph)$ci_lower
#' igraph::E(cg$graph)$ci_upper
#' }
#'
#' @seealso \code{\link{bootstrap_ci}}, \code{\link{catgraph}}
#' @importFrom igraph as_edgelist E
#' @importFrom utils flush.console
#' @export
catgraph_ci <- function(x,
                        R       = 1000L,
                        conf    = 0.95,
                        type    = c("percentile", "bca"),
                        seed    = NULL,
                        verbose = TRUE) {

  if (!inherits(x, "catgraph")) stop("`x` must be a catgraph object.")
  type <- match.arg(type)

  g  <- x$graph
  el <- igraph::as_edgelist(g, names = TRUE)
  ne <- nrow(el)

  ci_lo   <- numeric(ne)
  ci_hi   <- numeric(ne)

  for (i in seq_len(ne)) {
    if (verbose) {
      cat(sprintf("\r  Bootstrap CI: edge %d / %d", i, ne))
      utils::flush.console()
    }

    v1 <- el[i, 1]
    v2 <- el[i, 2]

    pair_seed <- if (!is.null(seed)) seed + i else NULL

    ci <- tryCatch(
      bootstrap_ci(
        x$data[[v1]], x$data[[v2]],
        R         = R,
        conf      = conf,
        type      = type,
        corrected = x$corrected,
        seed      = pair_seed
      ),
      error = function(e) list(ci_lower = NA_real_, ci_upper = NA_real_)
    )

    ci_lo[i] <- ci$ci_lower
    ci_hi[i] <- ci$ci_upper
  }

  if (verbose) cat("\n")

  igraph::E(g)$ci_lower <- ci_lo
  igraph::E(g)$ci_upper <- ci_hi
  igraph::E(g)$ci_conf  <- conf
  igraph::E(g)$ci_type  <- type

  x$graph <- g
  x
}
