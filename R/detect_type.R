#' Detect contingency table type for a pair of categorical variables
#'
#' For a given pair of columns, this function determines whether the
#' cross-tabulation yields a 2x2 or a larger (n x m) contingency table,
#' after removing rows with missing values in either column (pairwise
#' deletion). This drives the dispatch between the phi coefficient and
#' Cramer's V in \code{\link{effect_size}}.
#'
#' @param x A factor, character, or logical vector.
#' @param y A factor, character, or logical vector of the same length as
#'   \code{x}.
#'
#' @return A character string: \code{"2x2"} if both variables have exactly
#'   two levels observed in the pairwise-complete data, or \code{"RxC"}
#'   otherwise.
#'
#' @details
#' Level counts are computed on the pairwise-complete subset, not on the
#' full vector. This means that a variable with three levels in the full
#' data may be classified as having two levels for a specific pair if one
#' level is entirely missing when the other variable is observed.
#'
#' @examples
#' x <- c("a", "b", "a", "b", NA)
#' y <- c("yes", "no", "yes", "yes", "no")
#' detect_type(x, y)  # "2x2"
#'
#' z <- c("low", "mid", "high", "low", "mid")
#' detect_type(x[-5], z[-5])  # "RxC"
#'
#' @seealso \code{\link{effect_size}}, \code{\link{compute_assoc}}
#' @export
detect_type <- function(x, y) {
  stopifnot(length(x) == length(y))

  # Pairwise deletion: keep only rows where both are non-missing
  complete_idx <- !is.na(x) & !is.na(y)
  x_obs <- x[complete_idx]
  y_obs <- y[complete_idx]

  if (sum(complete_idx) < 2L) {
    warning("Fewer than 2 complete observations for this pair; returning 'RxC'.")
    return("RxC")
  }

  # Count distinct observed levels (drop unused factor levels)
  nx <- length(unique(x_obs))
  ny <- length(unique(y_obs))

  if (nx == 2L && ny == 2L) "2x2" else "RxC"
}
