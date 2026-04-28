#' Expand a contingency table or frequency data frame to observation-level format
#'
#' Converts multi-dimensional contingency tables (\code{table}, \code{array},
#' \code{ftable}), or data frames that contain a frequency/count column, into
#' a flat data frame where every row represents one observation. This is the
#' required input format for \code{\link{catgraph}}.
#'
#' @param tbl A \code{table}, \code{array}, \code{ftable}, or \code{data.frame}.
#'   All standard R contingency table formats are supported, including
#'   multi-dimensional arrays such as \code{Titanic} (4-D) and
#'   \code{HairEyeColor} (3-D).
#' @param freq_col Character or integer. Only used when \code{tbl} is a
#'   \code{data.frame}. The name or column index of the frequency/count column.
#'   If \code{NULL} (default), the function looks for a column named
#'   \code{"Freq"}, \code{"freq"}, \code{"n"}, \code{"count"}, or
#'   \code{"Count"} in that order. An error is raised if none is found.
#' @param as_factor Logical. If \code{TRUE} (default), all resulting columns
#'   are coerced to factors, preserving the level order from the original
#'   object's \code{dimnames}. Set to \code{FALSE} to return character columns.
#' @param drop_zero Logical. If \code{TRUE} (default), rows corresponding to
#'   zero-count cells are silently dropped. Set to \code{FALSE} to keep them
#'   (those rows will appear zero times and therefore not affect results, but
#'   the factor levels will still be present).
#'
#' @return A \code{data.frame} with one row per observation and one column per
#'   categorical variable. The column names are taken from \code{dimnames()}
#'   of the input object, or from the non-frequency columns of the input
#'   data frame. Row names are reset to \code{NULL}.
#'
#' @details
#' \strong{Accepted input formats:}
#' \describe{
#'   \item{One row per observation}{Already in the correct format — pass
#'     directly to \code{\link{catgraph}} without calling
#'     \code{expand_table()}.}
#'   \item{\code{table} / \code{array}}{The standard output of
#'     \code{table()}, \code{xtabs()}, or built-in datasets such as
#'     \code{Titanic}, \code{HairEyeColor}, and \code{UCBAdmissions}.}
#'   \item{\code{ftable}}{Converted to a \code{table} internally before
#'     expansion.}
#'   \item{\code{data.frame} with frequency column}{The output of
#'     \code{as.data.frame(some_table)}, which always contains a
#'     \code{Freq} column.}
#' }
#'
#' \strong{Not accepted:}
#' \describe{
#'   \item{Raw numeric matrices}{A plain matrix of counts without
#'     \code{dimnames} cannot be safely converted. Assign \code{dimnames}
#'     first.}
#'   \item{Numeric 0/1 columns}{Columns coded as integers or doubles are
#'     not treated as categorical. Coerce with \code{as.factor()} or
#'     \code{as.character()} first, or pass them through
#'     \code{\link{catgraph}} which will coerce and warn automatically.}
#' }
#'
#' @references
#' R Core Team (2024). \emph{R: A Language and Environment for Statistical
#'   Computing}. R Foundation for Statistical Computing, Vienna, Austria.
#'   \url{https://www.R-project.org/}
#'
#' @examples
#' # Built-in 4-D table
#' df <- expand_table(Titanic)
#' str(df)
#' nrow(df)  # 2201 passengers
#'
#' # Built-in 3-D table
#' df2 <- expand_table(HairEyeColor)
#' str(df2)
#'
#' # data.frame with Freq column (output of as.data.frame on a table)
#' tab_df <- as.data.frame(UCBAdmissions)
#' df3 <- expand_table(tab_df)
#' nrow(df3)  # 4526 applicants
#'
#' # Custom data frame with a count column
#' survey <- data.frame(
#'   gender   = c("M", "F", "M", "F"),
#'   smokes   = c("yes", "yes", "no", "no"),
#'   n        = c(23L, 15L, 48L, 61L)
#' )
#' df4 <- expand_table(survey, freq_col = "n")
#' nrow(df4)  # 147 observations
#'
#' # Use directly with catgraph
#' cg <- catgraph(expand_table(Titanic))
#' cg
#'
#' @seealso \code{\link{catgraph}}, \code{\link{assoc_matrix}}
#' @importFrom stats ftable
#' @export
expand_table <- function(tbl,
                         freq_col   = NULL,
                         as_factor  = TRUE,
                         drop_zero  = TRUE) {

  # ----------------------------------------------------------------- ftable
  if (inherits(tbl, "ftable")) {
    tbl <- stats::ftable(tbl)   # re-cast to get standard structure
    tbl <- as.table(tbl)
  }

  # ------------------------------------------------------------------ table / array
  if (is.table(tbl) || is.array(tbl)) {
    df <- as.data.frame(tbl, stringsAsFactors = FALSE)
    # as.data.frame on a table always produces a "Freq" column
    freq_col <- "Freq"
  } else if (is.data.frame(tbl)) {
    df <- tbl
    # Auto-detect frequency column if not supplied
    if (is.null(freq_col)) {
      candidates <- c("Freq", "freq", "n", "count", "Count")
      found <- intersect(candidates, names(df))
      if (length(found) == 0L) {
        stop(
          "Could not find a frequency column in the data frame. ",
          "Looked for: ", paste(candidates, collapse = ", "), ". ",
          "Please supply the column name via `freq_col`."
        )
      }
      freq_col <- found[1L]
    } else if (is.numeric(freq_col)) {
      freq_col <- names(df)[freq_col]
    }
  } else {
    stop(
      "`tbl` must be a table, array, ftable, or data.frame. ",
      "Got: ", class(tbl)[1L], "."
    )
  }

  # Validate frequency column
  if (!freq_col %in% names(df)) {
    stop("Frequency column '", freq_col, "' not found in data.")
  }
  freqs <- df[[freq_col]]
  if (!is.numeric(freqs)) {
    stop("Frequency column '", freq_col, "' must be numeric.")
  }
  if (any(freqs < 0, na.rm = TRUE)) {
    stop("Frequency column '", freq_col, "' contains negative values.")
  }
  freqs <- as.integer(round(freqs))

  # Drop zero-count rows
  if (drop_zero) {
    keep  <- !is.na(freqs) & freqs > 0L
    df    <- df[keep, , drop = FALSE]
    freqs <- freqs[keep]
  }

  if (sum(freqs, na.rm = TRUE) == 0L) {
    stop("All frequency counts are zero or NA after dropping. Nothing to expand.")
  }

  # Identify categorical columns (everything except the freq column)
  cat_cols <- setdiff(names(df), freq_col)

  if (length(cat_cols) < 2L) {
    stop(
      "After removing the frequency column, fewer than 2 categorical columns ",
      "remain. catgraph() requires at least 2 variables."
    )
  }

  # Expand: replicate each row according to its frequency
  idx    <- rep(seq_len(nrow(df)), times = freqs)
  result <- df[idx, cat_cols, drop = FALSE]
  rownames(result) <- NULL

  # Coerce to factor if requested
  if (as_factor) {
    result[] <- lapply(result, function(col) {
      if (is.factor(col)) col else factor(col)
    })
  }

  result
}
