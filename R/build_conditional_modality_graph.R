#' Build a modality graph for an observed subgroup
#'
#' Subsets \code{data} by one or more observed conditioning variables and
#' levels, then calls \code{\link{build_modality_graph}} on the remaining
#' variables. The conditioning specification is stored on the returned object
#' so that downstream analyses can report the subgroup definition.
#'
#' This is the subgroup-comparison entry point: for example, "what does the
#' modality association structure look like among women only?" or "among
#' respondents from wave 3 only?". The function conditions only on observed
#' variables supplied by the user. It does not estimate latent classes, infer
#' clusters of respondents, or perform causal adjustment.
#'
#' @param data A data frame of categorical variables.
#' @param given A named list specifying the conditioning. Each name is a
#'   column in \code{data}; each value is a length-1 character vector
#'   giving the level to condition on. Example:
#'   \code{list(sex = "female", wave = "2020")}. The conditioning
#'   columns are dropped from the modality graph by construction, so
#'   they do not appear as nodes.
#' @param drop_conditioning_vars Logical. If \code{TRUE} (default), the
#'   conditioning columns are removed before graph construction. If
#'   \code{FALSE}, they are retained, in which case they will appear
#'   as single-modality variables and contribute no edges (their only
#'   remaining level has probability 1 on the conditional subset).
#'   The default is almost always what you want; \code{FALSE} is
#'   provided for debugging and diagnostic use.
#' @param ... Further arguments passed to \code{\link{build_modality_graph}}
#'   (e.g., \code{remove_na}, \code{min_count}).
#'
#' @return A \code{catmodgraph} object as returned by
#'   \code{\link{build_modality_graph}}, augmented with a
#'   \code{conditioning} component: a list with elements
#'   \code{given} (the input specification), \code{n_original} (row
#'   count of the input data), and \code{n_conditional} (row count
#'   after subsetting). The S3 class is unchanged so that existing
#'   methods (plotting, pruning, clustering, comparison) work without
#'   modification.
#'
#' @details
#' \strong{Conditioning variables must be observed, not estimated.}
#' This function is designed for sociodemographic strata (sex,
#' country, wave, study site) --- variables whose value is recorded
#' directly in the data. Conditioning on a derived cluster label
#' reintroduces the circularity that the 0.6.0 refocus removed: the
#' cluster was itself estimated from the association structure, so
#' the "conditional" graph is not meaningfully separable from the
#' overall graph. If you need to analyse the association structure
#' within latent classes, use \pkg{poLCA} to fit the class model
#' first, then pass the external class assignment as observed input.
#'
#' \strong{Composition with other functions.} The returned object is
#' a plain \code{catmodgraph} and can be passed directly to
#' \code{\link{plot.catmodgraph}}, \code{\link{prune_modality_edges}},
#' \code{\link{cluster_modalities}}, and
#' \code{\link{compare_modality_graphs}}. When used in
#' \code{compare_modality_graphs()}, supply a named list where the
#' names are short descriptions of the conditioning (e.g.,
#' \code{list(women = mg_f, men = mg_m)}); the function does not
#' read the \code{conditioning} attribute to generate panel titles.
#'
#' @examples
#' \donttest{
#' data(survey_health)
#' mg_f <- build_conditional_modality_graph(
#'   survey_health, given = list(sex = "female")
#' )
#' mg_m <- build_conditional_modality_graph(
#'   survey_health, given = list(sex = "male")
#' )
#'
#' compare_modality_graphs(list(women = mg_f, men = mg_m),
#'                         restrict = "common")
#' }
#'
#' @seealso \code{\link{build_modality_graph}},
#'   \code{\link{compare_modality_graphs}},
#'   \code{\link{joint_balance}}
#' @export
build_conditional_modality_graph <- function(data,
                                             given,
                                             drop_conditioning_vars = TRUE,
                                             ...) {
  
  # ---- Validation -----------------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  
  if (!is.list(given) || length(given) == 0L ||
      is.null(names(given)) || any(names(given) == "")) {
    stop("`given` must be a non-empty named list, e.g. ",
         "list(sex = \"female\").", call. = FALSE)
  }
  if (anyDuplicated(names(given))) {
    stop("`given` must not contain duplicate names.", call. = FALSE)
  }
  
  missing_cols <- setdiff(names(given), names(data))
  if (length(missing_cols) > 0L) {
    stop("Conditioning columns not found in `data`: ",
         paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  
  bad_len <- vapply(given, function(v) length(v) != 1L, logical(1L))
  if (any(bad_len)) {
    stop("Each element of `given` must be a length-1 value. ",
         "Problem columns: ",
         paste(names(given)[bad_len], collapse = ", "),
         ". If you want to condition on multiple levels of a single ",
         "variable, pre-filter `data` and call build_modality_graph() ",
         "directly.",
         call. = FALSE)
  }
  
  if (!is.logical(drop_conditioning_vars) ||
      length(drop_conditioning_vars) != 1L ||
      is.na(drop_conditioning_vars)) {
    stop("`drop_conditioning_vars` must be TRUE or FALSE.", call. = FALSE)
  }
  
  # ---- Check each conditioning level exists --------------------------------
  for (col in names(given)) {
    col_vals <- unique(as.character(data[[col]]))
    if (!as.character(given[[col]]) %in% col_vals) {
      stop("Level '", given[[col]], "' not found in column '", col,
           "'. Observed levels: ",
           paste(utils::head(col_vals, 10), collapse = ", "),
           if (length(col_vals) > 10L) ", ..." else "",
           ".", call. = FALSE)
    }
  }
  
  # ---- Apply the conditioning ----------------------------------------------
  n_original <- nrow(data)
  mask <- rep(TRUE, n_original)
  for (col in names(given)) {
    mask <- mask & (!is.na(data[[col]]) &
                      as.character(data[[col]]) == as.character(given[[col]]))
  }
  sub_data <- data[mask, , drop = FALSE]
  n_cond   <- nrow(sub_data)
  
  if (n_cond == 0L) {
    stop("Conditioning yields an empty subset: no rows satisfy ",
         .format_given(given), ".", call. = FALSE)
  }
  if (n_cond < 10L) {
    warning("Conditional subset has only ", n_cond, " rows; ",
            "modality-graph estimates will be unstable.",
            call. = FALSE)
  }
  
  # ---- Drop conditioning columns (default) ---------------------------------
  if (drop_conditioning_vars) {
    sub_data <- sub_data[, setdiff(names(sub_data), names(given)),
                         drop = FALSE]
  }
  
  if (ncol(sub_data) < 2L) {
    stop("Fewer than 2 non-conditioning variables remain after ",
         "subsetting. Conditional modality graphs require at least ",
         "2 free variables.", call. = FALSE)
  }
  
  # ---- Build the graph and attach conditioning metadata --------------------
  mg <- build_modality_graph(sub_data, ...)
  
  mg$conditioning <- list(
    given         = given,
    n_original    = n_original,
    n_conditional = n_cond,
    dropped_vars  = if (drop_conditioning_vars) names(given) else character(0)
  )
  
  mg
}


# Internal helper: format a `given` list for error/warning messages
.format_given <- function(given) {
  parts <- paste0(names(given), " = \"",
                  vapply(given, as.character, character(1L)), "\"")
  paste(parts, collapse = ", ")
}