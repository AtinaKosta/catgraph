#' @describeIn joint_balance Print a concise diagnostic summary.
#' @param x A \code{jointbalance} object.
#' @param ... Ignored.
#' @export
print.jointbalance <- function(x, ...) {
  cat("jointbalance diagnostic\n")
  cat("  Group variable :", x$group, "\n")
  cat("  Levels         :", paste(x$group_levels, collapse = ", "), "\n")
  cat("  Variables      :", length(x$variables), "(",
      paste(utils::head(x$variables, 5), collapse = ", "),
      if (length(x$variables) > 5L) ", ..." else "", ")\n")
  cat("  Alpha          :", x$alpha, "\n\n")
  
  n_marg_sig <- sum(x$marginal$p_adjusted < x$alpha, na.rm = TRUE)
  cat("Marginal tests (BH-adjusted across variables):\n")
  cat("  Variables rejecting at alpha :",
      n_marg_sig, "/", nrow(x$marginal), "\n\n")
  
  n_pair_sig <- sum(x$pairwise_omnibus$p_bonferroni < x$alpha, na.rm = TRUE)
  cat("Pairwise omnibus tests (Bonferroni across group pairs):\n")
  cat("  Pairs rejecting at alpha     :",
      n_pair_sig, "/", nrow(x$pairwise_omnibus), "\n")
  if (n_pair_sig > 0L) {
    rej <- x$pairwise_omnibus[x$pairwise_omnibus$p_bonferroni < x$alpha, ]
    rej <- rej[order(rej$p_bonferroni), ]
    cat("  Rejecting pairs:\n")
    for (i in seq_len(nrow(rej))) {
      cat(sprintf("    %s vs %s : p_bonf = %.4g\n",
                  rej$group_x[i], rej$group_y[i], rej$p_bonferroni[i]))
    }
  }
  cat("\n")
  
  if (length(x$pairwise_edgewise) > 0L) {
    cat("Edge-wise post-hoc (BH across edges, within each pair):\n")
    for (nm in names(x$pairwise_edgewise)) {
      et <- x$pairwise_edgewise[[nm]]
      n_sig <- sum(et$edge_table$p_adjusted < x$alpha, na.rm = TRUE)
      cat(sprintf("  %-30s : %d of %d edges sig.\n",
                  nm, n_sig, nrow(et$edge_table)))
    }
    cat("\n")
  }
  
  cat("Three separate testing families; see ?joint_balance for details.\n")
  invisible(x)
}

#' @describeIn joint_balance Return the key diagnostic tables.
#' @param object A \code{jointbalance} object.
#' @export
summary.jointbalance <- function(object, ...) {
  list(
    marginal_table         = object$marginal,
    pairwise_omnibus_table = object$pairwise_omnibus,
    n_edgewise_pairs       = length(object$pairwise_edgewise),
    group                  = object$group,
    group_levels           = object$group_levels,
    alpha                  = object$alpha
  )
}

#' Plot a jointbalance diagnostic
#'
#' Two-panel layout: the left panel is a bar chart of adjusted marginal
#' p-values per variable (the "Table 1" view); the right panel is the
#' modality-difference graph for one selected group pair (from
#' \code{\link{plot_modality_difference}}). For \eqn{k > 2} groups the
#' default pair is the one with the smallest Bonferroni-adjusted
#' omnibus p-value.
#'
#' @param x A \code{jointbalance} object.
#' @param pair Optional length-2 character vector specifying which
#'   group pair to visualise. Default \code{NULL} picks the most
#'   significant pair.
#' @param layout Character. One of \code{"side_by_side"} (default) or
#'   \code{"marginal_only"}. Use \code{"marginal_only"} when no
#'   pairwise omnibus rejected and there is no edge-wise result to
#'   plot.
#' @param ... Further arguments passed to
#'   \code{\link{plot_modality_difference}}.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.jointbalance <- function(x,
                              pair   = NULL,
                              layout = c("side_by_side", "marginal_only"),
                              ...) {
  layout <- match.arg(layout)
  
  # -- Pick the pair to visualise --------------------------------------------
  if (layout == "side_by_side") {
    if (length(x$pairwise_edgewise) == 0L) {
      message("No pair rejected the omnibus at alpha = ", x$alpha,
              "; falling back to layout = 'marginal_only'.")
      layout <- "marginal_only"
    } else if (is.null(pair)) {
      # pick the most significant rejecting pair
      po <- x$pairwise_omnibus
      po <- po[po$p_bonferroni < x$alpha, ]
      po <- po[order(po$p_bonferroni), ]
      pair <- c(po$group_x[1L], po$group_y[1L])
    } else {
      if (!is.character(pair) || length(pair) != 2L) {
        stop("`pair` must be a length-2 character vector.", call. = FALSE)
      }
      if (!all(pair %in% x$group_levels)) {
        stop("Elements of `pair` must be in group levels: ",
             paste(x$group_levels, collapse = ", "), ".", call. = FALSE)
      }
    }
  }
  
  # -- Layout -----------------------------------------------------------------
  op <- if (layout == "side_by_side") {
    graphics::par(mfrow = c(1, 2), mar = c(4, 8, 3, 1))
  } else {
    graphics::par(mar = c(4, 8, 3, 1))
  }
  on.exit(graphics::par(op), add = TRUE)
  
  # -- Panel 1: marginal bar chart -------------------------------------------
  m <- x$marginal
  m <- m[order(m$p_adjusted), ]
  # -log10 with a floor for display
  neglog <- -log10(pmax(m$p_adjusted, 1e-10))
  cols   <- ifelse(m$p_adjusted < x$alpha, "#534AB7", "#D3D1C7")
  
  graphics::barplot(
    rev(neglog),
    horiz     = TRUE,
    names.arg = rev(m$variable),
    las       = 1,
    col       = rev(cols),
    border    = NA,
    cex.names = 0.75,
    xlab      = expression(-log[10](p[adj])),
    main      = "Marginal balance (per variable)"
  )
  graphics::abline(v = -log10(x$alpha), lty = 2, col = "grey30")
  
  # -- Panel 2: modality-difference graph ------------------------------------
  if (layout == "side_by_side") {
    pair_key <- paste0(pair[1L], "_vs_", pair[2L])
    # Be tolerant to argument order
    if (!pair_key %in% names(x$pairwise_edgewise)) {
      alt_key <- paste0(pair[2L], "_vs_", pair[1L])
      if (alt_key %in% names(x$pairwise_edgewise)) {
        pair_key <- alt_key
        pair <- rev(pair)
      } else {
        stop("No edge-wise result stored for pair ",
             pair[1L], " vs ", pair[2L],
             ". Available: ",
             paste(names(x$pairwise_edgewise), collapse = ", "), ".",
             call. = FALSE)
      }
    }
    et <- x$pairwise_edgewise[[pair_key]]
    
    plot_modality_difference(
      et,
      reference    = list(x$modality_graphs[[pair[1L]]],
                          x$modality_graphs[[pair[2L]]]),
      group_labels = c(paste0("stronger in ", pair[1L]),
                       paste0("stronger in ", pair[2L])),
      title        = sprintf("Joint differences: %s vs %s",
                             pair[1L], pair[2L]),
      ...
    )
  }
  
  invisible(x)
}