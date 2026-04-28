#' Plot a heatmap of pairwise association weights
#'
#' Produces a colour-coded heatmap of the dense all-pairs effect-size matrix
#' implied by a \code{catgraph} object. Heatmap fill values are computed from
#' the processed data stored in \code{x$data} via \code{\link{assoc_similarity}},
#' so true zero associations are shown as 0 rather than treated as absent
#' graph edges.
#'
#' @param x A \code{catgraph} object.
#' @param engine Character. \code{"ggplot2"} (default) or \code{"base"}.
#'   The \code{ggplot2} engine requires \pkg{ggplot2}; the \code{base} engine
#'   uses only \pkg{graphics} from base R.
#' @param show_values Logical. Whether to print effect-size values inside
#'   each cell. Default \code{TRUE}.
#' @param show_sig Logical. Whether to overlay significance stars
#'   (\code{***} p < 0.001, \code{**} p < 0.01, \code{*} p < 0.05,
#'   \code{.} p < 0.1) below each value. Default \code{FALSE}.
#' @param show_ci Logical. Whether to show bootstrapped confidence intervals
#'   as \code{[lo, hi]} text beneath each value. Requires that
#'   \code{\link{catgraph_ci}} has been called on \code{x}. Default
#'   \code{FALSE}.
#' @param palette Character vector of colours defining the gradient from
#'   low (weak association) to high (strong association). Default is a
#'   perceptually uniform purple ramp derived from the package colour system.
#'   Pass any vector of hex colours to override.
#' @param digits Integer. Number of decimal places for cell labels.
#'   Default \code{2L}.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param na_fill Character. Fill colour for cells that could not be
#'   computed (e.g. degenerate pairs). Default \code{"#D3D1C7"}
#'   (gray-100).
#' @param reorder Logical. Whether to reorder variables by hierarchical
#'   clustering of the effect-size matrix so that similar variables are
#'   adjacent. Default \code{TRUE}.
#'
#' @return For \code{engine = "ggplot2"}: a \code{ggplot} object (can be
#'   further customised with ggplot2 layers).
#'   For \code{engine = "base"}: \code{NULL}, invisibly, called for its
#'   side effect.
#'
#' @details
#' \strong{Colour palette}: the default palette is a five-stop sequence
#' from white (V = 0) through lilac to deep purple (V = 1), matching the
#' purple ramp used throughout the package. This choice avoids the
#' red/green palette that is problematic for colour-blind readers. Pass
#' \code{palette = c("#FFFFFF", "#5DCAA5", "#0F6E56")} for a teal ramp, for
#' example.
#'
#' \strong{Reordering}: when \code{reorder = TRUE}, the variables are
#' permuted by the first two components of an angular-order seriation of the
#' correlation matrix, following the \emph{corrplot} convention (Wei &
#' Simko, 2021). Because effect sizes are always non-negative, the
#' clustering uses \eqn{1 - V} as a dissimilarity measure, which groups
#' strongly associated variables together.
#'
#' @references
#' Wei, T., & Simko, V. (2021). \emph{R package corrplot: Visualization of
#'   a Correlation Matrix}. Version 0.92.
#'   \url{https://github.com/taiyun/corrplot}
#'
#' @examples
#' df <- as.data.frame(Titanic)
#' df_exp <- df[rep(seq_len(nrow(df)), df$Freq), -5]
#' cg <- catgraph(df_exp)
#' plot_heatmap(cg)
#' plot_heatmap(cg, show_sig = TRUE)
#' plot_heatmap(cg, engine = "base")
#'
#' @seealso \code{\link{assoc_matrix}}, \code{\link{catgraph_ci}},
#'   \code{\link{plot.catgraph}}
#' @importFrom igraph as_edgelist E V
#' @export
plot_heatmap <- function(x,
                         engine      = c("ggplot2", "base"),
                         show_values = TRUE,
                         show_sig    = FALSE,
                         show_ci     = FALSE,
                         palette     = c("#FFFFFF", "#EEEDFE", "#AFA9EC",
                                         "#7F77DD", "#534AB7", "#26215C"),
                         digits      = 2L,
                         title       = NULL,
                         na_fill     = "#D3D1C7",
                         reorder     = TRUE) {

  engine <- match.arg(engine)
  if (!inherits(x, "catgraph")) stop("`x` must be a catgraph object.")

  # ------------------------------------------------- build long data frame
  # Use the dense all-pairs similarity matrix for heatmap fill values.
  # This preserves true zero associations as 0 rather than treating them as
  # absent graph edges / NA cells.
  mat  <- assoc_similarity(
    x$data,
    corrected  = x$corrected,
    what       = "effect_size"
  )
  
  vars <- rownames(mat)
  p    <- length(vars)

  # Check CI availability
  has_ci <- !is.null(igraph::E(x$graph)$ci_lower) && show_ci
  if (show_ci && !has_ci) {
    warning(
      "No CI attributes found. Run catgraph_ci() first. Ignoring show_ci."
    )
    show_ci <- FALSE
  }

  # Significance stars helper
  sig_stars <- function(pv) {
    ifelse(is.na(pv),   "",
    ifelse(pv < 0.001, "***",
    ifelse(pv < 0.01,  "**",
    ifelse(pv < 0.05,  "*",
    ifelse(pv < 0.1,   ".", "")))))
  }

  # Build pval and CI matrices
  pmat   <- matrix(NA_real_,     nrow = p, ncol = p, dimnames = list(vars, vars))
  lo_mat <- matrix(NA_real_,     nrow = p, ncol = p, dimnames = list(vars, vars))
  hi_mat <- matrix(NA_real_,     nrow = p, ncol = p, dimnames = list(vars, vars))

  el <- igraph::as_edgelist(x$graph, names = TRUE)
  pv <- igraph::E(x$graph)$p_value

  for (i in seq_len(nrow(el))) {
    a <- el[i, 1]; b <- el[i, 2]
    pmat[a, b] <- pmat[b, a] <- pv[i]
  }

  if (has_ci) {
    lo <- igraph::E(x$graph)$ci_lower
    hi <- igraph::E(x$graph)$ci_upper
    for (i in seq_len(nrow(el))) {
      a <- el[i, 1]; b <- el[i, 2]
      lo_mat[a, b] <- lo_mat[b, a] <- lo[i]
      hi_mat[a, b] <- hi_mat[b, a] <- hi[i]
    }
  }

  # Reorder by hierarchical clustering on (1 - effect size)
  if (reorder && p > 2) {
    fill_mat <- mat
    fill_mat[is.na(fill_mat)] <- 0
    diss  <- stats::as.dist(1 - fill_mat)
    hc    <- stats::hclust(diss, method = "average")
    ord   <- hc$order
    vars  <- vars[ord]
    mat   <- mat[ord, ord]
    pmat  <- pmat[ord, ord]
    if (has_ci) {
      lo_mat <- lo_mat[ord, ord]
      hi_mat <- hi_mat[ord, ord]
    }
  }

  # Long format
  rows <- vector("list", p * p)
  k    <- 0L
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      k        <- k + 1L
      v        <- mat[i, j]
      pv_ij    <- pmat[i, j]
      rows[[k]] <- data.frame(
        row_var  = vars[i],
        col_var  = vars[j],
        value    = v,
        p_value  = pv_ij,
        is_diag  = (i == j),
        label    = if (i == j) "" else
                     if (is.na(v)) "NA" else
                     formatC(v, digits = digits, format = "f"),
        stars    = if (i == j || !show_sig) "" else sig_stars(pv_ij),
        ci_lo    = if (has_ci) lo_mat[i, j] else NA_real_,
        ci_hi    = if (has_ci) hi_mat[i, j] else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  df_long <- do.call(rbind, rows)

  # Factor ordering for axes
  df_long$row_var <- factor(df_long$row_var, levels = rev(vars))
  df_long$col_var <- factor(df_long$col_var, levels = vars)

  # ------------------------------------------------- ggplot2 engine
  if (engine == "ggplot2") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop(
        "Package 'ggplot2' is required for engine = 'ggplot2'. ",
        "Install with: install.packages('ggplot2')"
      )
    }

    # Build CI label if needed
    if (show_ci && has_ci) {
      df_long$ci_label <- ifelse(
        df_long$is_diag | is.na(df_long$ci_lo), "",
        paste0("[", formatC(df_long$ci_lo, digits = digits, format = "f"),
               ", ", formatC(df_long$ci_hi, digits = digits, format = "f"), "]")
      )
    }

    # Diagonal tiles get their own fill
    df_long$fill_val <- ifelse(df_long$is_diag, NA_real_, df_long$value)

    p_out <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = .data$col_var, y = .data$row_var)
    ) +
      ggplot2::geom_tile(
        data = df_long[!df_long$is_diag, ],
        ggplot2::aes(fill = .data$fill_val),
        colour = "white", linewidth = 0.5
      ) +
      ggplot2::geom_tile(
        data = df_long[df_long$is_diag, ],
        fill = "#F1EFE8", colour = "white", linewidth = 0.5
      ) +
      ggplot2::scale_fill_gradientn(
        colours   = palette,
        limits    = c(0, 1),
        na.value  = na_fill,
        name      = if (x$corrected) "Effect size\n(corrected)" else "Effect size"
      )

    if (show_values) {
      p_out <- p_out +
        ggplot2::geom_text(
          data = df_long[!df_long$is_diag, ],
          ggplot2::aes(label = .data$label),
          size  = 3,
          vjust = if (show_sig || (show_ci && has_ci)) 0 else 0.5
        )
    }

    if (show_sig) {
      p_out <- p_out +
        ggplot2::geom_text(
          data = df_long[!df_long$is_diag, ],
          ggplot2::aes(label = .data$stars),
          size  = 3,
          vjust = if (show_values) 1.8 else 0.5,
          color = "#3C3489"
        )
    }

    if (show_ci && has_ci) {
      p_out <- p_out +
        ggplot2::geom_text(
          data = df_long[!df_long$is_diag, ],
          ggplot2::aes(label = .data$ci_label),
          size  = 2,
          vjust = if (show_values) 2.4 else 0.5,
          color = "#5F5E5A"
        )
    }

    p_out <- p_out +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title   = ggplot2::element_blank(),
        panel.grid   = ggplot2::element_blank(),
        legend.position = "right"
      )

    if (!is.null(title)) {
      p_out <- p_out + ggplot2::ggtitle(title)
    }

    return(p_out)
  }

  # ------------------------------------------------- base engine
  # Build numeric matrix for image()
  num_mat <- mat
  num_mat[is.na(num_mat)] <- 0
  num_mat[cbind(seq_len(p), seq_len(p))] <- NA_real_  # blank diagonal

  # Colour ramp
  n_col  <- 100L
  pal_fn <- grDevices::colorRampPalette(palette)
  cols   <- pal_fn(n_col)

  breaks <- seq(0, 1, length.out = n_col + 1L)

  # Flip for image() which plots row 1 at the bottom
  mat_flip <- num_mat[rev(seq_len(p)), ]

  op <- graphics::par(mar = c(6, 6, 4, 5))
  on.exit(graphics::par(op), add = TRUE)

  graphics::image(
    x    = seq_len(p),
    y    = seq_len(p),
    z    = t(mat_flip),
    col  = cols,
    axes = FALSE,
    xlab = "", ylab = "",
    zlim = c(0, 1),
    main = if (!is.null(title)) title else ""
  )

  graphics::axis(1, at = seq_len(p), labels = vars,
                 las = 2, cex.axis = 0.8, tick = FALSE)
  graphics::axis(2, at = seq_len(p), labels = rev(vars),
                 las = 1, cex.axis = 0.8, tick = FALSE)

  if (show_values) {
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next
        v <- mat[rev(seq_len(p))[i], vars[j]]
        if (!is.na(v)) {
          lbl <- formatC(v, digits = digits, format = "f")
          if (show_sig) {
            pv_ij <- pmat[rev(seq_len(p))[i], vars[j]]
            lbl   <- paste0(lbl, sig_stars(pv_ij))
          }
          graphics::text(j, i, labels = lbl, cex = 0.7)
        }
      }
    }
  }

  # Colour legend
  leg_x <- p + 0.6
  graphics::image(
    x    = c(leg_x, leg_x + 0.4),
    y    = seq(0.5, p + 0.5, length.out = n_col + 1L),
    z    = matrix(seq(0, 1, length.out = n_col), nrow = 1),
    col  = cols,
    add  = TRUE
  )
  graphics::text(leg_x + 0.2, c(0.5, p / 2, p + 0.5),
                 labels = c("0", "0.5", "1"), cex = 0.65, xpd = TRUE)
  graphics::mtext("Effect size", side = 4, line = 3.5, cex = 0.8)

  invisible(NULL)
}
