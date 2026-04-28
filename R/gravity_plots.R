# =============================================================================
# gravity_plots.R
# Visualisation functions for modality gravity indices
#
# Functions:
#   .gravity_bars()        -- internal bar chart helper for plot_gravity()
#   plot_gravity()         -- 2x3 network panel comparison
#   compare_gravity()      -- dMGI comparison across two subgroups
#   print.modality_gravity -- formatted role-grouped console output
#   summary.modality_gravity -- role summary and diagnostics
#   plot_gravity_scatter() -- eigenvector vs dMGI scatter diagnostic
# =============================================================================

# Internal helper: 2x3 bar chart panel for plot_gravity(bars = TRUE)
# Not exported. All args are pre-validated by plot_gravity().
.gravity_bars <- function(grav, betw, eig, comm, pal,
                          n_top, att_col, sat_col, main_title) {
  
  .abbrev <- function(x, n = 36L)
    ifelse(nchar(x) > n, paste0(substr(x, 1L, n - 1L), "..."), x)
  
  .one_bar <- function(vals, ttl, signed = FALSE) {
    d <- data.frame(
      modality  = grav$modality,
      community = as.character(comm),
      value     = vals,
      stringsAsFactors = FALSE
    )
    d <- d[order(-d$value), ]
    d <- d[seq_len(min(n_top, nrow(d))), ]
    d <- d[rev(seq_len(nrow(d))), ]   # largest at top for horiz bars
    
    bar_col <- if (signed) {
      ifelse(d$value >= 0,
             grDevices::adjustcolor(att_col, 0.82),
             grDevices::adjustcolor(sat_col, 0.82))
    } else {
      vapply(d$community, function(k)
        if (k %in% names(pal)) pal[[k]] else "grey70",
        character(1L))
    }
    
    bp <- graphics::barplot(
      d$value,
      horiz     = TRUE,
      names.arg = .abbrev(d$modality),
      las       = 1L,
      col       = bar_col,
      border    = NA,
      cex.names = 0.55,
      cex.axis  = 0.65,
      main      = ttl,
      cex.main  = 0.80,
      xlab      = ""
    )
    graphics::text(
      x      = d$value,
      y      = bp,
      labels = formatC(d$value, format = "f", digits = 3L),
      pos    = ifelse(d$value >= 0, 4L, 2L),
      cex    = 0.52,
      col    = "grey30",
      xpd    = NA
    )
    if (signed) graphics::abline(v = 0, lty = 2L, col = "grey45")
    invisible(bp)
  }
  
  op <- graphics::par(
    mfrow = c(2L, 3L),
    mar   = c(3, 17, 2.5, 5),
    oma   = c(2, 0, 3, 0),
    xpd   = FALSE
  )
  on.exit(graphics::par(op), add = TRUE)
  
  .one_bar(grav$strength,   "Strength (sum of phi)")
  .one_bar(betw,             "Betweenness (bridge role)")
  .one_bar(eig,              "Eigenvector (core position)")
  .one_bar(grav$mgi_plus,   "MGI+ (gravitational mass)")
  .one_bar(grav$os,          "OS (orbital pull intensity)")
  .one_bar(grav$delta_mgi,  "dMGI (net gravity)", signed = TRUE)
  
  graphics::mtext(
    paste(main_title, "- top", n_top, "per measure"),
    outer = TRUE, line = 1.0, cex = 0.95, font = 2L
  )
  
  # Community legend
  if (length(pal) > 1L) {
    graphics::par(fig = c(0, 1, 0, 0.07), new = TRUE,
                  mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::legend("center",
                     legend = names(pal),
                     fill   = pal,
                     border = NA,
                     bty    = "n",
                     horiz  = TRUE,
                     cex    = 0.80,
                     title  = "Community")
  }
  
  invisible(NULL)
}


# -----------------------------------------------------------------------------
#' Plot gravity indices alongside traditional centrality for a catmodgraph
#'
#' Produces a 2 x 3 panel figure showing six structural measures for the
#' same modality network on a shared layout: strength, betweenness, and
#' eigenvector centrality (top row) alongside MGI+, OS, and dMGI (bottom
#' row).  Node size encodes the magnitude of each measure; colour encodes
#' community membership (top row) or attractor/satellite role (dMGI panel).
#'
#' @param x A \code{catmodgraph} object.
#' @param gravity A data frame returned by \code{\link{modality_gravity}}.
#'   If \code{NULL} (default), \code{modality_gravity(x)} is called
#'   internally.
#' @param layout_fn An igraph layout function. Defaults to
#'   \code{igraph::layout_with_fr}.
#' @param community_attr Character. Vertex attribute name for community
#'   membership (set by \code{\link{cluster_modalities}}). Defaults to
#'   \code{"cluster"}.
#' @param palette Character vector of colours for communities. If
#'   \code{NULL}, \code{grDevices::hcl.colors} with palette \code{"Dark 3"}
#'   is used.
#' @param seed Integer. Random seed for reproducible layouts. Default
#'   \code{1L}.
#' @param title Character. Overall figure title. Default
#'   \code{"Modality network: traditional centrality vs. gravity indices"}.
#' @param node_size_range Numeric vector of length 2. Minimum and maximum
#'   node sizes. Default \code{c(4, 22)}.
#' @param show_labels Logical. Whether to draw node labels. Default
#'   \code{FALSE}.
#' @param bars Logical. If \code{TRUE}, a second figure is produced showing
#'   a 2 x 3 bar chart panel with the top \code{bars_n} modalities per
#'   measure, coloured by community.  Default \code{FALSE}.
#' @param bars_n Integer. Number of top modalities to show per bar chart
#'   panel when \code{bars = TRUE}.  Default \code{12L}.
#' @param attractor_col Colour for attractor nodes in the dMGI panel.
#'   Default \code{"#1D9E75"}.
#' @param satellite_col Colour for satellite nodes in the dMGI panel.
#'   Default \code{"#D85A30"}.
#' @param ... Additional arguments passed to \code{\link{modality_gravity}}
#'   when \code{gravity = NULL}.
#'
#' @return Invisibly returns the \code{gravity} data frame used for
#'   plotting.  The primary effect is the figure drawn on the current
#'   graphics device.
#'
#' @seealso \code{\link{modality_gravity}}, \code{\link{compare_gravity}}
#'
#' @examples
#' data(survey_health)
#' mg  <- build_modality_graph(survey_health)
#' mg  <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
#' mg  <- cluster_modalities(mg, method = "louvain")
#' plot_gravity(mg)
#'
#' @export
plot_gravity <- function(x,
                         gravity        = NULL,
                         layout_fn      = igraph::layout_with_fr,
                         community_attr = "cluster",
                         palette        = NULL,
                         seed           = 1L,
                         title    = "Modality network: traditional centrality vs. gravity indices",
                         node_size_range = c(4, 22),
                         show_labels    = FALSE,
                         bars           = FALSE,
                         bars_n         = 12L,
                         attractor_col  = "#1D9E75",
                         satellite_col  = "#D85A30",
                         ...) {
  
  # ---- Validation ----------------------------------------------------------
  if (!inherits(x, "catmodgraph")) {
    stop("`x` must be a catmodgraph object.", call. = FALSE)
  }
  
  g <- x$graph
  
  if (igraph::vcount(g) == 0L) {
    stop("Graph has no vertices.", call. = FALSE)
  }
  
  # ---- Gravity scores ------------------------------------------------------
  if (is.null(gravity)) {
    # Capture dots but remove plot_gravity-specific args before forwarding
    dots <- list(...)
    dots[c("bars", "bars_n")] <- NULL
    gravity <- do.call(modality_gravity, c(list(x = x), dots))
  }
  
  nodes <- igraph::V(g)$name
  
  # Align gravity rows to vertex order
  grav_aligned <- gravity[match(nodes, gravity$node), ]
  
  # ---- Community colours ---------------------------------------------------
  has_comm <- community_attr %in% igraph::vertex_attr_names(g)
  if (has_comm) {
    comm <- igraph::vertex_attr(g, community_attr)
    u_comm <- sort(unique(comm[!is.na(comm)]))
    if (is.null(palette)) {
      palette <- grDevices::hcl.colors(max(length(u_comm), 3L),
                                       palette = "Dark 3")
    }
    pal <- setNames(palette[seq_along(u_comm)], as.character(u_comm))
    comm_col <- pal[as.character(comm)]
    comm_col[is.na(comm_col)] <- "grey70"
  } else {
    comm_col <- rep("grey70", igraph::vcount(g))
  }
  
  # ---- Layout --------------------------------------------------------------
  set.seed(seed)
  lay <- layout_fn(g)
  lay[, 1L] <- scales::rescale(lay[, 1L], to = c(-1, 1))
  lay[, 2L] <- scales::rescale(lay[, 2L], to = c(-1, 1))
  
  # ---- Betweenness & eigenvector -------------------------------------------
  w_abs <- abs(igraph::edge_attr(g, "weight"))
  betw  <- igraph::betweenness(g,
                               weights    = 1 / pmax(w_abs, 1e-9),
                               normalized = TRUE)
  eig   <- igraph::eigen_centrality(g, weights = w_abs)$vector
  
  # ---- Helper: scale values to node size range ----------------------------
  .scale_size <- function(v, log_compress = FALSE) {
    if (log_compress) v <- log1p(abs(v))
    rng <- range(v, na.rm = TRUE)
    if (diff(rng) < 1e-12) return(rep(mean(node_size_range), length(v)))
    node_size_range[1L] +
      (v - rng[1L]) / (rng[2L] - rng[1L]) *
      diff(node_size_range)
  }
  
  # ---- Helper: alpha scaling within a role group --------------------------
  .role_alpha <- function(v, role_mask, alpha_range = c(0.40, 0.95)) {
    sub_v <- abs(v[role_mask])
    rng   <- range(sub_v, na.rm = TRUE)
    if (diff(rng) < 1e-12) return(rep(mean(alpha_range), sum(role_mask)))
    alpha_range[1L] + (sub_v - rng[1L]) / diff(rng) * diff(alpha_range)
  }
  
  # ---- Helper: draw one panel ---------------------------------------------
  .one_panel <- function(sz, col, ttl, frame_col = "white") {
    plot(g,
         layout             = lay,
         vertex.color       = col,
         vertex.size        = sz,
         vertex.frame.color = frame_col,
         vertex.label       = if (show_labels) nodes else NA,
         vertex.label.cex   = 0.55,
         edge.width         = pmax(w_abs, 0.05) * 5,
         edge.color         = grDevices::adjustcolor("grey50", 0.28),
         edge.curved        = 0,
         main               = ttl,
         cex.main           = 0.82,
         rescale            = FALSE,
         xlim               = c(-1.05, 1.05),
         ylim               = c(-1.05, 1.05))
  }
  
  # ---- Build dMGI colours (signed, magnitude-alpha) -----------------------
  delta   <- grav_aligned$delta_mgi
  is_att  <- !is.na(delta) & delta > 0
  is_sat  <- !is.na(delta) & delta < 0
  alpha_v <- rep(0.55, length(delta))
  if (any(is_att)) alpha_v[is_att] <- .role_alpha(delta, is_att,
                                                  c(0.40, 0.95))
  if (any(is_sat)) alpha_v[is_sat] <- .role_alpha(delta, is_sat,
                                                  c(0.55, 0.90))
  dmgi_col <- mapply(function(role, a) {
    base <- if (role == "attractor") attractor_col
    else if (role == "satellite") satellite_col
    else "grey65"
    grDevices::adjustcolor(base, alpha.f = a)
  }, grav_aligned$role, alpha_v)
  
  dmgi_frame <- ifelse(grav_aligned$role == "satellite",
                       grDevices::adjustcolor(satellite_col, 0.7),
                       grDevices::adjustcolor("white", 0.4))
  
  # ---- Draw 2x3 panel ------------------------------------------------------
  op <- graphics::par(mfrow = c(2L, 3L),
                      mar   = c(1, 1, 2.8, 1),
                      oma   = c(if (has_comm) 4 else 1, 0, 3, 0))
  on.exit(graphics::par(op), add = TRUE)
  
  # Row 1 - traditional
  .one_panel(.scale_size(grav_aligned$strength),
             comm_col, "Strength (sum of phi)")
  
  .one_panel(.scale_size(betw),
             comm_col, "Betweenness (bridge role)")
  
  .one_panel(.scale_size(eig),
             comm_col, "Eigenvector (core position)")
  
  # Row 2 - gravity
  .one_panel(.scale_size(grav_aligned$mgi_plus, log_compress = TRUE),
             comm_col, "MGI+ (gravitational mass)")
  
  .one_panel(.scale_size(grav_aligned$os, log_compress = TRUE),
             comm_col, "OS (orbital pull intensity)")
  
  .one_panel(.scale_size(delta, log_compress = TRUE),
             dmgi_col,
             "dMGI (net gravity: green=attractor, red=satellite)",
             frame_col = dmgi_frame)
  
  graphics::mtext(title, outer = TRUE, line = 1.2,
                  cex = 1.0, font = 2L)
  
  # ---- Community legend at the bottom -------------------------------------
  if (has_comm) {
    graphics::par(fig     = c(0, 1, 0, 0.08),
                  new     = TRUE,
                  mar     = c(0, 0, 0, 0),
                  oma     = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::legend("center",
                     legend = as.character(u_comm),
                     fill   = pal,
                     border = NA,
                     bty    = "n",
                     horiz  = TRUE,
                     cex    = 0.82,
                     title  = "Community")
  }
  
  # ---- Optional bar chart panel --------------------------------------------
  if (isTRUE(bars)) {
    .gravity_bars(
      grav      = grav_aligned,
      betw      = betw,
      eig       = eig,
      comm      = if (has_comm) comm else rep(1L, igraph::vcount(g)),
      pal       = if (has_comm) pal else c("1" = "grey70"),
      n_top     = bars_n,
      att_col   = attractor_col,
      sat_col   = satellite_col,
      main_title = title
    )
  }
  
  invisible(gravity)
}
#'
#' Computes \code{\link{modality_gravity}} on two conditional modality graphs
#' and returns a side-by-side comparison of dMGI, OS, and role for every
#' modality present in either graph.  Optionally plots a dot-chart of dMGI
#' differences.
#'
#' @param x A named list of exactly two \code{catmodgraph} objects, e.g.
#'   \code{list(women = mg_f, men = mg_m)}.  Names are used as group labels.
#' @param plot Logical. Whether to draw a comparison dot-chart. Default
#'   \code{TRUE}.
#' @param top_n Integer. Number of modalities to show in the plot (those with
#'   the largest absolute dMGI difference). Default \code{20L}.
#' @param ... Additional arguments passed to \code{\link{modality_gravity}}.
#'
#' @return A data frame with columns \code{node}, \code{variable},
#'   \code{modality}, and for each group \code{g1}/\code{g2}:
#'   \code{prevalence_g1}, \code{delta_mgi_g1}, \code{os_g1},
#'   \code{role_g1} (and equivalents for g2), plus \code{delta_mgi_diff}
#'   (\code{g1 - g2}).  Rows are ordered by \code{abs(delta_mgi_diff)}
#'   descending.
#'
#' @seealso \code{\link{modality_gravity}}, \code{\link{plot_gravity}},
#'   \code{\link{compare_modality_graphs}}
#'
#' @examples
#' \donttest{
#' data(survey_health)
#' mg_f <- build_conditional_modality_graph(survey_health,
#'           given = list(sex = "female"))
#' mg_m <- build_conditional_modality_graph(survey_health,
#'           given = list(sex = "male"))
#' mg_f <- prune_modality_edges(mg_f, min_weight = 0.10, max_p = 0.05)
#' mg_m <- prune_modality_edges(mg_m, min_weight = 0.10, max_p = 0.05)
#' cmp  <- compare_gravity(list(female = mg_f, male = mg_m))
#' print(head(cmp, 10))
#' }
#'
#' @export
compare_gravity <- function(x, plot = TRUE, top_n = 20L, ...) {
  
  # ---- Validation ----------------------------------------------------------
  if (!is.list(x) || length(x) != 2L) {
    stop("`x` must be a named list of exactly two catmodgraph objects.",
         call. = FALSE)
  }
  if (is.null(names(x)) || any(names(x) == "")) {
    stop("`x` must be a named list, e.g. list(female = mg_f, male = mg_m).",
         call. = FALSE)
  }
  if (!all(vapply(x, inherits, logical(1L), "catmodgraph"))) {
    stop("Both elements of `x` must be catmodgraph objects.", call. = FALSE)
  }
  
  g1_name <- names(x)[1L]
  g2_name <- names(x)[2L]
  
  # ---- Compute gravity for each group -------------------------------------
  grav1 <- modality_gravity(x[[1L]], ...)
  grav2 <- modality_gravity(x[[2L]], ...)
  
  # ---- Merge on node -------------------------------------------------------
  cols_keep <- c("node", "variable", "modality",
                 "prevalence", "delta_mgi", "os", "role")
  
  m <- merge(
    grav1[, cols_keep],
    grav2[, cols_keep],
    by     = c("node", "variable", "modality"),
    all    = TRUE,
    suffixes = paste0("_", c(g1_name, g2_name))
  )
  
  # Replace NAs introduced by full join (node absent in one graph)
  num_cols <- grep("prevalence_|delta_mgi_|os_", names(m), value = TRUE)
  for (col in num_cols) {
    m[[col]][is.na(m[[col]])] <- 0
  }
  role_cols <- grep("^role_", names(m), value = TRUE)
  for (col in role_cols) {
    m[[col]][is.na(m[[col]])] <- "absent"
  }
  
  diff_col <- paste0("delta_mgi_", g1_name)
  diff_col2 <- paste0("delta_mgi_", g2_name)
  m$delta_mgi_diff <- m[[diff_col]] - m[[diff_col2]]
  
  m <- m[order(-abs(m$delta_mgi_diff)), ]
  rownames(m) <- NULL
  
  # ---- Optional plot -------------------------------------------------------
  if (plot) {
    n_show  <- min(top_n, nrow(m))
    m_plot  <- m[seq_len(n_show), ]
    m_plot  <- m_plot[order(m_plot$delta_mgi_diff), ]
    
    diff_vals <- m_plot$delta_mgi_diff
    bar_col   <- ifelse(diff_vals > 0,
                        grDevices::adjustcolor("#1D9E75", 0.80),
                        grDevices::adjustcolor("#D85A30", 0.80))
    
    op <- graphics::par(mar = c(5, 18, 4, 2))
    on.exit(graphics::par(op), add = TRUE)
    
    graphics::barplot(
      diff_vals,
      horiz     = TRUE,
      names.arg = m_plot$modality,
      las       = 1L,
      col       = bar_col,
      border    = NA,
      cex.names = 0.65,
      xlab      = "dMGI(g1) - dMGI(g2)",
      main      = sprintf(
        "dMGI difference: %s (green) vs %s (red)\ntop %d modalities by |diff|",
        g1_name, g2_name, n_show)
    )
    graphics::abline(v = 0, lty = 2L, col = "grey40")
  }
  
  invisible(m)
}


# -----------------------------------------------------------------------------
#' Print method for modality_gravity output
#'
#' Displays a formatted, role-grouped summary of the data frame returned by
#' \code{\link{modality_gravity}}.  Attractors are shown first in descending
#' \code{delta_mgi} order, then neutrals, then satellites in ascending order
#' (strongest satellite last).
#'
#' @param x A data frame returned by \code{\link{modality_gravity}}.
#' @param digits Integer. Number of decimal places for numeric columns.
#'   Default \code{3L}.
#' @param max_width Integer. Maximum characters for the node label column
#'   before truncation with \code{...}. Default \code{45L}.
#' @param ... Ignored. Present for S3 compatibility.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{modality_gravity}}, \code{\link{summary.modality_gravity}}
#'
#' @examples
#' data(survey_health)
#' mg   <- build_modality_graph(survey_health)
#' mg   <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
#' grav <- modality_gravity(mg)
#' print(grav)
#'
#' @export
print.modality_gravity <- function(x, digits = 3L, max_width = 45L, ...) {
  
  # If x has been subsetted and lost required columns, fall back to data.frame print
  required <- c("node", "prevalence", "mgi_plus", "mgi_minus",
                "delta_mgi", "os", "role")
  if (!all(required %in% names(x))) {
    class(x) <- "data.frame"
    print(x, ...)
    return(invisible(x))
  }
  
  .trunc <- function(s, w) {
    ifelse(nchar(s) > w, paste0(substr(s, 1L, w - 1L), "..."), s)
  }
  
  .fmt_num <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    formatC(round(v, digits), format = "f", digits = digits)
  }
  
  .section <- function(df, header, col) {
    if (nrow(df) == 0L) return(invisible(NULL))
    cat(header, "\n", sep = "")
    cat(sprintf("  %-*s  %8s  %8s  %8s  %8s  %6s\n",
                max_width, "Node",
                "Prev", "MGI+", "MGI-", "dMGI", "OS"))
    cat(paste(rep("-", max_width + 48L), collapse = ""), "\n")
    for (i in seq_len(nrow(df))) {
      cat(sprintf("  %-*s  %8s  %8s  %8s  %8s  %6s\n",
                  max_width,
                  .trunc(df$node[i], max_width),
                  .fmt_num(df$prevalence[i]),
                  .fmt_num(df$mgi_plus[i]),
                  .fmt_num(df$mgi_minus[i]),
                  .fmt_num(df$delta_mgi[i]),
                  .fmt_num(df$os[i])))
    }
    cat("\n")
  }
  
  att <- x[x$role == "attractor", ]
  neu <- x[x$role == "neutral",   ]
  sat <- x[x$role == "satellite", ]
  sat <- sat[order(sat$delta_mgi), ]   # strongest satellite last
  
  cat("\n")
  cat("=== Modality Gravity Index ===\n")
  cat(sprintf("(%d nodes: %d attractors, %d neutral, %d satellites)\n\n",
              nrow(x), nrow(att), nrow(neu), nrow(sat)))
  
  .section(att, "ATTRACTORS  [dMGI > 0]", "")
  .section(neu, "NEUTRAL     [dMGI = 0]", "")
  .section(sat, "SATELLITES  [dMGI < 0]", "")
  
  cat(paste0("Columns: Prev = prevalence  MGI+ = gravitational mass",
             "  MGI- = orbital drag  dMGI = net gravity  OS = orbital score\n"))
  
  invisible(x)
}


# -----------------------------------------------------------------------------
#' Summary method for modality_gravity output
#'
#' Returns a compact role-distribution table and per-variable breakdown,
#' and prints a readable synopsis to the console.
#'
#' @param object A data frame returned by \code{\link{modality_gravity}}.
#' @param ... Ignored.
#'
#' @return Invisibly returns a list with elements:
#'   \describe{
#'     \item{role_counts}{Named integer vector of attractor/neutral/satellite
#'       counts.}
#'     \item{by_variable}{Data frame with one row per variable showing the
#'       dominant role and mean \code{delta_mgi} across its modalities.}
#'     \item{top_attractor}{Single-row data frame: modality with highest
#'       \code{delta_mgi}.}
#'     \item{top_satellite}{Single-row data frame: modality with lowest
#'       \code{delta_mgi}.}
#'     \item{spearman_strength}{Spearman correlation between \code{strength}
#'       and \code{delta_mgi}.}
#'   }
#'
#' @seealso \code{\link{modality_gravity}}, \code{\link{print.modality_gravity}}
#'
#' @examples
#' data(survey_health)
#' mg   <- build_modality_graph(survey_health)
#' mg   <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
#' grav <- modality_gravity(mg)
#' summary(grav)
#'
#' @export
summary.modality_gravity <- function(object, ...) {
  
  role_counts <- table(object$role)
  
  by_var <- do.call(rbind, lapply(
    split(object, object$variable), function(df) {
      dominant <- names(sort(table(df$role), decreasing = TRUE))[1L]
      data.frame(
        variable     = df$variable[1L],
        n_modalities = nrow(df),
        mean_delta   = round(mean(df$delta_mgi, na.rm = TRUE), 3L),
        dominant_role = dominant,
        stringsAsFactors = FALSE
      )
    }
  ))
  by_var <- by_var[order(-by_var$mean_delta), ]
  rownames(by_var) <- NULL
  
  top_att <- object[which.max(object$delta_mgi), ]
  top_sat <- object[which.min(object$delta_mgi), ]
  
  rho_str <- if ("strength" %in% names(object))
    round(stats::cor(object$strength, object$delta_mgi,
                     method = "spearman", use = "complete.obs"), 3L)
  else NA_real_
  
  cat("\n=== Modality Gravity - Summary ===\n\n")
  
  cat("Role distribution:\n")
  print(role_counts)
  
  cat(sprintf("\nDominant attractor : %s  (dMGI = %.3f, prev = %.3f)\n",
              top_att$node, top_att$delta_mgi, top_att$prevalence))
  cat(sprintf("Strongest satellite: %s  (dMGI = %.3f, prev = %.3f)\n",
              top_sat$node, top_sat$delta_mgi, top_sat$prevalence))
  
  cat(sprintf("\nSpearman rho (strength vs dMGI): %.3f\n", rho_str))
  cat("(Values far from +/-1.0 indicate MGI captures structure",
      "beyond connectivity)\n")
  
  cat("\nPer-variable breakdown:\n")
  print(by_var, row.names = FALSE)
  
  invisible(list(
    role_counts      = role_counts,
    by_variable      = by_var,
    top_attractor    = top_att,
    top_satellite    = top_sat,
    spearman_strength = rho_str
  ))
}


# -----------------------------------------------------------------------------
#' Scatter plot of eigenvector centrality vs dMGI
#'
#' Plots eigenvector centrality on the x-axis against net gravity (dMGI) on
#' the y-axis for all modality nodes.  Points are coloured by role
#' (attractor/satellite/neutral) and contradiction cases - nodes where the two
#' metrics disagree most strongly - are automatically labelled.  The Spearman
#' correlation is annotated in the plot margin.
#'
#' This plot is the primary diagnostic for demonstrating that MGI captures
#' structural information not contained in standard centrality indices.
#'
#' @param x A data frame returned by \code{\link{modality_gravity}}.
#' @param catmodgraph A \code{catmodgraph} object matching \code{x}.  Used to
#'   compute eigenvector centrality if not already present as a column in
#'   \code{x}.
#' @param label_threshold Numeric in \[0, 1\]. Modalities are labelled if
#'   their eigenvector centrality exceeds this value and they are satellites,
#'   OR if their \code{delta_mgi} exceeds the 75th percentile and their
#'   eigenvector is below the median.  Default \code{0.25}.
#' @param attractor_col Colour for attractor nodes. Default \code{"#1D9E75"}.
#' @param satellite_col Colour for satellite nodes. Default \code{"#D85A30"}.
#' @param neutral_col Colour for neutral nodes. Default \code{"grey60"}.
#' @param point_size Numeric. Base point size. Default \code{1.8}.
#'
#' @return Invisibly returns a data frame with columns \code{node},
#'   \code{eigenvec}, \code{delta_mgi}, \code{role}, and
#'   \code{is_contradiction}.
#'
#' @seealso \code{\link{modality_gravity}}, \code{\link{plot_gravity}}
#'
#' @examples
#' data(survey_health)
#' mg   <- build_modality_graph(survey_health)
#' mg   <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
#' grav <- modality_gravity(mg)
#' plot_gravity_scatter(grav, mg)
#'
#' @export
plot_gravity_scatter <- function(x,
                                 catmodgraph,
                                 label_threshold = 0.25,
                                 attractor_col   = "#1D9E75",
                                 satellite_col   = "#D85A30",
                                 neutral_col     = "grey60",
                                 point_size      = 1.8) {
  
  if (!inherits(catmodgraph, "catmodgraph")) {
    stop("`catmodgraph` must be a catmodgraph object.", call. = FALSE)
  }
  
  g     <- catmodgraph$graph
  w_abs <- abs(igraph::edge_attr(g, "weight"))
  eig   <- igraph::eigen_centrality(g, weights = w_abs)$vector
  
  # Align to gravity row order
  eig_aligned <- eig[match(x$node, names(eig))]
  eig_aligned[is.na(eig_aligned)] <- 0
  
  df <- data.frame(
    node      = x$node,
    modality  = x$modality,
    eigenvec  = round(eig_aligned, 4L),
    delta_mgi = x$delta_mgi,
    os        = x$os,
    role      = x$role,
    stringsAsFactors = FALSE
  )
  
  # Spearman correlation
  rho <- stats::cor(df$eigenvec, df$delta_mgi,
                    method = "spearman", use = "complete.obs")
  
  # Contradiction flags
  med_eig  <- stats::median(df$eigenvec, na.rm = TRUE)
  q75_dmgi <- stats::quantile(df$delta_mgi, 0.75, na.rm = TRUE)
  
  df$is_contradiction <-
    (df$eigenvec > label_threshold & df$role == "satellite") |
    (df$delta_mgi > q75_dmgi       & df$eigenvec < med_eig)
  
  # Point colours
  pt_col <- ifelse(df$role == "attractor", attractor_col,
                   ifelse(df$role == "satellite", satellite_col,
                          neutral_col))
  
  # Plot
  op <- graphics::par(mar = c(5, 5, 4, 2))
  on.exit(graphics::par(op), add = TRUE)
  
  graphics::plot(
    df$eigenvec, df$delta_mgi,
    col  = pt_col,
    pch  = 19,
    cex  = point_size,
    xlab = "Eigenvector centrality",
    ylab = "dMGI (net gravity)",
    main = "Eigenvector centrality vs. dMGI"
  )
  
  graphics::abline(h = 0,
                   v = stats::median(df$eigenvec, na.rm = TRUE),
                   lty = 2L, col = "grey55")
  
  # Annotate contradictions
  contr <- df[df$is_contradiction, ]
  if (nrow(contr) > 0L) {
    graphics::text(
      contr$eigenvec, contr$delta_mgi,
      labels = contr$modality,
      cex    = 0.62,
      pos    = 3L,
      col    = "grey20"
    )
  }
  
  # Spearman annotation - ASCII safe
  graphics::mtext(
    sprintf("Spearman rho = %.3f   R2 = %.3f", rho, rho^2),
    side = 3L, line = 0.3, cex = 0.85, col = "grey30"
  )
  
  graphics::legend(
    "topright",
    legend  = c("attractor", "neutral", "satellite"),
    col     = c(attractor_col, neutral_col, satellite_col),
    pch     = 19L,
    pt.cex  = 1.3,
    bty     = "n",
    cex     = 0.85
  )
  
  invisible(df)
}