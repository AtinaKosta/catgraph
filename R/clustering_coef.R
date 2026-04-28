#' Compute weighted clustering coefficients for all variables in a catgraph
#'
#' Returns a data frame of five local clustering coefficients for every node
#' in a \code{catgraph} object. Four are established weighted extensions of
#' the Watts-Strogatz coefficient; the fifth (redundancy) is a measure
#' specific to effect-size graphs that quantifies how much of a pairwise
#' association is explained by indirect paths through other variables.
#'
#' @param x A \code{catgraph} object.
#' @param method Character vector. Which coefficient(s) to compute. Any
#'   subset of \code{c("watts_strogatz", "barrat", "onnela", "zhang",
#'   "redundancy")} or \code{"all"} (default).
#' @param normalize Logical. If \code{TRUE} (default), each coefficient is
#'   normalised to [0, 1] by dividing by its theoretical maximum. For the
#'   four classical coefficients the theoretical maximum is already 1.
#'   For redundancy the observed maximum is used.
#'
#' @return A \code{data.frame} with one row per variable. Columns are
#'   \code{variable} plus one column per requested method. The data frame
#'   is sorted by the mean clustering coefficient across all requested
#'   methods in descending order.
#'
#' @details
#' All five coefficients measure the extent to which the neighbours of a
#' node are also connected to each other — the triangle-closing tendency —
#' but they differ in how they incorporate edge weights.
#'
#' \strong{Watts-Strogatz (unweighted; Watts & Strogatz, 1998):}
#'
#' \deqn{C_i^{WS} = \frac{t_i}{k_i(k_i-1)}}
#'
#' where \eqn{t_i} is the number of triangles through node \eqn{i} and
#' \eqn{k_i} is its degree. This is the classical unweighted coefficient
#' included as a baseline.
#'
#' \strong{Barrat et al. (2004):}
#'
#' \deqn{C_i^{B} = \frac{1}{s_i(k_i-1)}
#'   \sum_{j,h} \frac{w_{ij}+w_{ih}}{2} a_{ij} a_{ih} a_{jh}}
#'
#' where \eqn{s_i = \sum_j w_{ij}} is the node strength and \eqn{a_{ij}}
#' is the binary adjacency. Weights the contribution of each triangle by
#' the average weight of the two edges incident to the focal node. Nodes
#' that close triangles through their strongest edges score highly.
#'
#' \strong{Onnela et al. (2005):}
#'
#' \deqn{C_i^{O} = \frac{1}{k_i(k_i-1)}
#'   \sum_{j,h} (\hat{w}_{ij} \hat{w}_{ih} \hat{w}_{jh})^{1/3}}
#'
#' where \eqn{\hat{w} = w / \max(w)} are normalised weights. Uses the
#' geometric mean of triangle edge weights, making it sensitive to weak
#' ties: a triangle with one very weak edge scores low.
#'
#' \strong{Zhang & Horvath (2005):}
#'
#' \deqn{C_i^{ZH} = \frac{\sum_{j,h} w_{ij} w_{ih} w_{jh}}
#'   {\left(\sum_{j \neq i} w_{ij}\right)^2 -
#'    \sum_{j \neq i} w_{ij}^2}}
#'
#' The numerator sums products of three edge weights around all triangles;
#' the denominator is the squared strength minus the sum of squared weights.
#' Originally proposed for co-expression networks (WGCNA). Amplifies
#' strong-edge triangles via cubic weighting.
#'
#' \strong{Redundancy (original):}
#'
#' For each edge \eqn{(i,j)} in the graph, the redundancy ratio is defined
#' as the direct weight divided by the maximum indirect path weight:
#'
#' \deqn{R_{ij} = \frac{w_{ij}}{\max_{k \neq i,j}
#'   \min(w_{ik}, w_{kj})}}
#'
#' The node-level redundancy coefficient is the mean of \eqn{R_{ij}} over
#' all edges incident to node \eqn{i}. A value near 1 indicates that the
#' direct association between two variables is no stronger than what their
#' shared neighbours would predict — the edge may be mediated. A value
#' well above 1 indicates a genuine direct association that exceeds
#' indirect explanation. When no indirect path exists (k < 3), the edge
#' is assigned \eqn{R_{ij} = \infty} and excluded from the node mean.
#'
#' @references
#' Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A. (2004).
#'   The architecture of complex weighted networks.
#'   \emph{Proceedings of the National Academy of Sciences}, 101(11),
#'   3747--3752. \doi{10.1073/pnas.0400087101}
#'
#' Onnela, J.-P., Saramaki, J., Kertesz, J., & Kaski, K. (2005).
#'   Intensity and coherence of motifs in weighted complex networks.
#'   \emph{Physical Review E}, 71(6), 065103.
#'   \doi{10.1103/PhysRevE.71.065103}
#'
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
#'   'small-world' networks. \emph{Nature}, 393(6684), 440--442.
#'   \doi{10.1038/30918}
#'
#' Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
#'   co-expression network analysis. \emph{Statistical Applications in
#'   Genetics and Molecular Biology}, 4(1), Article 17.
#'   \doi{10.2202/1544-6115.1128}
#'
#' @examples
#' df <- expand_table(Titanic)
#' cg <- catgraph(df)
#'
#' cc <- clustering_coef(cg)
#' cc
#'
#' # Single method
#' clustering_coef(cg, method = "barrat")
#'
#' # Compare all methods
#' compare_clustering(cg)
#'
#' @seealso \code{\link{compare_clustering}}, \code{\link{plot_clustering}},
#'   \code{\link{node_centrality}}
#' @importFrom igraph as_adjacency_matrix V E vcount
#' @export
clustering_coef <- function(x,
                             method    = "all",
                             normalize = TRUE) {

  if (!inherits(x, "catgraph")) stop("`x` must be a catgraph object.")

  valid <- c("watts_strogatz", "barrat", "onnela", "zhang", "redundancy")
  if (identical(method, "all")) method <- valid
  bad <- setdiff(method, valid)
  if (length(bad) > 0L) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "),
         ". Choose from: ", paste(valid, collapse = ", "))
  }

  g    <- x$graph
  p    <- igraph::vcount(g)
  vars <- igraph::V(g)$name

  # Weight matrix (symmetric, zero diagonal)
  W <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight",
                                              sparse = FALSE))
  W[is.na(W)] <- 0
  diag(W)     <- 0

  # Binary adjacency: 1 where an edge exists, 0 otherwise. Post v0.4.0,
  # absent edges are stored as exactly 0 in W (not as a synthetic epsilon),
  # so the threshold here is purely a float-safe "strictly positive" test.
  A <- (W > 0) * 1.0

  # Degree vector
  k <- rowSums(A)

  # Strength vector
  s <- rowSums(W)

  out <- data.frame(variable = vars, stringsAsFactors = FALSE)

  # ------------------------------------------------------ Watts-Strogatz
  if ("watts_strogatz" %in% method) {
    cc_ws <- numeric(p)
    for (i in seq_len(p)) {
      if (k[i] < 2) { cc_ws[i] <- NA_real_; next }
      nb   <- which(A[i, ] > 0)
      denom <- k[i] * (k[i] - 1)
      numer <- sum(A[nb, nb])
      cc_ws[i] <- numer / denom
    }
    out$watts_strogatz <- cc_ws
  }

  # ------------------------------------------------------ Barrat
  if ("barrat" %in% method) {
    cc_b <- numeric(p)
    for (i in seq_len(p)) {
      if (k[i] < 2) { cc_b[i] <- NA_real_; next }
      nb    <- which(A[i, ] > 0)
      denom <- s[i] * (k[i] - 1)
      if (denom < .Machine$double.eps) { cc_b[i] <- 0; next }
      numer <- 0
      for (j in nb) {
        for (h in nb) {
          if (j >= h) next
          numer <- numer + ((W[i,j] + W[i,h]) / 2) * A[j,h]
        }
      }
      cc_b[i] <- numer / denom
    }
    out$barrat <- cc_b
  }

  # ------------------------------------------------------ Onnela
  if ("onnela" %in% method) {
    w_max  <- max(W, na.rm = TRUE)
    W_norm <- if (w_max > 0) W / w_max else W

    cc_o <- numeric(p)
    for (i in seq_len(p)) {
      if (k[i] < 2) { cc_o[i] <- NA_real_; next }
      nb    <- which(A[i, ] > 0)
      denom <- k[i] * (k[i] - 1)
      numer <- 0
      for (j in nb) {
        for (h in nb) {
          if (j >= h) next
          if (A[j,h] > 0) {
            numer <- numer +
              (W_norm[i,j] * W_norm[i,h] * W_norm[j,h])^(1/3)
          }
        }
      }
      cc_o[i] <- numer / denom
    }
    out$onnela <- cc_o
  }

  # ------------------------------------------------------ Zhang-Horvath
  # Vectorised: the triple-loop numerator sum_{j != i, h != i, h != j}
  # W_ij W_jh W_hi equals (W %*% W %*% W)[i, i] when diag(W) = 0, which is
  # enforced above. This turns an O(p^3) R loop into a single O(p^3) BLAS
  # call — asymptotically the same but orders of magnitude faster in
  # practice for p > ~15.
  if ("zhang" %in% method) {
    W2   <- W %*% W
    W3   <- W2 %*% W
    numerator_vec <- diag(W3)
    
    denominator_vec <- s^2 - rowSums(W^2)
    
    cc_zh <- ifelse(
      k < 2,
      NA_real_,
      ifelse(abs(denominator_vec) < .Machine$double.eps,
             0,
             numerator_vec / denominator_vec)
    )
    out$zhang <- cc_zh
  }

  # ------------------------------------------------------ Redundancy
  if ("redundancy" %in% method) {
    red_node <- numeric(p)
    for (i in seq_len(p)) {
      nb <- which(A[i, ] > 0)
      if (length(nb) == 0L) { red_node[i] <- NA_real_; next }

      ratios <- numeric(length(nb))
      for (idx in seq_along(nb)) {
        j       <- nb[idx]
        w_direct <- W[i, j]

        # Find all indirect paths i -> k -> j (k != i, k != j)
        third_nodes <- setdiff(which(A[i, ] > 0 & A[j, ] > 0), c(i, j))

        if (length(third_nodes) == 0L) {
          # No indirect path — exclude from mean (NA)
          ratios[idx] <- NA_real_
        } else {
          max_indirect <- max(vapply(third_nodes, function(k) {
            min(W[i, k], W[k, j])
          }, numeric(1L)))
          ratios[idx] <- if (max_indirect < .Machine$double.eps) NA_real_
                         else w_direct / max_indirect
        }
      }
      valid_ratios    <- ratios[!is.na(ratios)]
      red_node[i] <- if (length(valid_ratios) == 0L) NA_real_
                     else mean(valid_ratios)
    }

    # Normalise by observed max (redundancy has no fixed upper bound)
    if (normalize) {
      mx <- max(red_node, na.rm = TRUE)
      if (!is.na(mx) && mx > 0) red_node <- red_node / mx
    }
    out$redundancy <- red_node
  }

  # Sort by mean CC across all computed methods
  cc_cols <- intersect(method, names(out))
  if (length(cc_cols) > 0L) {
    row_means <- rowMeans(out[, cc_cols, drop = FALSE], na.rm = TRUE)
    out <- out[order(row_means, decreasing = TRUE), ]
  }
  rownames(out) <- NULL
  out
}


#' Compare all weighted clustering coefficients side by side
#'
#' A convenience wrapper around \code{\link{clustering_coef}} that always
#' computes all five methods and returns them in a single wide data frame,
#' making it easy to identify variables that score consistently high (robust
#' cluster hubs) versus those that vary across methods (structurally
#' ambiguous variables).
#'
#' @param x A \code{catgraph} object.
#' @param normalize Logical. Default \code{TRUE}.
#'
#' @return A \code{data.frame} with columns \code{variable},
#'   \code{watts_strogatz}, \code{barrat}, \code{onnela}, \code{zhang},
#'   \code{redundancy}, and \code{mean_cc} (the row mean across all five
#'   methods, excluding \code{NA}s). Sorted by \code{mean_cc} descending.
#'
#' @examples
#' df <- expand_table(Titanic)
#' cg <- catgraph(df)
#' compare_clustering(cg)
#'
#' @seealso \code{\link{clustering_coef}}, \code{\link{plot_clustering}}
#' @export
compare_clustering <- function(x, normalize = TRUE) {
  cc  <- clustering_coef(x, method = "all", normalize = normalize)
  cc_cols <- c("watts_strogatz", "barrat", "onnela", "zhang", "redundancy")
  cc_cols <- intersect(cc_cols, names(cc))
  cc$mean_cc <- rowMeans(cc[, cc_cols, drop = FALSE], na.rm = TRUE)
  cc <- cc[order(cc$mean_cc, decreasing = TRUE), ]
  rownames(cc) <- NULL
  cc
}


#' Plot weighted clustering coefficients for a catgraph
#'
#' Produces a bar chart (single method) or a grouped/faceted comparison
#' chart (multiple methods) of clustering coefficients.
#'
#' @param x A \code{catgraph} object.
#' @param method Character. One method name or \code{"all"}. Default
#'   \code{"barrat"}.
#' @param normalize Logical. Default \code{TRUE}.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param engine Character. \code{"ggplot2"} (default) or \code{"base"}.
#'
#' @return For \code{engine = "ggplot2"}: a \code{ggplot} object.
#'   For \code{engine = "base"}: \code{NULL}, invisibly.
#'
#' @examples
#' df <- expand_table(Titanic)
#' cg <- catgraph(df)
#' plot_clustering(cg)
#' plot_clustering(cg, method = "all")
#'
#' @seealso \code{\link{clustering_coef}}, \code{\link{compare_clustering}}
#' @importFrom graphics par barplot
#' @importFrom stats reshape
#' @export
plot_clustering <- function(x,
                             method    = "barrat",
                             normalize = TRUE,
                             title     = NULL,
                             engine    = c("ggplot2", "base")) {

  engine <- match.arg(engine)
  valid  <- c("watts_strogatz", "barrat", "onnela", "zhang", "redundancy")

  if (method == "all") {
    cc      <- compare_clustering(x, normalize = normalize)
    cc_cols <- intersect(valid, names(cc))
  } else {
    if (!method %in% valid) {
      stop("`method` must be one of: ", paste(c(valid, "all"), collapse = ", "))
    }
    cc      <- clustering_coef(x, method = method, normalize = normalize)
    cc_cols <- method
  }

  if (engine == "base") {
    m <- if (method == "all") "barrat" else method
    vals <- cc[[m]]
    ord  <- order(vals, na.last = TRUE)
    graphics::par(mar = c(4, 8, 3, 2))
    graphics::barplot(
      vals[ord],
      names.arg = cc$variable[ord],
      horiz     = TRUE,
      las       = 1,
      col       = "#D85A30",
      border    = NA,
      xlab      = if (normalize) paste(m, "(normalised)") else m,
      main      = if (!is.null(title)) title else paste("Clustering coefficient:", m)
    )
    return(invisible(NULL))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  pal <- c(
    watts_strogatz = "#888780",
    barrat         = "#D85A30",
    onnela         = "#7F77DD",
    zhang          = "#1D9E75",
    redundancy     = "#BA7517"
  )

  if (method != "all") {
    cc$variable <- factor(cc$variable,
                           levels = cc$variable[order(cc[[method]])])
    p <- ggplot2::ggplot(cc,
           ggplot2::aes(x = .data[[method]], y = .data$variable)) +
      ggplot2::geom_col(fill = pal[method], width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = formatC(.data[[method]], digits = 3, format = "f")),
        hjust = -0.1, size = 3, color = "#444441", na.rm = TRUE
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.15))
      ) +
      ggplot2::labs(
        x     = if (normalize) paste(method, "(normalised)") else method,
        y     = NULL,
        title = title
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  } else {
    # Long format for faceted comparison
    long <- stats::reshape(
      cc[, c("variable", cc_cols)],
      varying   = cc_cols,
      v.names   = "value",
      timevar   = "method",
      times     = cc_cols,
      direction = "long"
    )
    long$method   <- factor(long$method,   levels = cc_cols)
    long$variable <- factor(long$variable,
                             levels = cc$variable[order(cc$mean_cc)])
    long$bar_col  <- pal[as.character(long$method)]

    p <- ggplot2::ggplot(long,
           ggplot2::aes(x = .data$value, y = .data$variable,
                        fill = .data$method)) +
      ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::facet_wrap(~ .data$method, scales = "free_x", nrow = 2) +
      ggplot2::labs(
        x     = if (normalize) "Normalised value" else "Value",
        y     = NULL,
        title = title,
        fill  = "Method"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position    = "none"
      )
  }

  p
}
