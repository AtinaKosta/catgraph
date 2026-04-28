# tests/testthat/test-catgraph-utils.R
# Merged from test-detect_clusters.R, test-clustering_coef.R,
# and test-compute_assoc-warnings.R (each had only 1 test).

library(testthat)
library(catgraph)

# ---- detect_clusters --------------------------------------------------------

test_that("detect_clusters dispatches correctly for all five methods", {
  set.seed(1)
  df <- data.frame(
    A = sample(c("x", "y"),         250, replace = TRUE),
    B = sample(c("yes", "no"),      250, replace = TRUE),
    C = sample(c("lo", "mi", "hi"), 250, replace = TRUE),
    D = sample(c("red", "blue"),    250, replace = TRUE),
    stringsAsFactors = FALSE
  )
  cg <- suppressWarnings(catgraph(df))
  
  methods <- c("louvain", "walktrap", "fast_greedy",
               "label_prop", "edge_betweenness")
  
  for (m in methods) {
    set.seed(1)
    cg_c <- suppressWarnings(detect_clusters(cg, method = m))
    expect_true(
      !is.null(cg_c$clustering) &&
        inherits(cg_c, "catgraph") &&
        is.numeric(cg_c$clustering$modularity) &&
        cg_c$clustering$n_clusters >= 1L,
      label = paste0("detect_clusters dispatch for method = '", m, "'")
    )
  }
})

# ---- clustering_coef --------------------------------------------------------

test_that("vectorised Zhang-Horvath matches the original triple-loop formula", {
  set.seed(42)
  n <- 250
  df <- data.frame(
    A = sample(c("x", "y"),           n, replace = TRUE),
    B = sample(c("yes", "no"),        n, replace = TRUE),
    C = sample(c("lo", "mi", "hi"),   n, replace = TRUE),
    D = sample(c("red", "blue"),      n, replace = TRUE),
    E = sample(c("p", "q", "r", "s"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  cg <- suppressWarnings(catgraph(df))
  cc <- clustering_coef(cg, method = "zhang", normalize = FALSE)
  
  W <- as.matrix(igraph::as_adjacency_matrix(cg$graph, attr = "weight",
                                             sparse = FALSE))
  W[is.na(W)] <- 0
  diag(W)     <- 0
  A <- (W > 0) * 1.0
  k <- rowSums(A)
  s <- rowSums(W)
  p <- nrow(W)
  
  ref <- numeric(p)
  for (i in seq_len(p)) {
    if (k[i] < 2) { ref[i] <- NA_real_; next }
    numer <- 0
    for (j in seq_len(p)) {
      if (j == i) next
      for (h in seq_len(p)) {
        if (h == i || h == j) next
        numer <- numer + W[i, j] * W[j, h] * W[h, i]
      }
    }
    denom  <- s[i]^2 - sum(W[i, ]^2)
    ref[i] <- if (abs(denom) < .Machine$double.eps) 0 else numer / denom
  }
  
  cc_aligned <- cc$zhang[match(igraph::V(cg$graph)$name, cc$variable)]
  expect_equal(cc_aligned, ref, tolerance = 1e-10)
})

# ---- compute_assoc warnings -------------------------------------------------

test_that("sparse-table warnings include real variable names from build_graph", {
  set.seed(1)
  df <- data.frame(
    USENOW3 = sample(c("1", "2", "3"), 500, replace = TRUE,
                     prob = c(0.49, 0.49, 0.02)),
    ALCDAY4 = sample(as.character(c(101:130, 201:208)), 500, replace = TRUE),
    Z       = sample(c("a", "b"), 500, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  warns <- character()
  withCallingHandlers(
    catgraph(df, corrected = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  expect_true(any(grepl("USENOW3", warns, fixed = TRUE)))
  expect_true(any(grepl("ALCDAY4", warns, fixed = TRUE)))
  expect_false(any(grepl("pair \\(x, y\\)", warns)))
})