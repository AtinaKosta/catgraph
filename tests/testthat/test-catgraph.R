library(testthat)
library(catgraph)

# ------------------------------------------------------------------ helpers
make_df <- function(n = 200, seed = 1) {
  set.seed(seed)
  data.frame(
    A = sample(c("x", "y"),       n, replace = TRUE),
    B = sample(c("yes", "no"),    n, replace = TRUE),
    C = sample(c("lo", "mi", "hi"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------- detect_type tests
test_that("detect_type returns '2x2' for binary pairs", {
  x <- c("a", "b", "a", "b")
  y <- c("1", "0", "1", "0")
  expect_equal(detect_type(x, y), "2x2")
})

test_that("detect_type returns 'RxC' for multinomial pairs", {
  x <- c("a", "b", "c", "a")
  y <- c("1", "0", "1", "0")
  expect_equal(detect_type(x, y), "RxC")
})

test_that("detect_type handles NAs via pairwise deletion", {
  x <- c("a", "b", NA, "a")
  y <- c("1", "0",  "1", "1")
  # After pairwise deletion, x has levels a, b — still 2x2
  expect_equal(detect_type(x, y), "2x2")
})

# ------------------------------------------------------ compute_assoc tests
test_that("compute_assoc returns a list with expected names", {
  df <- make_df()
  res <- compute_assoc(df$A, df$B)
  expect_named(res, c("statistic", "p_value", "df", "n", "table", "type"))
})

test_that("compute_assoc statistic is non-negative", {
  df <- make_df()
  res <- compute_assoc(df$A, df$C)
  expect_gte(res$statistic, 0)
})

test_that("compute_assoc returns NA for all-NA input", {
  x <- rep(NA_character_, 5)
  y <- rep(NA_character_, 5)
  suppressWarnings(res <- compute_assoc(x, y))
  expect_true(is.na(res$statistic))
})

# ------------------------------------------------------- effect_size tests
test_that("phi is in [0, 1] for 2x2 tables", {
  df <- make_df()
  res <- effect_size(df$A, df$B)
  expect_equal(res$metric, "phi")
  expect_gte(res$effect_size, 0)
  expect_lte(res$effect_size, 1)
})

test_that("cramers_v is in [0, 1] for nxm tables", {
  df <- make_df()
  res <- effect_size(df$A, df$C)
  expect_equal(res$metric, "cramers_v")
  expect_gte(res$effect_size, 0)
  expect_lte(res$effect_size, 1)
})

test_that("bias-corrected V is >= 0", {
  df <- make_df()
  res <- effect_size(df$A, df$C, corrected = TRUE)
  expect_gte(res$effect_size, 0)
  expect_true(res$corrected)
})

test_that("corrected phi <= classical phi (on average)", {
  # With n=200 independent draws, corrected should be <= classical
  df <- make_df()
  cl <- effect_size(df$A, df$B, corrected = FALSE)$effect_size
  co <- effect_size(df$A, df$B, corrected = TRUE)$effect_size
  expect_lte(co, cl + 1e-10)   # allow floating point tolerance
})

# --------------------------------------------------------- build_graph tests
test_that("build_graph returns an igraph object", {
  df <- make_df()
  g  <- build_graph(df)
  expect_true(igraph::is_igraph(g))
})

test_that("build_graph has correct edge count (C(p,2))", {
  df <- make_df()  # 3 columns
  g  <- build_graph(df)
  expect_equal(igraph::ecount(g), 3L)  # C(3,2)=3
})

test_that("build_graph edge weights are non-negative", {
  df <- make_df()
  g  <- build_graph(df)
  expect_true(all(igraph::E(g)$weight >= 0, na.rm = TRUE))
})

# ----------------------------------------------------------- catgraph tests
test_that("catgraph() returns a catgraph object", {
  df <- make_df()
  cg <- catgraph(df)
  expect_s3_class(cg, "catgraph")
})

test_that("catgraph slots are consistent", {
  df <- make_df()
  cg <- catgraph(df)
  expect_equal(cg$n_vars, 3L)
  expect_equal(cg$n_pairs, 3L)
  expect_false(cg$corrected)
})

test_that("print.catgraph runs without error", {
  df <- make_df()
  cg <- catgraph(df)
  expect_output(print(cg), "catgraph")
})

test_that("summary.catgraph returns a data frame", {
  df <- make_df()
  cg <- catgraph(df)
  out <- summary(cg)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
})

# --------------------------------------------------------- prune_edges tests
test_that("prune_edges removes edges below threshold", {
  df <- make_df()
  cg <- catgraph(df)
  # Use a threshold above all expected effect sizes for independent data
  cg_p <- prune_edges(cg, min_weight = 0.99)
  expect_lte(cg_p$n_pairs, cg$n_pairs)
})

test_that("prune_edges with max_p = 0 removes all edges", {
  df <- make_df()
  cg <- catgraph(df)
  cg_p <- prune_edges(cg, max_p = 0)
  expect_equal(cg_p$n_pairs, 0L)
})

# ------------------------------------------------------ assoc_matrix tests
test_that("assoc_matrix returns a symmetric matrix", {
  df <- make_df()
  cg <- catgraph(df)
  m  <- assoc_matrix(cg)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(3L, 3L))
  # Symmetry
  expect_equal(m["A", "C"], m["C", "A"])
})

test_that("assoc_matrix tidy format returns a data frame", {
  df <- make_df()
  cg <- catgraph(df)
  td <- assoc_matrix(cg, format = "tidy")
  expect_s3_class(td, "data.frame")
  expect_true("effect_size" %in% names(td))
})

# --------------------------------------------------- detect_clusters tests
test_that("detect_clusters adds clustering slot", {
  df <- make_df()
  cg <- catgraph(df)
  set.seed(1)
  cg <- detect_clusters(cg)
  expect_true(!is.null(cg$clustering))
  expect_true("n_clusters" %in% names(cg$clustering))
})

test_that("vertex cluster attribute is set", {
  df <- make_df()
  cg <- catgraph(df)
  set.seed(1)
  cg <- detect_clusters(cg)
  expect_true(!is.null(igraph::V(cg$graph)$cluster))
})

test_that("print.catgraph reports metric mix", {
  # Mixed 2x2 and nxm edges should produce a "Metric mix" line
  set.seed(1)
  df <- data.frame(
    A = sample(c("x", "y"), 200, replace = TRUE),
    B = sample(c("yes", "no"), 200, replace = TRUE),
    C = sample(c("lo", "mi", "hi"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )
  cg <- suppressWarnings(catgraph(df))
  out <- capture.output(print(cg))
  expect_true(any(grepl("Metric mix", out)))
  expect_true(any(grepl("phi|cramers_v", out)))
})
