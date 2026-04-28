library(testthat)
library(catgraph)

# Verifies the v0.4.0 fix: build_graph() no longer forces zero weights to
# .Machine$double.eps. True zero associations must be absent edges, not
# near-zero edges.

make_independent_df <- function(n = 400, seed = 42) {
  # Three columns drawn independently — all pairwise associations should be
  # statistically indistinguishable from zero.
  set.seed(seed)
  data.frame(
    A = sample(c("x", "y"),         n, replace = TRUE),
    B = sample(c("yes", "no"),      n, replace = TRUE),
    C = sample(c("lo", "mi", "hi"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

make_degenerate_df <- function(seed = 1) {
  # Constructed so that A vs B has literally zero chi-square (perfectly
  # uniform 2x2 table): we want to see that edge dropped from the graph.
  set.seed(seed)
  x <- rep(c("a", "b"), each = 50)
  y <- rep(c("1", "0"), times = 50)
  z <- sample(c("p", "q"), 100, replace = TRUE)
  data.frame(A = x, B = y, C = z, stringsAsFactors = FALSE)
}

test_that("build_graph does NOT inject epsilon edges for zero weights", {
  df <- make_degenerate_df()
  g  <- build_graph(df)

  w <- igraph::E(g)$weight
  # No edge should sit at machine-epsilon. Either it's a real positive
  # weight or the edge is absent.
  expect_false(any(abs(w - .Machine$double.eps) < 1e-15))
})

test_that("build_graph produces a sparse graph when some pairs have w=0", {
  df <- make_degenerate_df()
  g  <- build_graph(df)

  p <- igraph::vcount(g)
  complete_edges <- p * (p - 1) / 2
  # If the A-B pair has weight exactly 0, we should have strictly fewer
  # edges than the complete graph.
  expect_lt(igraph::ecount(g), complete_edges)
})

test_that("build_graph preserves all vertices even when some are isolated", {
  df <- make_degenerate_df()
  g  <- build_graph(df)
  expect_equal(igraph::vcount(g), ncol(df))
  expect_setequal(igraph::V(g)$name, names(df))
})

test_that("catgraph$data contains the PROCESSED data, not the raw input", {
  # Give build_graph() a numeric column so it must coerce, plus a constant
  # column that must be dropped. $data should reflect both changes.
  df <- data.frame(
    A = sample(c("x", "y"), 100, replace = TRUE),
    B = sample(c("yes", "no"), 100, replace = TRUE),
    # numeric: must be coerced to character
    D = sample(1:3, 100, replace = TRUE),
    # constant: must be dropped
    K = rep("only", 100),
    stringsAsFactors = FALSE
  )
  cg <- suppressMessages(suppressWarnings(catgraph(df)))

  # Processed data should NOT include the dropped constant column.
  expect_false("K" %in% names(cg$data))
  # The coerced D column should now be character, not integer.
  expect_type(cg$data$D, "character")
  # Raw data should still preserve the original structure.
  expect_true("K" %in% names(cg$raw_data))
  expect_type(cg$raw_data$D, "integer")
})

test_that("processed data has exactly the same columns build_graph used", {
  df <- make_independent_df()
  cg <- suppressWarnings(catgraph(df))
  # Every graph vertex corresponds to a column in cg$data
  expect_setequal(names(cg$data), igraph::V(cg$graph)$name)
})

test_that("build_graph attribute is detached in the catgraph object", {
  df <- make_independent_df()
  cg <- suppressWarnings(catgraph(df))
  # The processed_data attribute should not leak into the user-facing graph.
  expect_null(attr(cg$graph, "processed_data"))
})
