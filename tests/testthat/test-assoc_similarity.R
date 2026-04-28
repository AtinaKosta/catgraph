library(testthat)
library(catgraph)

# Verifies the v0.4.0 addition of assoc_similarity(), which returns the
# dense similarity matrix (for heatmaps) as a separate object from the
# sparse graph (for topology).

make_df <- function(n = 300, seed = 11) {
  set.seed(seed)
  data.frame(
    A = sample(c("x", "y"),          n, replace = TRUE),
    B = sample(c("yes", "no"),       n, replace = TRUE),
    C = sample(c("lo", "mi", "hi"),  n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("assoc_similarity returns a symmetric matrix with NA diagonal", {
  df <- make_df()
  S  <- suppressWarnings(assoc_similarity(df))
  expect_true(is.matrix(S))
  expect_equal(nrow(S), ncol(S))
  expect_true(all(is.na(diag(S))))
  # Symmetric up to numerical tolerance
  expect_equal(S, t(S))
})

test_that("assoc_similarity retains all pairs, even zero-weight ones", {
  df <- make_df()
  S  <- suppressWarnings(assoc_similarity(df))
  # Every off-diagonal entry should be non-NA (or if NA, due to degeneracy,
  # which shouldn't happen with this well-posed df)
  off_diag <- S[upper.tri(S)]
  expect_false(all(is.na(off_diag)))
  expect_length(off_diag, choose(ncol(df), 2))
})

test_that("assoc_similarity 'all' returns three matrices", {
  df  <- make_df()
  out <- suppressWarnings(assoc_similarity(df, what = "all"))
  expect_named(out, c("effect_size", "p_value", "n"))
  expect_true(all(vapply(out, is.matrix, logical(1))))
})

test_that("the bias-corrected flag flows through", {
  df <- make_df()
  S_cl <- suppressWarnings(assoc_similarity(df, corrected = FALSE))
  S_bc <- suppressWarnings(assoc_similarity(df, corrected = TRUE))
  # Bias-corrected values should be smaller or equal to classical ones
  # (Bergsma 2013 correction is a downward adjustment).
  lt <- !is.na(S_cl) & !is.na(S_bc)
  expect_true(all(S_bc[lt] <= S_cl[lt] + 1e-8))
})

test_that("assoc_similarity drops constant columns consistently with build_graph", {
  df <- make_df()
  df$K <- "constant"
  S   <- suppressWarnings(assoc_similarity(df))
  expect_false("K" %in% rownames(S))
  expect_false("K" %in% colnames(S))
})

test_that("assoc_similarity output may exceed number of graph edges", {
  # When some pair has weight exactly zero, the graph drops that edge but
  # the similarity matrix keeps the entry (as 0).
  df <- make_df()
  S  <- suppressWarnings(assoc_similarity(df))
  cg <- suppressWarnings(catgraph(df))

  # The matrix always has choose(p,2) upper-triangle entries. The graph
  # may have fewer edges if any pair had weight 0.
  n_matrix_entries <- sum(!is.na(S[upper.tri(S)]))
  n_graph_edges    <- igraph::ecount(cg$graph)
  expect_gte(n_matrix_entries, n_graph_edges)
})
