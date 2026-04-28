library(testthat)
library(catgraph)

# Verifies the v0.4.0 addition of multiple-testing control to prune_edges().

make_random_wide_df <- function(p = 12, n = 250, seed = 7) {
  # Many independent categorical variables -> many simultaneous tests.
  # With p = 12 variables, there are choose(12, 2) = 66 edges.
  set.seed(seed)
  as.data.frame(
    replicate(p, sample(c("a", "b", "c"), n, replace = TRUE),
              simplify = FALSE),
    col.names = paste0("V", seq_len(p))
  )
}

test_that("prune_edges adds p_value_adj and p_adjust_method edge attributes", {
  df <- make_random_wide_df()
  cg <- suppressWarnings(catgraph(df))
  cp <- prune_edges(cg, p_adjust = "BH")

  expect_true("p_value_adj"     %in% igraph::edge_attr_names(cp$graph))
  expect_true("p_adjust_method" %in% igraph::edge_attr_names(cp$graph))
  expect_true(all(igraph::E(cp$graph)$p_adjust_method == "BH"))
})

test_that("BH adjustment yields fewer surviving edges than raw at same threshold", {
  # All variables are independent -> raw p<0.05 will yield roughly 5% false
  # positives; BH should suppress most of them.
  df <- make_random_wide_df()
  cg <- suppressWarnings(catgraph(df))

  n_raw <- suppressWarnings(
    prune_edges(cg, max_p = 0.05, p_adjust = "none")
  )$n_pairs
  n_bh  <- suppressWarnings(
    prune_edges(cg, max_p = 0.05, p_adjust = "BH")
  )$n_pairs

  # BH should keep no more edges than raw at the same nominal threshold,
  # and on pure noise data usually keeps strictly fewer.
  expect_lte(n_bh, n_raw)
})

test_that("Bonferroni is the most conservative adjustment", {
  df <- make_random_wide_df()
  cg <- suppressWarnings(catgraph(df))

  n_bh   <- prune_edges(cg, max_p = 0.05, p_adjust = "BH")$n_pairs
  n_holm <- prune_edges(cg, max_p = 0.05, p_adjust = "holm")$n_pairs
  n_bonf <- prune_edges(cg, max_p = 0.05, p_adjust = "bonferroni")$n_pairs

  # BH >= Holm is NOT generally guaranteed edge-count-wise, but both
  # should be at least as strict as raw p-values. Bonferroni must be
  # <= Holm.
  expect_lte(n_bonf, n_holm)
})

test_that("p_adjust = 'none' reproduces pre-0.4.0 behaviour", {
  df <- make_random_wide_df()
  cg <- suppressWarnings(catgraph(df))
  cp <- prune_edges(cg, max_p = 0.05, p_adjust = "none")

  # When no adjustment is applied, p_value_adj should equal the raw p-values
  # for all retained edges.
  if (igraph::ecount(cp$graph) > 0L) {
    expect_equal(igraph::E(cp$graph)$p_value_adj,
                 igraph::E(cp$graph)$p_value)
  }
})

test_that("invalid p_adjust string is rejected", {
  df <- make_random_wide_df(p = 4)
  cg <- suppressWarnings(catgraph(df))
  expect_error(prune_edges(cg, p_adjust = "fdr_custom"),
               regexp = "should be one of")
})

test_that("prune_edges handles graphs with no edges gracefully", {
  df <- make_random_wide_df(p = 3)
  cg <- suppressWarnings(catgraph(df))
  cp <- prune_edges(cg, min_weight = 1.1)   # impossible threshold
  # The first prune_edges() above stamped p_value_adj on cp's (now-empty)
  # edge set. A second adjusted prune should warn that adjustment cannot
  # be re-applied on an empty graph.
  expect_warning(
    cp2 <- prune_edges(cp, p_adjust = "BH"),
    regexp = "no edges"
  )
  expect_equal(cp2$n_pairs, 0L)
  
  # With p_adjust = "none", the empty-graph path is silent.
  expect_silent(
    cp3 <- prune_edges(cp, p_adjust = "none")
  )
  expect_equal(cp3$n_pairs, 0L)
})

test_that("prune_edges with remove_isolates drops orphan vertices", {
  df <- make_random_wide_df(p = 6)
  cg <- suppressWarnings(catgraph(df))
  cp <- prune_edges(cg, min_weight = 1.1, remove_isolates = TRUE)
  # All variables should now be removed because every edge was dropped.
  expect_equal(cp$n_vars, 0L)
})


test_that("chained adjusted prunes trigger the re-adjustment warning", {
  set.seed(1)
  df <- as.data.frame(
    replicate(8, sample(c("a", "b", "c"), 250, replace = TRUE),
              simplify = FALSE),
    col.names = paste0("V", 1:8)
  )
  cg <- suppressWarnings(catgraph(df))
  
  # First BH prune with max_p = 1: stamps p_value_adj on every edge but
  # drops nothing, so cg_once is guaranteed non-empty regardless of the
  # random seed or BH denominators.
  cg_once <- suppressWarnings(
    prune_edges(cg, max_p = 1, p_adjust = "BH")
  )
  
  # Sanity: cg_once must have edges AND carry the prior adjustment stamp,
  # otherwise the warning detection in the next call cannot fire and the
  # test would be a false pass.
  stopifnot(igraph::ecount(cg_once$graph) > 0L)
  stopifnot("p_value_adj" %in% igraph::edge_attr_names(cg_once$graph))
  
  # Second BH prune: prior adjustment detected, MUST warn.
  expect_warning(
    prune_edges(cg_once, max_p = 0.05, p_adjust = "BH"),
    regexp = "already been adjusted for multiple testing"
  )
  
  # Explicit 'none' bypasses the warning: user opting out.
  expect_silent(
    suppressMessages(prune_edges(cg_once, max_p = 1, p_adjust = "none"))
  )
})

test_that("single prune_edges call with both thresholds does not warn", {
  # This is the correct usage pattern: one call, both filters.
  set.seed(1)
  df <- as.data.frame(
    replicate(6, sample(c("a", "b", "c"), 250, replace = TRUE),
              simplify = FALSE),
    col.names = paste0("V", 1:6)
  )
  cg <- suppressWarnings(catgraph(df))
  
  expect_silent(
    suppressMessages(
      prune_edges(cg, min_weight = 0.05, max_p = 0.05, p_adjust = "BH")
    )
  )
})