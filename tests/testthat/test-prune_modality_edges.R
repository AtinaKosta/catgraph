library(testthat)
library(catgraph)

test_that("prune_modality_edges returns a catmodgraph object", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg2 <- prune_modality_edges(mg, min_weight = 0.1)
  
  expect_s3_class(mg2, "catmodgraph")
})

test_that("prune_modality_edges reduces edge count when threshold is applied", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg2 <- prune_modality_edges(mg, min_weight = 0.1)
  
  expect_lte(igraph::ecount(mg2$graph), igraph::ecount(mg$graph))
})

test_that("prune_modality_edges can remove all edges", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg2 <- prune_modality_edges(mg, min_weight = 1.1)
  
  expect_equal(igraph::ecount(mg2$graph), 0)
  expect_equal(igraph::vcount(mg2$graph), igraph::vcount(mg$graph))
})

test_that("prune_modality_edges can remove isolates", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg2 <- prune_modality_edges(mg, min_weight = 1.1, remove_isolates = TRUE)
  
  expect_equal(igraph::ecount(mg2$graph), 0)
  expect_equal(igraph::vcount(mg2$graph), 0)
  expect_equal(nrow(mg2$modalities), 0)
  expect_equal(ncol(mg2$indicator_matrix), 0)
})


test_that("prune_modality_edges keeps membership aligned when present", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg)
  
  mg2 <- prune_modality_edges(mg, min_weight = 0.1)
  
  if (!is.null(mg2$membership)) {
    expect_equal(length(mg2$membership), igraph::vcount(mg2$graph))
    expect_true(all(names(mg2$membership) %in% igraph::V(mg2$graph)$name))
  } else {
    expect_equal(igraph::vcount(mg2$graph), 0)
  }
})

test_that("prune_modality_edges validates input class", {
  expect_error(
    prune_modality_edges(list()),
    "`x` must be a catmodgraph object."
  )
})

test_that("prune_modality_edges validates min_weight", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(
    prune_modality_edges(mg, min_weight = -0.1),
    "`min_weight` must be a single non-negative number."
  )
})

test_that("prune_modality_edges validates max_p", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(
    prune_modality_edges(mg, max_p = 2),
    "`max_p` must be a single number in \\[0, 1\\]."
  )
})

test_that("prune_modality_edges validates remove_isolates", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(
    prune_modality_edges(mg, remove_isolates = "yes"),
    "`remove_isolates` must be TRUE or FALSE."
  )
})