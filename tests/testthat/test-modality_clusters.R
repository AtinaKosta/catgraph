library(testthat)
library(catgraph)

test_that("cluster_modalities adds membership and vertex cluster attribute", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg)
  
  expect_s3_class(mg, "catmodgraph")
  expect_false(is.null(mg$membership))
  expect_equal(length(mg$membership), igraph::vcount(mg$graph))
  expect_true("cluster" %in% igraph::vertex_attr_names(mg$graph))
  expect_equal(length(igraph::V(mg$graph)$cluster), igraph::vcount(mg$graph))
})

test_that("cluster_modalities produces at least two clusters on Titanic", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg)
  
  expect_gte(length(unique(mg$membership)), 2)
})

test_that("cluster_modalities errors on non-catmodgraph input", {
  expect_error(
    cluster_modalities(list()),
    "`x` must be a catmodgraph object."
  )
})

test_that("cluster_modalities works with walktrap method", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg, method = "walktrap")
  
  expect_false(is.null(mg$membership))
  expect_gte(length(unique(mg$membership)), 1)
})