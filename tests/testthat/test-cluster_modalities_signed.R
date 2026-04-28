library(testthat)
library(catgraph)

test_that("cluster_modalities accepts signed = TRUE with default louvain", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(cluster_modalities(mg, signed = TRUE), NA)
})

test_that("cluster_modalities(signed = TRUE) produces membership for every vertex", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg, signed = TRUE)
  
  expect_false(is.null(mg$membership))
  expect_equal(length(mg$membership), igraph::vcount(mg$graph))
  expect_false(any(is.na(mg$membership)))
  expect_true("cluster" %in% igraph::vertex_attr_names(mg$graph))
})

test_that("cluster_modalities(signed = TRUE) preserves original graph edges", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  n_edges_before <- igraph::ecount(mg$graph)
  
  mg <- cluster_modalities(mg, signed = TRUE)
  
  # The original graph should still have all its edges, including repulsion ones
  expect_equal(igraph::ecount(mg$graph), n_edges_before)
})

test_that("cluster_modalities validates signed argument", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(cluster_modalities(mg, signed = "yes"),
               "`signed` must be TRUE or FALSE.")
  expect_error(cluster_modalities(mg, signed = NA),
               "`signed` must be TRUE or FALSE.")
  expect_error(cluster_modalities(mg, signed = c(TRUE, FALSE)),
               "`signed` must be TRUE or FALSE.")
})

test_that("cluster_modalities errors on signed + walktrap combination", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(
    cluster_modalities(mg, method = "walktrap", signed = TRUE),
    "supported only for method = 'louvain'"
  )
})

test_that("cluster_modalities(signed = TRUE) errors when std_resid missing", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg$graph <- igraph::delete_edge_attr(mg$graph, "std_resid")
  
  expect_error(
    cluster_modalities(mg, signed = TRUE),
    "requires the `std_resid` edge attribute"
  )
})

test_that("signed and unsigned clustering run without error", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  set.seed(1)
  mg_u <- cluster_modalities(mg, signed = FALSE)
  set.seed(1)
  mg_s <- cluster_modalities(mg, signed = TRUE)
  
  expect_false(is.null(mg_u$membership))
  expect_false(is.null(mg_s$membership))
  expect_true(is.numeric(mg_u$membership))
  expect_true(is.numeric(mg_s$membership))
})