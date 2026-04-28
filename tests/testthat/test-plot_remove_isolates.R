library(testthat)
library(catgraph)

test_that("plot.catmodgraph removes isolates by default", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  # Heavy pruning to create isolates
  mg <- prune_modality_edges(mg, min_weight = 0.20)
  
  # Should not error even though isolates exist
  expect_error(plot(mg), NA)
})

test_that("plot.catmodgraph respects remove_isolates = FALSE", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- prune_modality_edges(mg, min_weight = 0.20)
  
  expect_error(plot(mg, remove_isolates = FALSE), NA)
})

test_that("plot.catmodgraph validates remove_isolates argument", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(plot(mg, remove_isolates = "yes"),
               "`remove_isolates` must be TRUE or FALSE.")
  expect_error(plot(mg, remove_isolates = NA),
               "`remove_isolates` must be TRUE or FALSE.")
})

test_that("plot.catmodgraph does not modify the input catmodgraph", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- prune_modality_edges(mg, min_weight = 0.20)
  
  n_before <- igraph::vcount(mg$graph)
  invisible(plot(mg, remove_isolates = TRUE))
  n_after <- igraph::vcount(mg$graph)
  
  expect_equal(n_before, n_after)
})

test_that("plot.catmodgraph errors gracefully on fully-isolated graph", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  # Extreme pruning: remove all edges, leaving only isolates
  mg <- prune_modality_edges(mg, min_weight = 0.999)
  
  expect_error(
    plot(mg, remove_isolates = TRUE),
    "All vertices are isolates"
  )
})