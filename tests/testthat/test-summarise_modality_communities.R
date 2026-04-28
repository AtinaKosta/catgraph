library(testthat)
library(catgraph)

test_that("summarise_modality_communities returns a catmodcommunity object", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg)
  comm <- summarise_modality_communities(mg)
  
  expect_s3_class(comm, "catmodcommunity")
  expect_true(all(c("community_summary", "community_members",
                    "variable_composition") %in% names(comm)))
})

test_that("summarise_modality_communities errors if clustering not run", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(
    summarise_modality_communities(mg),
    "Run `cluster_modalities\\(\\)` before summarising"
  )
})

test_that("community_summary has expected columns", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg)
  comm <- summarise_modality_communities(mg)
  
  expect_true(all(c("community", "n_modalities", "n_variables",
                    "variables", "mean_internal_weight") %in%
                    names(comm$community_summary)))
})

test_that("print.catmodcommunity runs without error", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg)
  comm <- summarise_modality_communities(mg)
  
  expect_output(print(comm), "catmodcommunity")
})

test_that("summary.catmodcommunity returns the community summary data frame", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg)
  comm <- summarise_modality_communities(mg)
  
  expect_s3_class(summary(comm), "data.frame")
})