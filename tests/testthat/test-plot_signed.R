library(testthat)
library(catgraph)

test_that("plot.catmodgraph accepts signed = TRUE", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(plot(mg, signed = TRUE), NA)
})

test_that("plot.catmodgraph validates signed argument", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(plot(mg, signed = "yes"),
               "`signed` must be TRUE or FALSE.")
  expect_error(plot(mg, signed = NA),
               "`signed` must be TRUE or FALSE.")
  expect_error(plot(mg, signed = c(TRUE, FALSE)),
               "`signed` must be TRUE or FALSE.")
})

test_that("plot.catmodgraph signed = TRUE with missing std_resid warns and falls back", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  # Strip std_resid using delete_edge_attr for a proper removal
  mg$graph <- igraph::delete_edge_attr(mg$graph, "std_resid")
  
  # Confirm the attribute is really gone before testing the plot
  expect_false("std_resid" %in% igraph::edge_attr_names(mg$graph))
  
  expect_warning(
    plot(mg, signed = TRUE),
    "`std_resid` edge attribute is missing"
  )
})

test_that("plot.catmodgraph signed = FALSE is unchanged (backward compatibility)", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(plot(mg, signed = FALSE), NA)
  expect_error(plot(mg), NA)  # default behavior
})