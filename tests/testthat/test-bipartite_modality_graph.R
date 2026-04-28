library(testthat)
library(catgraph)

test_that("bipartite_modality_graph returns a catbipartite object", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  
  expect_s3_class(bg, "catbipartite")
  expect_true(all(c("graph", "modalities", "n_rows",
                    "n_modalities", "data") %in% names(bg)))
})

test_that("bipartite graph has correct vertex count", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  
  # Total vertices = respondents + modalities
  expect_equal(igraph::vcount(bg$graph), bg$n_rows + bg$n_modalities)
  
  # `type` vertex attribute follows igraph's bipartite convention
  vtype <- igraph::V(bg$graph)$type
  expect_equal(sum(!vtype), bg$n_rows)
  expect_equal(sum(vtype), bg$n_modalities)
})

test_that("every respondent has exactly one edge per variable", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  
  # Each respondent endorses exactly one modality per variable,
  # so respondent degree = number of variables.
  vtype <- igraph::V(bg$graph)$type
  resp_degrees <- igraph::degree(bg$graph)[!vtype]
  
  expect_true(all(resp_degrees == ncol(bg$data)))
})

test_that("bipartite_modality_graph errors on non-data.frame", {
  expect_error(
    bipartite_modality_graph(list(a = 1:3)),
    "`data` must be a data frame."
  )
})

test_that("bipartite_modality_graph errors on single-column data", {
  df <- data.frame(a = c("x", "y", "z"))
  expect_error(
    bipartite_modality_graph(df),
    "at least two variables"
  )
})

test_that("bipartite_modality_graph validates remove_na", {
  df <- expand_table(Titanic)
  expect_error(bipartite_modality_graph(df, remove_na = "yes"),
               "`remove_na` must be TRUE or FALSE.")
})

test_that("bipartite_modality_graph validates row_prefix", {
  df <- expand_table(Titanic)
  expect_error(bipartite_modality_graph(df, row_prefix = ""),
               "`row_prefix` must be a single non-empty character string.")
})

test_that("print.catbipartite runs without error", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  expect_output(print(bg), "catbipartite")
})

test_that("summary.catbipartite returns expected fields", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  s <- summary(bg)
  
  expect_true(all(c("n_rows", "n_modalities", "n_edges",
                    "variables", "modalities_per_variable") %in% names(s)))
})

test_that("plot.catbipartite runs without error", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  
  expect_error(plot(bg), NA)
  expect_error(plot(bg, show_respondents = FALSE), NA)
})

test_that("plot.catbipartite respects max_respondents", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  
  # 2201 respondents; with max_respondents = 100 the plot should subsample
  expect_error(plot(bg, max_respondents = 100), NA)
})

test_that("plot.catbipartite validates arguments", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  
  expect_error(plot(bg, show_respondents = "yes"),
               "`show_respondents` must be TRUE or FALSE.")
  expect_error(plot(bg, max_respondents = -5),
               "`max_respondents` must be NULL or a single positive integer.")
})

test_that("projection onto modality side works", {
  df <- expand_table(Titanic)
  bg <- bipartite_modality_graph(df)
  
  # igraph's bipartite_projection should succeed on our `type` attribute
  proj <- igraph::bipartite_projection(bg$graph, which = "true")
  expect_true(igraph::is_igraph(proj))
  expect_equal(igraph::vcount(proj), bg$n_modalities)
})