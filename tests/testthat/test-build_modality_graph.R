library(testthat)
library(catgraph)

test_that("build_modality_graph returns a catmodgraph object", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_s3_class(mg, "catmodgraph")
  expect_true(is.list(mg))
  expect_true("graph" %in% names(mg))
  expect_true("modalities" %in% names(mg))
  expect_true("indicator_matrix" %in% names(mg))
  expect_true("data" %in% names(mg))
})

test_that("build_modality_graph creates the correct number of modality nodes", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_equal(nrow(mg$modalities), 10)
  expect_equal(igraph::vcount(mg$graph), 10)
  
  expect_equal(
    as.integer(table(mg$modalities$variable)),
    c(2L, 4L, 2L, 2L)
  )
  
  expect_equal(
    names(table(mg$modalities$variable)),
    c("Age", "Class", "Sex", "Survived")
  )
})

test_that("build_modality_graph excludes same-variable modality pairs", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  ed <- igraph::as_data_frame(mg$graph, what = "edges")
  verts <- igraph::as_data_frame(mg$graph, what = "vertices")
  
  var_map <- stats::setNames(verts$variable, verts$name)
  
  edge_var1 <- unname(var_map[ed$from])
  edge_var2 <- unname(var_map[ed$to])
  
  expect_true(all(edge_var1 != edge_var2))
})

test_that("build_modality_graph creates the expected number of edges for Titanic", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  # Titanic modalities:
  # Class = 4, Sex = 2, Age = 2, Survived = 2
  # Total node pairs = choose(10, 2) = 45
  # Same-variable pairs excluded = choose(4,2)+choose(2,2)+choose(2,2)+choose(2,2) = 9
  # Expected cross-variable pairs = 36
  expect_equal(igraph::ecount(mg$graph), 36)
})

test_that("build_modality_graph edge attributes are present and valid", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  w <- igraph::E(mg$graph)$weight
  p <- igraph::E(mg$graph)$p_value
  r <- igraph::E(mg$graph)$std_resid
  
  expect_true(is.numeric(w))
  expect_true(is.numeric(p))
  expect_true(is.numeric(r))
  
  expect_true(all(!is.na(w)))
  expect_true(all(w >= 0))
  expect_true(all(w <= 1))
  
  expect_true(all(is.na(p) | (p >= 0 & p <= 1)))
})

test_that("build_modality_graph indicator matrix dimensions match processed data and modalities", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_equal(nrow(mg$indicator_matrix), nrow(mg$data))
  expect_equal(ncol(mg$indicator_matrix), nrow(mg$modalities))
})

test_that("build_modality_graph errors on non-data-frame input", {
  expect_error(
    build_modality_graph(list(a = 1:3, b = 4:6)),
    "`data` must be a data.frame."
  )
})

test_that("build_modality_graph errors when fewer than two variables are provided", {
  df <- data.frame(A = c("x", "y", "x"), stringsAsFactors = FALSE)
  
  expect_error(
    build_modality_graph(df),
    "`data` must contain at least two variables."
  )
})

test_that("build_modality_graph removes incomplete rows when remove_na = TRUE", {
  df <- data.frame(
    A = c("x", "y", NA, "x"),
    B = c("u", "v", "u", NA),
    C = c("m", "m", "n", "n"),
    stringsAsFactors = FALSE
  )
  
  mg <- build_modality_graph(df, remove_na = TRUE)
  
  expect_equal(nrow(mg$data), 2)
  expect_equal(nrow(mg$indicator_matrix), 2)
})

test_that("build_modality_graph keeps complete structure with min_count filtering", {
  df <- data.frame(
    A = c("x", "x", "x", "y"),
    B = c("u", "u", "v", "v"),
    C = c("m", "n", "m", "n"),
    stringsAsFactors = FALSE
  )
  
  mg <- build_modality_graph(df, min_count = 2)
  
  expect_true(all(colSums(mg$indicator_matrix) >= 2))
  expect_equal(ncol(mg$indicator_matrix), nrow(mg$modalities))
})