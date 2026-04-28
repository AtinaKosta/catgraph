library(testthat)
library(catgraph)

make_paired_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  mkdf <- function() {
    data.frame(
      a = sample(letters[1:3], n, replace = TRUE),
      b = sample(letters[1:3], n, replace = TRUE),
      c = sample(letters[1:3], n, replace = TRUE),
      stringsAsFactors = TRUE
    )
  }
  list(x = mkdf(), y = mkdf())
}

test_that("test_modality_edge_differences runs with defaults", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  res <- test_modality_edge_differences(
    mg_x, mg_y, n_perm = 100, seed = 1, verbose = FALSE
  )
  
  expect_s3_class(res, "catmodedgetest")
  expect_true(is.data.frame(res$edge_table))
  expect_true(all(c("from", "to", "weight_x", "weight_y",
                    "obs_diff", "p_empirical",
                    "p_adjusted") %in% names(res$edge_table)))
})

test_that("test_modality_edge_differences runs with edges = 'union'", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  res <- test_modality_edge_differences(
    mg_x, mg_y, n_perm = 50, edges = "union",
    seed = 1, verbose = FALSE
  )
  expect_equal(res$edges, "union")
})

test_that("test_modality_edge_differences p-values are in [0,1]", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  res <- test_modality_edge_differences(
    mg_x, mg_y, n_perm = 50, seed = 1, verbose = FALSE
  )
  
  expect_true(all(res$edge_table$p_empirical >= 0 &
                    res$edge_table$p_empirical <= 1))
  expect_true(all(res$edge_table$p_adjusted >= 0 &
                    res$edge_table$p_adjusted <= 1))
})

test_that("test_modality_edge_differences validates inputs", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  expect_error(test_modality_edge_differences(list(), mg_y),
               "must be a catmodgraph")
  expect_error(test_modality_edge_differences(mg_x, mg_y, n_perm = 5),
               "must be a single integer >= 10")
})

test_that("print, summary, plot methods work", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  res <- test_modality_edge_differences(
    mg_x, mg_y, n_perm = 50, seed = 1, verbose = FALSE
  )
  
  expect_output(print(res), "Edge-wise permutation test")
  s <- summary(res)
  expect_true(is.list(s))
  expect_error(plot(res), NA)
})