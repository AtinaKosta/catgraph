library(testthat)
library(catgraph)

# Helper: two small datasets drawn from the same distribution
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

test_that("test_modality_graph_equality runs with defaults", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  res <- test_modality_graph_equality(
    mg_x, mg_y, n_perm = 100, seed = 1, verbose = FALSE
  )
  
  expect_s3_class(res, "catmodtest")
  expect_true(is.numeric(res$observed))
  expect_true(is.numeric(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_equal(length(res$null_distribution), 100)
})

test_that("test_modality_graph_equality runs with each statistic", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  for (stat in c("frobenius", "jaccard", "max")) {
    res <- test_modality_graph_equality(
      mg_x, mg_y, n_perm = 50, statistic = stat,
      seed = 1, verbose = FALSE
    )
    expect_s3_class(res, "catmodtest")
    expect_equal(res$statistic, stat)
  }
})

test_that("test_modality_graph_equality runs in pipeline mode", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  res <- test_modality_graph_equality(
    mg_x, mg_y, n_perm = 50, test_type = "pipeline",
    seed = 1, verbose = FALSE
  )
  expect_s3_class(res, "catmodtest")
  expect_equal(res$test_type, "pipeline")
})

test_that("test_modality_graph_equality validates inputs", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  expect_error(test_modality_graph_equality(list(), mg_y),
               "must be a catmodgraph")
  expect_error(test_modality_graph_equality(mg_x, list()),
               "must be a catmodgraph")
  expect_error(test_modality_graph_equality(mg_x, mg_y, n_perm = 5),
               "must be a single integer >= 10")
  expect_error(test_modality_graph_equality(mg_x, mg_y,
                                            statistic = "weird"),
               "should be one of")
})

test_that("test_modality_graph_equality enforces shared variable set", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  # Build y with a different variable set
  d_y2 <- d$y
  colnames(d_y2) <- c("a", "b", "z")
  mg_y2 <- build_modality_graph(d_y2)
  
  expect_error(
    test_modality_graph_equality(mg_x, mg_y2, n_perm = 50,
                                 verbose = FALSE),
    "same variables"
  )
})

test_that("stratified permutation accepts strata argument", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  strata <- sample(c("A", "B"), nrow(d$x) + nrow(d$y), replace = TRUE)
  
  res <- test_modality_graph_equality(
    mg_x, mg_y, n_perm = 50, strata = strata,
    seed = 1, verbose = FALSE
  )
  expect_true(res$strata_used)
})

test_that("strata length mismatch errors", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  expect_error(
    test_modality_graph_equality(mg_x, mg_y, n_perm = 50,
                                 strata = c("A", "B"),
                                 verbose = FALSE),
    "must have length"
  )
})

test_that("print, summary, plot methods work", {
  d <- make_paired_data()
  mg_x <- build_modality_graph(d$x)
  mg_y <- build_modality_graph(d$y)
  
  res <- test_modality_graph_equality(
    mg_x, mg_y, n_perm = 50, seed = 1, verbose = FALSE
  )
  
  expect_output(print(res), "Permutation test")
  s <- summary(res)
  expect_true(is.list(s))
  expect_true("decision" %in% names(s))
  expect_error(plot(res), NA)
})