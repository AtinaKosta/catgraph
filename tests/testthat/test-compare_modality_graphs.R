library(testthat)
library(catgraph)

make_two_modgraphs <- function() {
  df <- expand_table(HairEyeColor)
  # Sparse-cell warnings are expected; suppress for clean test output.
  suppressWarnings(list(
    Female = build_modality_graph(df[df$Sex == "Female", c("Hair", "Eye")]),
    Male   = build_modality_graph(df[df$Sex == "Male",   c("Hair", "Eye")])
  ))
}

test_that("compare_modality_graphs runs with defaults", {
  mgs <- make_two_modgraphs()
  expect_error(compare_modality_graphs(mgs), NA)
})

test_that("compare_modality_graphs runs with restrict = 'union'", {
  mgs <- make_two_modgraphs()
  expect_error(compare_modality_graphs(mgs, restrict = "union"), NA)
})

test_that("compare_modality_graphs runs with pruning = 'none'", {
  mgs <- make_two_modgraphs()
  expect_error(compare_modality_graphs(mgs, pruning = "none"), NA)
})

test_that("compare_modality_graphs runs with signed = TRUE", {
  mgs <- make_two_modgraphs()
  expect_error(compare_modality_graphs(mgs, signed = TRUE), NA)
})

test_that("compare_modality_graphs validates input list", {
  expect_error(compare_modality_graphs("not a list"),
               "must be a list of at least two")
  expect_error(compare_modality_graphs(list()),
               "at least two")
})

test_that("compare_modality_graphs requires named list", {
  mgs <- make_two_modgraphs()
  names(mgs) <- NULL
  expect_error(compare_modality_graphs(mgs),
               "named list with unique, non-empty names")
})

test_that("compare_modality_graphs validates signed", {
  mgs <- make_two_modgraphs()
  expect_error(compare_modality_graphs(mgs, signed = "yes"),
               "`signed` must be TRUE or FALSE.")
})

test_that("compare_modality_graphs returns the union graph invisibly", {
  mgs <- make_two_modgraphs()
  ug <- compare_modality_graphs(mgs)
  expect_true(igraph::is_igraph(ug))
})