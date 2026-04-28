library(testthat)
library(catgraph)

make_two_catgraphs <- function() {
  df <- expand_table(HairEyeColor)
  # Sparse-cell warnings are expected on these small subsets; suppress in
  # the helper so test output stays focused on the behaviour under test.
  suppressWarnings(list(
    Female = catgraph(df[df$Sex == "Female", c("Hair", "Eye")],
                      corrected = TRUE),
    Male   = catgraph(df[df$Sex == "Male",   c("Hair", "Eye")],
                      corrected = TRUE)
  ))
}

test_that("compare_catgraphs runs with default pruning", {
  cgs <- make_two_catgraphs()
  expect_error(compare_catgraphs(cgs), NA)
})

test_that("compare_catgraphs runs with each pruning mode", {
  cgs <- make_two_catgraphs()
  for (mode in c("pooled", "individual", "overlay", "none")) {
    expect_error(compare_catgraphs(cgs, pruning = mode), NA,
                 info = paste("mode:", mode))
  }
})

test_that("compare_catgraphs validates input list", {
  expect_error(compare_catgraphs("not a list"),
               "must be a list of at least two")
  expect_error(compare_catgraphs(list()),
               "at least two")
})

test_that("compare_catgraphs requires named list", {
  cgs <- make_two_catgraphs()
  names(cgs) <- NULL
  expect_error(compare_catgraphs(cgs),
               "named list with unique, non-empty names")
})

test_that("compare_catgraphs requires all elements to be catgraphs", {
  cgs <- make_two_catgraphs()
  cgs$Male <- list()
  expect_error(compare_catgraphs(cgs),
               "must be catgraph objects")
})

test_that("compare_catgraphs requires matching variable sets", {
  df <- expand_table(HairEyeColor)
  suppressWarnings({
    cg1 <- catgraph(df[, c("Hair", "Eye")], corrected = TRUE)
    cg2 <- catgraph(df[, c("Hair", "Sex")], corrected = TRUE)
  })
  expect_error(
    compare_catgraphs(list(a = cg1, b = cg2)),
    "same variable set"
  )
})

test_that("compare_catgraphs returns the union graph invisibly", {
  cgs <- make_two_catgraphs()
  ug <- compare_catgraphs(cgs)
  expect_true(igraph::is_igraph(ug))
})