library(testthat)
library(catgraph)

make_edge_test <- function() {
  df <- expand_table(HairEyeColor)
  suppressWarnings({
    mg_f <- build_modality_graph(df[df$Sex == "Female", c("Hair", "Eye")])
    mg_m <- build_modality_graph(df[df$Sex == "Male",   c("Hair", "Eye")])
    test_modality_edge_differences(
      mg_f, mg_m, n_perm = 50, edges = "union",
      seed = 1, verbose = FALSE
    )
  })
}

test_that("plot_modality_difference runs with defaults", {
  et <- make_edge_test()
  expect_error(plot_modality_difference(et), NA)
})

test_that("plot_modality_difference returns an igraph invisibly", {
  et <- make_edge_test()
  g <- plot_modality_difference(et)
  expect_true(igraph::is_igraph(g))
  expect_true("obs_diff"   %in% igraph::edge_attr_names(g))
  expect_true("p_adjusted" %in% igraph::edge_attr_names(g))
})

test_that("plot_modality_difference rejects non-catmodedgetest input", {
  expect_error(plot_modality_difference(list()),
               "must be a catmodedgetest object")
})

test_that("plot_modality_difference validates alpha_fdr", {
  et <- make_edge_test()
  expect_error(plot_modality_difference(et, alpha_fdr = 0),
               "`alpha_fdr` must be a single number in \\(0, 1\\]\\.")
  expect_error(plot_modality_difference(et, alpha_fdr = 1.5),
               "`alpha_fdr` must be a single number")
})

test_that("plot_modality_difference validates group_labels", {
  et <- make_edge_test()
  expect_error(plot_modality_difference(et, group_labels = "only one"),
               "length-2 character vector")
})

test_that("plot_modality_difference respects show_nonsig = FALSE", {
  et <- make_edge_test()
  # with n_perm = 50 and HairEyeColor, it's possible nothing is significant;
  # so we test the degenerate message path
  et_fake <- et
  et_fake$edge_table$p_adjusted <- 0.9
  expect_error(
    plot_modality_difference(et_fake, show_nonsig = FALSE),
    "No edges survive alpha_fdr"
  )
})

test_that("plot_modality_difference accepts reference as list of modgraphs", {
  df <- expand_table(HairEyeColor)
  suppressWarnings({
    mg_f <- build_modality_graph(df[df$Sex == "Female", c("Hair", "Eye")])
    mg_m <- build_modality_graph(df[df$Sex == "Male",   c("Hair", "Eye")])
  })
  et <- suppressWarnings(
    test_modality_edge_differences(mg_f, mg_m, n_perm = 50, edges = "union",
                                   seed = 1, verbose = FALSE)
  )
  expect_error(
    plot_modality_difference(et, reference = list(mg_f, mg_m)),
    NA
  )
})