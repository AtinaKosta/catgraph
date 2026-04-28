library(testthat)
library(catgraph)

make_df <- function(n = 200, seed = 1) {
  set.seed(seed)
  data.frame(
    A = sample(c("x", "y"), n, replace = TRUE),
    B = sample(c("yes", "no"), n, replace = TRUE),
    C = sample(c("lo", "mi", "hi"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------- plot_heatmap: ggplot2
test_that("plot_heatmap ggplot2 engine returns a ggplot when ggplot2 available", {
  skip_if_not_installed("ggplot2")
  df <- make_df()
  cg <- catgraph(df)
  p  <- plot_heatmap(cg)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap errors on non-catgraph input", {
  expect_error(plot_heatmap(list()), "`x` must be a catgraph object")
})

test_that("plot_heatmap show_values = FALSE still returns ggplot", {
  skip_if_not_installed("ggplot2")
  df <- make_df()
  cg <- catgraph(df)
  p  <- plot_heatmap(cg, show_values = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap show_sig = TRUE runs without error", {
  skip_if_not_installed("ggplot2")
  df <- make_df()
  cg <- catgraph(df)
  expect_silent(plot_heatmap(cg, show_sig = TRUE))
})

test_that("plot_heatmap show_ci warns when no CI present", {
  skip_if_not_installed("ggplot2")
  df <- make_df()
  cg <- catgraph(df)
  expect_warning(plot_heatmap(cg, show_ci = TRUE), "catgraph_ci")
})

test_that("plot_heatmap with CI renders when CI is present", {
  skip_if_not_installed("ggplot2")
  df <- make_df()
  cg <- catgraph(df)
  cg <- suppressWarnings(catgraph_ci(cg, R = 100, seed = 1, verbose = FALSE))
  p  <- plot_heatmap(cg, show_ci = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap reorder = FALSE runs without error", {
  skip_if_not_installed("ggplot2")
  df <- make_df()
  cg <- catgraph(df)
  p  <- plot_heatmap(cg, reorder = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap custom palette is accepted", {
  skip_if_not_installed("ggplot2")
  df <- make_df()
  cg <- catgraph(df)
  p  <- plot_heatmap(cg, palette = c("#FFFFFF", "#1D9E75", "#04342C"))
  expect_s3_class(p, "ggplot")
})

# -------------------------------------------------- plot_heatmap: base
test_that("plot_heatmap base engine runs without error", {
  df  <- make_df()
  cg  <- catgraph(df)
  # base engine uses graphics device; wrap in a null device to avoid window
  grDevices::pdf(tempfile())
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- plot_heatmap(cg, engine = "base")
  expect_null(result)
})


test_that("plot_heatmap uses dense all-pairs matrix, not sparse graph adjacency", {
  skip_if_not_installed("ggplot2")
  
  # Construct a dataset where A and B are exactly independent, so their
  # association is truly zero and should be absent from the graph but present
  # in the dense similarity matrix as 0.
  x <- rep(c("a", "b"), each = 50)
  y <- rep(c("0", "1"), times = 50)
  z <- sample(c("p", "q"), 100, replace = TRUE)
  
  df <- data.frame(A = x, B = y, C = z, stringsAsFactors = FALSE)
  
  cg <- suppressWarnings(catgraph(df))
  
  # The sparse graph should drop the zero A-B edge
  td <- assoc_matrix(cg, format = "tidy")
  expect_false(any((td$var1 == "A" & td$var2 == "B") |
                     (td$var1 == "B" & td$var2 == "A")))
  
  # The dense similarity matrix should retain A-B as 0
  S <- suppressWarnings(assoc_similarity(df, what = "effect_size"))
  expect_equal(S["A", "B"], 0)
  expect_equal(S["B", "A"], 0)
  
  # Heatmap should still render without error from the catgraph object
  p <- plot_heatmap(cg)
  expect_s3_class(p, "ggplot")
})
