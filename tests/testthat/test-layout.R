library(testthat)
library(catgraph)

# Verifies the v0.4.0 fix: plot.catgraph() now respects the `layout` argument
# in the igraph branch (it was previously hard-coded to Fruchterman-Reingold).

make_small_cg <- function(seed = 3) {
  set.seed(seed)
  df <- data.frame(
    A = sample(c("x", "y"),          200, replace = TRUE),
    B = sample(c("yes", "no"),       200, replace = TRUE),
    C = sample(c("lo", "mi", "hi"),  200, replace = TRUE),
    D = sample(c("p", "q"),          200, replace = TRUE),
    stringsAsFactors = FALSE
  )
  suppressWarnings(catgraph(df))
}

test_that("plot.catgraph runs for every advertised layout", {
  cg <- make_small_cg()
  # Invisible NULL opens a null device so we don't litter PDF files.
  op <- grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (lay in c("fr", "kk", "circle", "grid", "graphopt", "nicely", "random")) {
    expect_silent(plot(cg, layout = lay))
  }
})

test_that("unknown layout string triggers an informative error", {
  cg <- make_small_cg()
  op <- grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(plot(cg, layout = "does_not_exist"),
               regexp = "Unknown `layout`")
})

test_that("internal layout dispatcher returns a matrix with 2 columns", {
  # .catgraph_igraph_layout is un-exported; use ::: equivalent via getFromNamespace
  f <- getFromNamespace(".catgraph_igraph_layout", "catgraph")
  cg <- make_small_cg()

  for (lay in c("fr", "kk", "circle")) {
    m <- f(cg$graph, lay)
    expect_true(is.matrix(m))
    expect_equal(ncol(m), 2L)
    expect_equal(nrow(m), igraph::vcount(cg$graph))
  }
})

test_that("circle layout differs from FR layout (proving layout arg is used)", {
  f <- getFromNamespace(".catgraph_igraph_layout", "catgraph")
  cg <- make_small_cg()

  set.seed(1); m_fr     <- f(cg$graph, "fr")
  set.seed(1); m_circle <- f(cg$graph, "circle")

  # Coordinates from a circular layout lie on a circle; FR does not. If the
  # layout argument were silently ignored both would be identical.
  expect_false(isTRUE(all.equal(m_fr, m_circle)))
})
