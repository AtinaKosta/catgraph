library(testthat)
library(catgraph)

test_that("plot.catmodgraph runs with default settings", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_invisible(plot(mg))
})

test_that("plot.catmodgraph runs with cluster colouring after clustering", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  mg <- cluster_modalities(mg)
  
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_invisible(plot(mg, color_by = "cluster"))
})

test_that("plot.catmodgraph errors if cluster colouring requested before clustering", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_error(
    plot(mg, color_by = "cluster"),
    "`color_by = \"cluster\"` requires `cluster_modalities\\(\\)` to be run first."
  )
})

test_that("plot.catmodgraph accepts show_labels = FALSE", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_invisible(plot(mg, show_labels = FALSE))
})

test_that("plot.catmodgraph accepts kk layout", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  
  expect_invisible(plot(mg, layout = "kk"))
})

test_that("plot.catmodgraph errors on non-catmodgraph input", {
  expect_error(
    plot.catmodgraph(list()),
    "`x` must be a catmodgraph object."
  )
})

test_that("plot.catmodgraph validates show_labels", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(
    plot(mg, show_labels = "yes"),
    "`show_labels` must be TRUE or FALSE."
  )
})

test_that("plot.catmodgraph validates vertex_size", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(
    plot(mg, vertex_size = 0),
    "`vertex_size` must be a single positive number."
  )
})

test_that("plot.catmodgraph validates edge_scale", {
  df <- expand_table(Titanic)
  mg <- build_modality_graph(df)
  
  expect_error(
    plot(mg, edge_scale = 0),
    "`edge_scale` must be a single positive number."
  )
})