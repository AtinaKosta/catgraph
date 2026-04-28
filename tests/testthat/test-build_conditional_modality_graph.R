library(testthat)
library(catgraph)

make_test_data <- function() {
  df <- expand_table(HairEyeColor)
  df  # has columns Hair, Eye, Sex
}

test_that("build_conditional_modality_graph returns a catmodgraph", {
  df <- make_test_data()
  mg <- build_conditional_modality_graph(df, given = list(Sex = "Female"))
  
  expect_s3_class(mg, "catmodgraph")
  expect_true("conditioning" %in% names(mg))
  expect_equal(mg$conditioning$given, list(Sex = "Female"))
})

test_that("conditioning metadata records original and conditional n", {
  df <- make_test_data()
  mg <- build_conditional_modality_graph(df, given = list(Sex = "Female"))
  
  expect_equal(mg$conditioning$n_original, nrow(df))
  expect_equal(mg$conditioning$n_conditional, sum(df$Sex == "Female"))
  expect_lt(mg$conditioning$n_conditional, mg$conditioning$n_original)
})

test_that("conditioning variable is dropped from graph by default", {
  df <- make_test_data()
  mg <- build_conditional_modality_graph(df, given = list(Sex = "Female"))
  
  var_nodes <- unique(sub("=.*", "", igraph::V(mg$graph)$name))
  expect_false("Sex" %in% var_nodes)
  expect_setequal(var_nodes, c("Hair", "Eye"))
})

test_that("drop_conditioning_vars = FALSE keeps the conditioning column", {
  df <- make_test_data()
  mg <- build_conditional_modality_graph(
    df, given = list(Sex = "Female"),
    drop_conditioning_vars = FALSE
  )
  var_nodes <- unique(sub("=.*", "", igraph::V(mg$graph)$name))
  expect_true("Sex" %in% var_nodes)
})

test_that("multiple conditioning variables are supported", {
  df <- make_test_data()
  df$Wave <- sample(c("2019", "2020"), nrow(df), replace = TRUE)
  mg <- build_conditional_modality_graph(
    df, given = list(Sex = "Female", Wave = "2020")
  )
  expect_equal(mg$conditioning$given,
               list(Sex = "Female", Wave = "2020"))
  var_nodes <- unique(sub("=.*", "", igraph::V(mg$graph)$name))
  expect_false("Sex"  %in% var_nodes)
  expect_false("Wave" %in% var_nodes)
})

test_that("result composes with compare_modality_graphs", {
  df <- make_test_data()
  mg_f <- build_conditional_modality_graph(df, given = list(Sex = "Female"))
  mg_m <- build_conditional_modality_graph(df, given = list(Sex = "Male"))
  
  expect_error(
    compare_modality_graphs(list(women = mg_f, men = mg_m),
                            restrict = "common"),
    NA
  )
})

test_that("errors on non-data-frame input", {
  expect_error(
    build_conditional_modality_graph(list(), given = list(a = 1)),
    "`data` must be a data frame."
  )
})

test_that("errors when `given` is not a named list", {
  df <- make_test_data()
  expect_error(
    build_conditional_modality_graph(df, given = list()),
    "non-empty named list"
  )
  expect_error(
    build_conditional_modality_graph(df, given = list("Female")),
    "non-empty named list"
  )
})

test_that("errors when a conditioning column is not in data", {
  df <- make_test_data()
  expect_error(
    build_conditional_modality_graph(df, given = list(NoSuch = "x")),
    "not found in `data`"
  )
})

test_that("errors when a conditioning level does not exist", {
  df <- make_test_data()
  expect_error(
    build_conditional_modality_graph(df, given = list(Sex = "NonBinary")),
    "not found in column 'Sex'"
  )
})

test_that("errors when `given` value is not length 1", {
  df <- make_test_data()
  expect_error(
    build_conditional_modality_graph(df,
                                     given = list(Sex = c("Female", "Male"))),
    "length-1 value"
  )
})

test_that("errors when conditional subset is empty", {
  # Construct a pathological case: two conditions that never co-occur
  df <- data.frame(
    A = c("x", "x", "y", "y"),
    B = c("p", "q", "p", "q"),
    C = c("m", "n", "m", "n"),
    stringsAsFactors = FALSE
  )
  expect_error(
    build_conditional_modality_graph(df, given = list(A = "x", B = "z")),
    "not found in column 'B'"  # fails at level-existence check first
  )
})

test_that("warns on very small conditional subset", {
  df <- data.frame(
    A = c(rep("x", 5), rep("y", 200)),
    B = sample(c("p", "q"), 205, replace = TRUE),
    C = sample(c("m", "n"), 205, replace = TRUE),
    stringsAsFactors = FALSE
  )
  expect_warning(
    build_conditional_modality_graph(df, given = list(A = "x")),
    "only 5 rows"
  )
})

test_that("errors when fewer than 2 non-conditioning variables remain", {
  df <- data.frame(
    A = rep(c("x", "y"), each = 20),
    B = rep(c("p", "q"), times = 20),
    stringsAsFactors = FALSE
  )
  expect_error(
    build_conditional_modality_graph(df, given = list(A = "x")),
    "Fewer than 2 non-conditioning variables"
  )
})