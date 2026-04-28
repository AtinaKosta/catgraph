library(testthat)
library(catgraph)

# ------------------------------------------------------------------ helpers
make_associated <- function(n = 150, seed = 42) {
  set.seed(seed)
  x <- sample(c("A", "B"), n, replace = TRUE, prob = c(0.4, 0.6))
  y <- ifelse(
    x == "A",
    sample(c("yes", "no"), n, replace = TRUE, prob = c(0.75, 0.25)),
    sample(c("yes", "no"), n, replace = TRUE, prob = c(0.25, 0.75))
  )
  list(x = x, y = y)
}

make_df <- function(n = 200, seed = 1) {
  set.seed(seed)
  data.frame(
    A = sample(c("x", "y"), n, replace = TRUE),
    B = sample(c("yes", "no"), n, replace = TRUE),
    C = sample(c("lo", "mi", "hi"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# Wrapper: tests use R < 500 for speed, which legitimately triggers the
# "R < 500 may yield unstable confidence intervals" warning. That warning
# is useful for end users but noise in a test suite, so we silence it here.
boot_q <- function(...) suppressWarnings(bootstrap_ci(...))
cgci_q <- function(...) suppressWarnings(catgraph_ci(...))

# --------------------------------------------- bootstrap_ci: basic output
test_that("bootstrap_ci returns list with expected names", {
  d  <- make_associated()
  ci <- boot_q(d$x, d$y, R = 200, seed = 1)
  expect_named(
    ci,
    c("estimate", "ci_lower", "ci_upper", "conf", "type",
      "R", "metric", "corrected", "n", "boot_dist")
  )
})

test_that("bootstrap_ci lower <= estimate <= upper", {
  d  <- make_associated()
  ci <- boot_q(d$x, d$y, R = 200, seed = 1)
  expect_lte(ci$ci_lower, ci$estimate + 1e-10)
  expect_gte(ci$ci_upper, ci$estimate - 1e-10)
})

test_that("bootstrap_ci bounds are in [0, 1]", {
  d  <- make_associated()
  ci <- boot_q(d$x, d$y, R = 200, seed = 1)
  expect_gte(ci$ci_lower, 0)
  expect_lte(ci$ci_upper, 1)
})

test_that("bootstrap_ci conf level is preserved", {
  d  <- make_associated()
  ci <- boot_q(d$x, d$y, R = 200, conf = 0.90, seed = 1)
  expect_equal(ci$conf, 0.90)
})

test_that("bootstrap_ci BCa type runs without error", {
  d  <- make_associated()
  ci <- suppressWarnings(
    boot_q(d$x, d$y, R = 200, type = "bca", seed = 1)
  )
  expect_equal(ci$type, "bca")
  expect_true(is.numeric(ci$ci_lower))
})

test_that("bootstrap_ci corrected flag is passed through", {
  d  <- make_associated()
  ci <- boot_q(d$x, d$y, R = 200, corrected = TRUE, seed = 1)
  expect_true(ci$corrected)
})

test_that("bootstrap_ci seed ensures reproducibility", {
  d   <- make_associated()
  ci1 <- boot_q(d$x, d$y, R = 200, seed = 99)
  ci2 <- boot_q(d$x, d$y, R = 200, seed = 99)
  expect_equal(ci1$ci_lower, ci2$ci_lower)
  expect_equal(ci1$ci_upper, ci2$ci_upper)
})

test_that("bootstrap_ci errors on too few complete observations", {
  x <- c(NA, NA, "a", NA)
  y <- c("b", NA, "c", NA)
  expect_error(boot_q(x, y, R = 100))
})

test_that("bootstrap_ci boot_dist has approximately R entries", {
  d  <- make_associated()
  ci <- boot_q(d$x, d$y, R = 200, seed = 1)
  # May be slightly fewer than R if degenerate resamples are dropped
  expect_lte(length(ci$boot_dist), 200)
  expect_gte(length(ci$boot_dist), 180)  # allow up to 10% degenerate
})

# --------------------------------------------- bootstrap_ci: Cramer's V
test_that("bootstrap_ci uses cramers_v for nxm pair", {
  set.seed(7)
  x  <- sample(c("A", "B"), 120, replace = TRUE)
  z  <- sample(c("lo", "mi", "hi"), 120, replace = TRUE)
  ci <- boot_q(x, z, R = 200, seed = 1)
  expect_equal(ci$metric, "cramers_v")
})

# --------------------------------------------- catgraph_ci integration
test_that("catgraph_ci adds ci edge attributes", {
  df <- make_df()
  cg <- catgraph(df)
  cg <- cgci_q(cg, R = 100, seed = 1, verbose = FALSE)
  expect_false(is.null(igraph::E(cg$graph)$ci_lower))
  expect_false(is.null(igraph::E(cg$graph)$ci_upper))
  expect_false(is.null(igraph::E(cg$graph)$ci_conf))
  expect_false(is.null(igraph::E(cg$graph)$ci_type))
})

test_that("catgraph_ci ci_lower <= weight <= ci_upper for each edge", {
  df <- make_df()
  cg <- catgraph(df)
  cg <- cgci_q(cg, R = 100, seed = 1, verbose = FALSE)
  w  <- igraph::E(cg$graph)$weight
  lo <- igraph::E(cg$graph)$ci_lower
  hi <- igraph::E(cg$graph)$ci_upper
  # Percentile CI should bracket the estimate (with small float tolerance)
  expect_true(all(lo <= w + 1e-10, na.rm = TRUE))
  expect_true(all(hi >= w - 1e-10, na.rm = TRUE))
})
