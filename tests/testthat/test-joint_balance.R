library(testthat)
library(catgraph)

make_small_data <- function(seed = 1) {
  # A small data frame with a known group difference in the AB association:
  # group "a" has A and B positively associated; group "b" has them independent
  set.seed(seed)
  n  <- 400
  grp <- sample(c("a", "b"), n, replace = TRUE)
  A  <- character(n); B <- character(n); C <- character(n)
  for (i in seq_len(n)) {
    if (grp[i] == "a") {
      A[i] <- sample(c("x", "y"), 1L, prob = c(0.6, 0.4))
      B[i] <- if (A[i] == "x") sample(c("p", "q"), 1L, prob = c(0.8, 0.2))
      else sample(c("p", "q"), 1L, prob = c(0.2, 0.8))
    } else {
      A[i] <- sample(c("x", "y"), 1L, prob = c(0.5, 0.5))
      B[i] <- sample(c("p", "q"), 1L, prob = c(0.5, 0.5))
    }
    C[i] <- sample(c("m", "n"), 1L)
  }
  data.frame(grp = grp, A = A, B = B, C = C, stringsAsFactors = FALSE)
}

test_that("joint_balance returns a jointbalance object with the right shape", {
  df <- make_small_data()
  jb <- joint_balance(df, group = "grp", n_perm = 50, n_perm_edge = 50,
                      seed = 1, verbose = FALSE)
  
  expect_s3_class(jb, "jointbalance")
  expect_named(jb, c("group", "group_levels", "variables",
                     "marginal", "pairwise_omnibus", "pairwise_edgewise",
                     "modality_graphs", "alpha", "call"))
  expect_equal(jb$group, "grp")
  expect_setequal(jb$group_levels, c("a", "b"))
  expect_setequal(jb$variables, c("A", "B", "C"))
})

test_that("marginal table has one row per variable and valid columns", {
  df <- make_small_data()
  jb <- joint_balance(df, group = "grp", n_perm = 50, seed = 1,
                      verbose = FALSE)
  m <- jb$marginal
  expect_equal(nrow(m), 3L)
  expect_setequal(m$variable, c("A", "B", "C"))
  expect_true(all(!is.na(m$p_adjusted)))
  expect_true(all(m$p_adjusted >= 0 & m$p_adjusted <= 1))
})

test_that("pairwise_omnibus has choose(k,2) rows with Bonferroni", {
  df <- make_small_data()
  jb <- joint_balance(df, group = "grp", n_perm = 50, seed = 1,
                      verbose = FALSE)
  po <- jb$pairwise_omnibus
  expect_equal(nrow(po), 1L)  # k = 2 -> one pair
  expect_true(all(po$p_bonferroni >= po$p_value))
  expect_true(all(po$p_bonferroni <= 1))
})

test_that("edge-wise post-hoc runs when omnibus rejects, skipped otherwise", {
  df <- make_small_data()
  jb <- joint_balance(df, group = "grp", n_perm = 200, n_perm_edge = 200,
                      seed = 1, verbose = FALSE)
  if (jb$pairwise_omnibus$p_bonferroni[1L] < jb$alpha) {
    expect_length(jb$pairwise_edgewise, 1L)
    expect_s3_class(jb$pairwise_edgewise[[1L]], "catmodedgetest")
  } else {
    expect_length(jb$pairwise_edgewise, 0L)
  }
})

test_that("run_edgewise = FALSE skips the post-hoc step", {
  df <- make_small_data()
  jb <- joint_balance(df, group = "grp", n_perm = 50,
                      run_edgewise = FALSE, seed = 1, verbose = FALSE)
  expect_length(jb$pairwise_edgewise, 0L)
})

test_that("joint_balance errors on bad group argument", {
  df <- make_small_data()
  expect_error(joint_balance(df, group = "not_a_column"),
               "must be the name of a column")
})

test_that("joint_balance errors when group has < 2 levels", {
  df <- make_small_data()
  df$grp <- "only_one_level"
  expect_error(
    joint_balance(df, group = "grp"),
    "at least 2 levels"
  )
})

test_that("joint_balance warns on numeric variable coercion", {
  df <- make_small_data()
  df$num_var <- rnorm(nrow(df))
  expect_warning(
    joint_balance(df, group = "grp", variables = c("A", "num_var"),
                  n_perm = 50, seed = 1, verbose = FALSE),
    "numeric; coercing to factor"
  )
})

test_that("joint_balance handles k > 2 groups", {
  df <- make_small_data()
  df$grp3 <- sample(c("x", "y", "z"), nrow(df), replace = TRUE)
  jb <- joint_balance(df, group = "grp3", variables = c("A", "B", "C"),
                      n_perm = 50, seed = 1, verbose = FALSE)
  expect_equal(nrow(jb$pairwise_omnibus), 3L)  # choose(3, 2)
  expect_equal(length(jb$modality_graphs), 3L)
})

test_that("print method runs without error", {
  df <- make_small_data()
  jb <- joint_balance(df, group = "grp", n_perm = 50, seed = 1,
                      verbose = FALSE)
  expect_output(print(jb), "jointbalance diagnostic")
  expect_output(print(jb), "Marginal tests")
  expect_output(print(jb), "Pairwise omnibus tests")
})

test_that("summary method returns expected list", {
  df <- make_small_data()
  jb <- joint_balance(df, group = "grp", n_perm = 50, seed = 1,
                      verbose = FALSE)
  s <- summary(jb)
  expect_named(s, c("marginal_table", "pairwise_omnibus_table",
                    "n_edgewise_pairs", "group", "group_levels", "alpha"))
})

test_that("plot method runs with side_by_side and marginal_only", {
  df <- make_small_data()
  jb <- joint_balance(df, group = "grp", n_perm = 200, n_perm_edge = 200,
                      seed = 1, verbose = FALSE)
  # marginal_only always works
  expect_error(plot(jb, layout = "marginal_only"), NA)
  # side_by_side only works if we have an edgewise result
  if (length(jb$pairwise_edgewise) > 0L) {
    expect_error(plot(jb, layout = "side_by_side"), NA)
  }
})