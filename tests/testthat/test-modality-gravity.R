# tests/testthat/test-modality-gravity.R
# Tests for modality_gravity(), print/summary methods, plot_gravity(),
# plot_gravity_scatter(), and compare_gravity().

library(catgraph)

# ---- Shared fixtures --------------------------------------------------------

.make_mg <- function(pruned = TRUE) {
  data(survey_health, envir = environment())
  mg <- build_modality_graph(survey_health)
  if (pruned) mg <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
  mg
}

.make_mg_clustered <- function() {
  mg <- .make_mg()
  cluster_modalities(mg, method = "louvain")
}

.make_conditional <- function(sex_val) {
  data(survey_health, envir = environment())
  mg <- build_conditional_modality_graph(survey_health,
                                         given = list(sex = sex_val))
  prune_modality_edges(mg, min_weight = 0.05, max_p = 0.10)
}


# =============================================================================
# modality_gravity()
# =============================================================================

test_that("modality_gravity returns a data frame with required columns", {
  mg   <- .make_mg()
  grav <- modality_gravity(mg)
  
  expect_s3_class(grav, "data.frame")
  
  required_cols <- c("node", "variable", "modality", "prevalence",
                     "degree", "strength",
                     "mgi_plus", "mgi_minus", "delta_mgi",
                     "mgi_plus_norm", "os", "role")
  expect_true(all(required_cols %in% names(grav)),
              info = paste("Missing columns:",
                           paste(setdiff(required_cols, names(grav)),
                                 collapse = ", ")))
})

test_that("modality_gravity row count matches vertex count", {
  mg   <- .make_mg()
  grav <- modality_gravity(mg)
  expect_equal(nrow(grav), igraph::vcount(mg$graph))
})

test_that("modality_gravity prevalences sum to 1 within each variable", {
  mg   <- .make_mg()
  grav <- modality_gravity(mg)
  
  by_var <- split(grav$prevalence, grav$variable)
  sums   <- vapply(by_var, sum, numeric(1L))
  expect_true(all(abs(sums - 1) < 1e-3),   # rounding at 4dp allows ~0.001
              info = paste("Variable sums:", paste(round(sums, 4), collapse = " ")))
})

test_that("modality_gravity mgi_plus is non-negative", {
  grav <- modality_gravity(.make_mg())
  expect_true(all(grav$mgi_plus >= 0))
})

test_that("modality_gravity mgi_minus is non-negative", {
  grav <- modality_gravity(.make_mg())
  expect_true(all(grav$mgi_minus >= 0))
})

test_that("modality_gravity delta_mgi equals mgi_plus minus mgi_minus", {
  grav <- modality_gravity(.make_mg())
  # Both sides are rounded to 4dp independently so allow 1e-3 tolerance
  expect_equal(grav$delta_mgi,
               round(grav$mgi_plus - grav$mgi_minus, 4L),
               tolerance = 1e-3)
})

test_that("modality_gravity mgi_plus_norm is in [0, 1]", {
  grav <- modality_gravity(.make_mg())
  expect_true(all(grav$mgi_plus_norm >= 0 & grav$mgi_plus_norm <= 1))
})

test_that("modality_gravity mgi_plus_norm maximum is 1", {
  grav <- modality_gravity(.make_mg())
  expect_equal(max(grav$mgi_plus_norm), 1)
})

test_that("modality_gravity role is consistent with delta_mgi sign", {
  grav <- modality_gravity(.make_mg())
  expect_true(all(grav$role[grav$delta_mgi > 0]  == "attractor"))
  expect_true(all(grav$role[grav$delta_mgi < 0]  == "satellite"))
  expect_true(all(grav$role[grav$delta_mgi == 0] == "neutral"))
})

test_that("modality_gravity os is non-negative", {
  grav <- modality_gravity(.make_mg())
  expect_true(all(grav$os >= 0))
})

test_that("modality_gravity is ordered by delta_mgi descending", {
  grav <- modality_gravity(.make_mg())
  expect_true(all(diff(grav$delta_mgi) <= 0))
})

test_that("modality_gravity errors on non-catmodgraph input", {
  expect_error(modality_gravity(list(a = 1)), "catmodgraph")
})

test_that("modality_gravity errors when data slot is missing", {
  mg      <- .make_mg()
  mg$data <- NULL
  expect_error(modality_gravity(mg), "data")
})

test_that("modality_gravity works on unpruned graph", {
  data(survey_health, envir = environment())
  mg   <- build_modality_graph(survey_health)
  grav <- modality_gravity(mg)
  expect_s3_class(grav, "data.frame")
  expect_equal(nrow(grav), igraph::vcount(mg$graph))
})

test_that("modality_gravity isolated nodes (degree 0) have mgi and os of zero", {
  mg   <- .make_mg()
  grav <- modality_gravity(mg)
  iso  <- grav[grav$degree == 0L, ]
  if (nrow(iso) > 0L) {
    expect_true(all(iso$mgi_plus  == 0))
    expect_true(all(iso$mgi_minus == 0))
    expect_true(all(iso$os        == 0))
    expect_true(all(iso$role      == "neutral"))
  }
})

test_that("modality_gravity dominant node has highest mgi_plus", {
  # lung_disease=no is 85% prevalent and should dominate
  grav <- modality_gravity(.make_mg())
  top  <- grav$node[1L]
  expect_equal(grav$mgi_plus[1L], max(grav$mgi_plus))
})


# =============================================================================
# print.modality_gravity
# =============================================================================

test_that("print.modality_gravity outputs without error", {
  grav <- modality_gravity(.make_mg())
  expect_true(inherits(grav, "modality_gravity"))
  expect_output(print(grav), "Modality Gravity Index")
})

test_that("print.modality_gravity shows ATTRACTORS section when present", {
  grav <- modality_gravity(.make_mg())
  expect_output(print(grav), "ATTRACTORS")
})

test_that("print.modality_gravity shows SATELLITES section when present", {
  grav <- modality_gravity(.make_mg())
  if (any(grav$role == "satellite")) {
    expect_output(print(grav), "SATELLITES")
  }
})

test_that("print.modality_gravity returns invisibly", {
  grav <- modality_gravity(.make_mg())
  ret  <- withVisible(print(grav))
  expect_false(ret$visible)
})


# =============================================================================
# summary.modality_gravity
# =============================================================================

test_that("summary.modality_gravity returns a list with required elements", {
  grav <- modality_gravity(.make_mg())
  s    <- suppressMessages(suppressWarnings(
    capture.output(s <- summary(grav))
  ))
  # summary() prints and returns invisibly — call directly
  s <- summary(grav)
  expect_type(s, "list")
  expect_true(all(c("role_counts", "by_variable",
                    "top_attractor", "top_satellite",
                    "spearman_strength") %in% names(s)))
})

test_that("summary role_counts sums to nrow of gravity df", {
  grav <- modality_gravity(.make_mg())
  s    <- summary(grav)
  expect_equal(sum(s$role_counts), nrow(grav))
})

test_that("summary by_variable has one row per unique variable", {
  grav <- modality_gravity(.make_mg())
  s    <- summary(grav)
  expect_equal(nrow(s$by_variable), length(unique(grav$variable)))
})

test_that("summary spearman_strength is in [-1, 1]", {
  grav <- modality_gravity(.make_mg())
  s    <- summary(grav)
  expect_true(s$spearman_strength >= -1 & s$spearman_strength <= 1)
})

test_that("summary top_attractor matches row with max delta_mgi", {
  grav <- modality_gravity(.make_mg())
  s    <- summary(grav)
  expect_equal(s$top_attractor$delta_mgi, max(grav$delta_mgi))
})

test_that("summary top_satellite matches row with min delta_mgi", {
  grav <- modality_gravity(.make_mg())
  s    <- summary(grav)
  expect_equal(s$top_satellite$delta_mgi, min(grav$delta_mgi))
})


# =============================================================================
# plot_gravity()
# =============================================================================

test_that("plot_gravity runs without error on clustered catmodgraph", {
  mg <- .make_mg_clustered()
  expect_silent(plot_gravity(mg))
})

test_that("plot_gravity runs without error when cluster attr absent", {
  mg <- .make_mg()
  expect_silent(plot_gravity(mg))
})

test_that("plot_gravity accepts pre-computed gravity", {
  mg   <- .make_mg_clustered()
  grav <- modality_gravity(mg)
  expect_silent(plot_gravity(mg, gravity = grav))
})

test_that("plot_gravity returns gravity df invisibly", {
  mg  <- .make_mg_clustered()
  ret <- plot_gravity(mg)
  expect_s3_class(ret, "data.frame")
  expect_true("delta_mgi" %in% names(ret))
})

test_that("plot_gravity bars=TRUE runs without error", {
  mg <- .make_mg_clustered()
  # bars=TRUE opens a new plot window which may fail in headless environments
  result <- tryCatch(
    plot_gravity(mg, bars = TRUE, bars_n = 8L),
    error = function(e) {
      # Accept "figure margins too large" as a graphics device issue not a code bug
      if (grepl("figure margins|plot\\.new", conditionMessage(e))) {
        "skipped_headless"
      } else {
        stop(e)
      }
    }
  )
  # Either it ran successfully (returns df) or was skipped due to graphics
  expect_true(is.data.frame(result) || identical(result, "skipped_headless"))
})

test_that("plot_gravity errors on non-catmodgraph", {
  expect_error(plot_gravity(list()), "catmodgraph")
})


# =============================================================================
# plot_gravity_scatter()
# =============================================================================

test_that("plot_gravity_scatter runs without error", {
  mg   <- .make_mg()
  grav <- modality_gravity(mg)
  expect_silent(plot_gravity_scatter(grav, mg))
})

test_that("plot_gravity_scatter returns a data frame invisibly", {
  mg   <- .make_mg()
  grav <- modality_gravity(mg)
  ret  <- plot_gravity_scatter(grav, mg)
  expect_s3_class(ret, "data.frame")
  expect_true(all(c("node", "eigenvec", "delta_mgi",
                    "role", "is_contradiction") %in% names(ret)))
})

test_that("plot_gravity_scatter eigenvec values in [0, 1]", {
  mg   <- .make_mg()
  grav <- modality_gravity(mg)
  ret  <- plot_gravity_scatter(grav, mg)
  expect_true(all(ret$eigenvec >= 0 & ret$eigenvec <= 1))
})

test_that("plot_gravity_scatter errors on non-catmodgraph", {
  grav <- modality_gravity(.make_mg())
  expect_error(plot_gravity_scatter(grav, list()), "catmodgraph")
})


# =============================================================================
# compare_gravity()
# =============================================================================

test_that("compare_gravity returns a data frame", {
  mg_f <- .make_conditional("female")
  mg_m <- .make_conditional("male")
  cmp  <- compare_gravity(list(female = mg_f, male = mg_m), plot = FALSE)
  expect_s3_class(cmp, "data.frame")
})

test_that("compare_gravity output has delta_mgi_diff column", {
  mg_f <- .make_conditional("female")
  mg_m <- .make_conditional("male")
  cmp  <- compare_gravity(list(female = mg_f, male = mg_m), plot = FALSE)
  expect_true("delta_mgi_diff" %in% names(cmp))
})

test_that("compare_gravity delta_mgi_diff equals g1 minus g2", {
  mg_f <- .make_conditional("female")
  mg_m <- .make_conditional("male")
  cmp  <- compare_gravity(list(female = mg_f, male = mg_m), plot = FALSE)
  
  expect_equal(
    cmp$delta_mgi_diff,
    cmp$delta_mgi_female - cmp$delta_mgi_male,
    tolerance = 1e-6
  )
})

test_that("compare_gravity is ordered by abs(delta_mgi_diff) descending", {
  mg_f <- .make_conditional("female")
  mg_m <- .make_conditional("male")
  cmp  <- compare_gravity(list(female = mg_f, male = mg_m), plot = FALSE)
  expect_true(all(diff(abs(cmp$delta_mgi_diff)) <= 0))
})

test_that("compare_gravity errors on non-list input", {
  expect_error(compare_gravity("not_a_list"), "list")
})

test_that("compare_gravity errors on list of length != 2", {
  mg <- .make_mg()
  expect_error(compare_gravity(list(a = mg)), "exactly two")
})

test_that("compare_gravity errors on unnamed list", {
  mg_f <- .make_conditional("female")
  mg_m <- .make_conditional("male")
  expect_error(compare_gravity(list(mg_f, mg_m)), "named list")
})

test_that("compare_gravity errors if elements are not catmodgraph", {
  mg_f <- .make_conditional("female")
  expect_error(compare_gravity(list(a = mg_f, b = list())), "catmodgraph")
})

test_that("compare_gravity plot=TRUE runs without error", {
  mg_f <- .make_conditional("female")
  mg_m <- .make_conditional("male")
  expect_silent(compare_gravity(list(female = mg_f, male = mg_m),
                                plot = TRUE, top_n = 10L))
})


# =============================================================================
# node_centrality.catmodgraph()
# =============================================================================

test_that("node_centrality dispatches to catmodgraph method", {
  mg  <- .make_mg()
  nc  <- node_centrality(mg)
  expect_s3_class(nc, "data.frame")
  expect_true("delta_mgi" %in% names(nc))
  expect_true("os"        %in% names(nc))
  expect_true("role"      %in% names(nc))
})

test_that("node_centrality.catmodgraph has correct column set", {
  mg  <- .make_mg()
  nc  <- node_centrality(mg)
  expected <- c("node", "variable", "modality", "prevalence", "degree",
                "strength", "w_betweenness", "w_closeness",
                "w_eigenvector", "w_pagerank",
                "mgi_plus", "mgi_minus", "delta_mgi",
                "mgi_plus_norm", "os", "role")
  expect_true(all(expected %in% names(nc)),
              info = paste("Missing:", paste(setdiff(expected, names(nc)),
                                             collapse = ", ")))
})

test_that("node_centrality.catmodgraph row count matches vertex count", {
  mg  <- .make_mg()
  nc  <- node_centrality(mg)
  expect_equal(nrow(nc), igraph::vcount(mg$graph))
})

test_that("node_centrality.catmodgraph normalized values in [0, 1]", {
  mg  <- .make_mg()
  nc  <- node_centrality(mg, normalize = TRUE)
  trad_cols <- c("strength", "w_betweenness", "w_closeness",
                 "w_eigenvector", "w_pagerank")
  for (col in trad_cols) {
    v <- nc[[col]]
    v <- v[!is.na(v)]   # closeness may be NA for disconnected nodes
    expect_true(all(v >= 0 & v <= 1),
                info = paste(col, "out of [0,1]"))
  }
})

test_that("node_centrality.catmodgraph unnormalized strength is non-negative", {
  mg  <- .make_mg()
  nc  <- node_centrality(mg, normalize = FALSE)
  expect_true(all(nc$strength >= 0))
})

test_that("node_centrality.catmodgraph ordered by delta_mgi descending", {
  mg  <- .make_mg()
  nc  <- node_centrality(mg)
  expect_true(all(diff(nc$delta_mgi) <= 0))
})