# ============================================================================
# Rerun of Experiments 2 and 3 with FIXED DGP
# ============================================================================

library(catgraph)
library(igraph)
library(ggplot2)

# Corrected DGP
gen_null <- function(n) {
  z <- sample(1:3, n, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  data.frame(
    a = factor(z),
    b = factor(ifelse(runif(n) < 0.7, z, sample(1:3, n, replace = TRUE))),
    c = factor(sample(1:3, n, replace = TRUE))
  )
}
gen_alt <- function(n, strength) {
  z <- sample(1:3, n, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  a <- factor(z)
  b <- factor(ifelse(runif(n) < 0.7, z, sample(1:3, n, replace = TRUE)))
  # CORRECTED: at strength=0 this matches gen_null exactly
  c_vals <- ifelse(runif(n) < strength, z, sample(1:3, n, replace = TRUE))
  data.frame(a = a, b = b, c = factor(c_vals))
}

run_one_test <- function(df_x, df_y, n_perm = 100, statistic = "frobenius") {
  mg_x <- build_modality_graph(df_x)
  mg_y <- build_modality_graph(df_y)
  res <- tryCatch(
    test_modality_graph_equality(
      mg_x, mg_y, n_perm = n_perm,
      statistic = statistic, verbose = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(res)) NA_real_ else res$p_value
}

set.seed(2026)
n_perm <- 100


# ---- EXPERIMENT 2: Power curve ----
cat("\n#### EXPERIMENT 2: POWER CURVE (corrected) ####\n")

strengths    <- c(0, 0.2, 0.4, 0.6, 0.8)
n_reps_power <- 50
n_power      <- 300

power_results <- data.frame(
  strength = strengths, reject_rate = NA_real_, se_reject = NA_real_
)

t_start <- Sys.time()
for (i in seq_along(strengths)) {
  s_i <- strengths[i]
  cat(sprintf("  strength=%.1f ... ", s_i))
  t0 <- Sys.time()
  pvals <- replicate(n_reps_power, {
    run_one_test(gen_null(n_power), gen_alt(n_power, s_i),
                 n_perm = n_perm)
  })
  pvals <- pvals[!is.na(pvals)]
  rr <- mean(pvals < 0.05)
  se <- sqrt(rr * (1 - rr) / length(pvals))
  power_results$reject_rate[i] <- rr
  power_results$se_reject[i]   <- se
  cat(sprintf("reject = %.3f (%.0fs)\n",
              rr, as.numeric(Sys.time() - t0, units = "secs")))
}
cat(sprintf("Experiment 2 total: %.1f min\n\n",
            as.numeric(Sys.time() - t_start, units = "mins")))
print(power_results, digits = 3)


# ---- EXPERIMENT 3: Sample-size scaling ----
cat("\n#### EXPERIMENT 3: SAMPLE-SIZE SCALING (corrected) ####\n")

n_scales       <- c(150, 400, 800)
n_reps_ss      <- 40
strength_fixed <- 0.4

scaling_results <- data.frame(
  n = n_scales, reject_rate = NA_real_, se_reject = NA_real_
)

t_start <- Sys.time()
for (i in seq_along(n_scales)) {
  n_i <- n_scales[i]
  cat(sprintf("  n=%d ... ", n_i))
  t0 <- Sys.time()
  pvals <- replicate(n_reps_ss, {
    run_one_test(gen_null(n_i), gen_alt(n_i, strength_fixed),
                 n_perm = n_perm)
  })
  pvals <- pvals[!is.na(pvals)]
  rr <- mean(pvals < 0.05)
  se <- sqrt(rr * (1 - rr) / length(pvals))
  scaling_results$reject_rate[i] <- rr
  scaling_results$se_reject[i]   <- se
  cat(sprintf("reject = %.3f (%.0fs)\n",
              rr, as.numeric(Sys.time() - t0, units = "secs")))
}
cat(sprintf("Experiment 3 total: %.1f min\n\n",
            as.numeric(Sys.time() - t_start, units = "mins")))
print(scaling_results, digits = 3)


# ---- Plots ----
p2 <- ggplot(power_results, aes(x = strength, y = reject_rate)) +
  geom_line(color = "#534AB7", linewidth = 1.2) +
  geom_point(color = "#534AB7", size = 3) +
  geom_errorbar(aes(ymin = reject_rate - 1.96 * se_reject,
                    ymax = reject_rate + 1.96 * se_reject),
                width = 0.03, color = "#534AB7") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "#D85A30") +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "grey30") +
  labs(title = "Power curve (corrected DGP)",
       subtitle = sprintf("n = %d per group, statistic = frobenius", n_power),
       x = "Effect strength",
       y = "Rejection rate") +
  ylim(0, 1) +
  theme_minimal(base_size = 11)
print(p2)

p3 <- ggplot(scaling_results, aes(x = n, y = reject_rate)) +
  geom_line(color = "#534AB7", linewidth = 1.2) +
  geom_point(color = "#534AB7", size = 3) +
  geom_errorbar(aes(ymin = reject_rate - 1.96 * se_reject,
                    ymax = reject_rate + 1.96 * se_reject),
                width = 20, color = "#534AB7") +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "grey30") +
  scale_x_log10() +
  labs(title = "Sample-size scaling (corrected DGP)",
       subtitle = sprintf("Fixed strength = %.1f, statistic = frobenius",
                          strength_fixed),
       x = "Sample size (log scale)",
       y = "Rejection rate") +
  ylim(0, 1) +
  theme_minimal(base_size = 11)
print(p3)


# ---- Assessment ----
cat("\n#### ASSESSMENT (corrected) ####\n")

# At strength=0, reject rate should be near 0.05 (new sanity check)
null_in_alt <- power_results$reject_rate[power_results$strength == 0]
cat(sprintf("Reject rate at strength=0 (should be ~0.05): %.3f %s\n",
            null_in_alt,
            ifelse(abs(null_in_alt - 0.05) < 0.08, "OK", "PROBLEM")))

pwr_monotone <- all(diff(power_results$reject_rate) >= -0.05)
cat("Power rises monotonically with strength? ",
    ifelse(pwr_monotone, "YES", "NO"), "\n")

pwr_max <- tail(power_results$reject_rate, 1) > 0.70
cat("Power > 70% at max strength, n=300?      ",
    ifelse(pwr_max, "YES", "NO"), "\n")

scale_rising <- tail(scaling_results$reject_rate, 1) >
  head(scaling_results$reject_rate, 1)
cat("Power rises with n at fixed effect?      ",
    ifelse(scale_rising, "YES", "NO"), "\n")