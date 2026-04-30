# Bootstrap confidence intervals for phi or Cramer's V

Estimates a confidence interval for a pairwise effect size (phi or
Cramer's V) using a non-parametric percentile bootstrap. Rows are
resampled with replacement from the pairwise-complete observations, and
the chosen effect size estimator is recomputed on each resample. The
resulting empirical distribution is used to derive the interval.

## Usage

``` r
bootstrap_ci(
  x,
  y,
  R = 1000L,
  conf = 0.95,
  type = c("percentile", "bca"),
  corrected = FALSE,
  correct = FALSE,
  seed = NULL
)
```

## Arguments

- x:

  A factor, character, or logical vector.

- y:

  A factor, character, or logical vector of the same length as `x`.

- R:

  Integer. Number of bootstrap resamples. Default `1000L`. Values below
  500 are accepted with a warning.

- conf:

  Numeric in (0, 1). Confidence level. Default `0.95`.

- type:

  Character. Bootstrap interval type: `"percentile"` (default) or
  `"bca"` (bias-corrected and accelerated, Efron 1987). The BCa interval
  generally has better coverage but requires the jackknife influence
  values and is slower.

- corrected:

  Logical. Whether to use the bias-corrected estimator (Bergsma, 2013).
  Default `FALSE`.

- correct:

  Logical. Yates' continuity correction. Default `FALSE`.

- seed:

  Integer or `NULL`. Optional random seed for reproducibility. Default
  `NULL` (no seed set).

## Value

A named list with:

- `estimate`:

  The point estimate of phi or Cramer's V on the original data.

- `ci_lower`:

  Lower confidence bound.

- `ci_upper`:

  Upper confidence bound.

- `conf`:

  The requested confidence level.

- `type`:

  The interval type used.

- `R`:

  Number of resamples actually used (may be lower than requested if
  degenerate resamples were removed).

- `metric`:

  Character: `"phi"` or `"cramers_v"`.

- `corrected`:

  Logical.

- `n`:

  Number of pairwise-complete observations.

- `boot_dist`:

  Numeric vector of length `R` containing the full bootstrap
  distribution (useful for plotting).

## Details

**Percentile interval** (Efron & Tibshirani, 1993, Ch. 13):

The \\\alpha/2\\ and \\1-\alpha/2\\ quantiles of the bootstrap
distribution \\\hat{\theta}^\*\_1, \ldots, \hat{\theta}^\*\_R\\ are used
directly as the lower and upper bounds.

**BCa interval** (Efron, 1987):

Adjusts the quantiles used for the interval by a bias-correction term
\\\hat{z}\_0\\ (estimated from the proportion of bootstrap resamples
below the point estimate) and an acceleration term \\\hat{a}\\
(estimated from the jackknife influence values). This gives second-order
accurate intervals without requiring a transformation. The formulas
follow DiCiccio & Efron (1996).

**Degenerate resamples**: Bootstrap resamples that yield a single-level
variable (all observations in one category) produce `NA` effect sizes
and are discarded. A warning is issued if more than 5% of resamples are
discarded.

**Effect size on a boundary**: Both phi and Cramer's V are bounded at 0.
Bootstrap distributions for weak associations are therefore right-skewed
and pile up at 0. The BCa correction partially addresses this; for
near-zero associations the percentile lower bound may be 0, which is
correct.

## References

DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
*Statistical Science*, 11(3), 189–212.
[doi:10.1214/ss/1032280214](https://doi.org/10.1214/ss/1032280214)

Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the
American Statistical Association*, 82(397), 171–185.
[doi:10.1080/01621459.1987.10478410](https://doi.org/10.1080/01621459.1987.10478410)

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall.
[doi:10.1007/978-1-4899-4541-9](https://doi.org/10.1007/978-1-4899-4541-9)

## See also

[`effect_size`](https://atinakosta.github.io/catgraph/reference/effect_size.md),
[`catgraph_ci`](https://atinakosta.github.io/catgraph/reference/catgraph_ci.md)

## Examples

``` r
set.seed(42)
x <- sample(c("A", "B"), 120, replace = TRUE,
            prob = c(0.4, 0.6))
y <- ifelse(x == "A",
            sample(c("yes", "no"), 120, replace = TRUE, prob = c(0.7, 0.3)),
            sample(c("yes", "no"), 120, replace = TRUE, prob = c(0.3, 0.7)))

# Small R for a fast example; use R >= 1000 in real work.
ci <- bootstrap_ci(x, y, R = 200, seed = 1)
#> Warning: R < 500 may yield unstable confidence intervals.
ci$estimate
#> [1] 0.4323126
ci$ci_lower
#> [1] 0.2452065
ci$ci_upper
#> [1] 0.5829789

# \donttest{
# BCa interval (slower: adds a leave-one-out jackknife)
ci_bca <- bootstrap_ci(x, y, R = 500, type = "bca", seed = 1)
# }
```
