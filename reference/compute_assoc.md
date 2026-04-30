# Compute chi-square association between two categorical variables

Performs a Pearson chi-square test of independence on the
pairwise-complete observations of two categorical variables. Returns the
test statistic, p-value, degrees of freedom, and the contingency table
used.

## Usage

``` r
compute_assoc(
  x,
  y,
  correct = FALSE,
  simulate_p = FALSE,
  B = 2000L,
  x_name = NULL,
  y_name = NULL
)
```

## Arguments

- x:

  A factor, character, or logical vector.

- y:

  A factor, character, or logical vector of the same length as `x`.

- correct:

  Logical. Whether to apply Yates' continuity correction for 2x2 tables.
  Default `FALSE` to keep results comparable with effect size formulas
  that do not incorporate the correction (Agresti, 2002, p. 77).

- simulate_p:

  Logical. If `TRUE`, the p-value is estimated by Monte Carlo
  simulation. Recommended when expected cell frequencies are small.
  Default `FALSE`.

- B:

  Integer. Monte Carlo resamples when `simulate_p = TRUE`. Default
  `2000L`.

- x_name, y_name:

  Optional character. Variable names used in warning messages to
  identify which pair triggered a sparsity warning. If `NULL` (default),
  the names are deparsed from the call. Primarily used internally by
  [`build_graph`](https://atinakosta.github.io/catgraph/reference/build_graph.md);
  users do not normally supply these.

## Value

A named list with:

- statistic:

  The chi-square test statistic.

- p_value:

  The p-value of the test.

- df:

  Degrees of freedom (NA under simulate_p = TRUE).

- n:

  Number of pairwise-complete observations.

- table:

  The contingency table.

- type:

  Table type from
  [`detect_type`](https://atinakosta.github.io/catgraph/reference/detect_type.md).

## Details

Rows with `NA` in either `x` or `y` are removed before tabulation
(pairwise deletion). Two warning conditions are checked:

1.  **Any expected count \< 5** (Cochran, 1954): a minor sparsity
    warning. Consider `simulate_p = TRUE`.

2.  **Severe sparsity** - either more than 20% of cells have expected
    count \< 5, or in a table of minimum dimension \\\geq 3\\ the
    observations-to-cells ratio is below 5. In this regime the
    chi-square approximation is unreliable *and* Cramer's V becomes
    numerically unstable, especially in its classical (uncorrected)
    form. Consider collapsing categories, enabling bias correction
    (`corrected = TRUE` in
    [`effect_size`](https://atinakosta.github.io/catgraph/reference/effect_size.md)),
    or using `simulate_p = TRUE`.

**Comparability across table dimensions.** Both phi and Cramer's V are
normalised to \\\[0, 1\]\\, but this does *not* make them directly
comparable across tables of different dimensions. Under independence,
the sampling distribution of V depends on table dimension and sample
size: for the same true association strength, V from a \\5 \times 5\\
table is not exchangeable with V from a \\2 \times 2\\ table. A V of
0.25 observed on a \\5 \times 5\\ table and a V of 0.25 observed on a
\\2 \times 2\\ table represent comparable values on the normalised scale
but may correspond to different sampling-distribution quantiles under
independence. In a `catgraph` object, variables with many more
categories than their neighbours may therefore score higher on Cramer's
V at any given level of true dependence, which can inflate their
apparent centrality. Bias correction via `corrected = TRUE` (Bergsma,
2013) partially mitigates this. See the package vignette,
"Methodological caveats", item 2 (mixed metrics).

## References

Agresti, A. (2002). *Categorical Data Analysis* (2nd ed.).
[doi:10.1002/0471249688](https://doi.org/10.1002/0471249688)

Cochran, W. G. (1954). Some methods for strengthening the common
chi-square tests. *Biometrics*, 10(4), 417–451.
[doi:10.2307/3001616](https://doi.org/10.2307/3001616)

## See also

[`effect_size`](https://atinakosta.github.io/catgraph/reference/effect_size.md),
[`detect_type`](https://atinakosta.github.io/catgraph/reference/detect_type.md)

## Examples

``` r
set.seed(42)
x <- sample(c("A", "B"), 100, replace = TRUE)
y <- sample(c("yes", "no"), 100, replace = TRUE)
result <- compute_assoc(x, y)
result$statistic
#> [1] 0.5312869
result$p_value
#> [1] 0.4660663
```
