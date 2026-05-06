# Bayesian Cramér's V for a pair of categorical variables

Computes a Bayesian estimate of Cramér's V by applying a symmetric
Dirichlet prior to the contingency table cell counts before computing
the association measure. This shrinks edge weights toward zero for
sparse tables, producing more stable estimates than classical or
bias-corrected Cramér's V when expected cell frequencies are small.

## Usage

``` r
bayesian_cramers_v(x, y, alpha = 0.5)
```

## Arguments

- x:

  A factor, character, or logical vector.

- y:

  A factor, character, or logical vector of the same length as `x`.

- alpha:

  Numeric. Dirichlet prior concentration parameter added to each cell
  count before computing the association. Must be \> 0. Default `0.5`
  (Jeffreys prior). Use `alpha = 1` for the Laplace (uniform) prior.

## Value

A named list with:

- `effect_size`:

  Bayesian Cramér's V, numeric in \\\[0, 1\]\\.

- `metric`:

  Character: `"bayesian_cramers_v"`.

- `alpha`:

  The prior concentration used.

- `type`:

  Contingency-table type: `"2x2"`, `"RxC"`, or `"degenerate"`.

- `statistic`:

  The chi-square statistic computed on the **smoothed** table, for
  reference.

- `p_value`:

  Chi-square p-value from the **original** (unsmoothed) table. Smoothing
  is applied only to the effect-size estimate, not to the test.

- `df`:

  Degrees of freedom.

- `n`:

  Number of pairwise-complete observations (unsmoothed).

## Details

**Dirichlet smoothing.** Under a symmetric Dirichlet(\\\alpha\\) prior
on the \\r \times c\\ cell probability vector, the posterior mean of
each cell probability is:

\$\$\hat{p}\_{ij} = \frac{n\_{ij} + \alpha}{n + \alpha \cdot r \cdot
c}\$\$

Cramér's V is then computed from the chi-square statistic derived from
these smoothed proportions rather than from the raw counts. This is
equivalent to computing Cramér's V on a pseudo-count table
\\\tilde{n}\_{ij} = n\_{ij} + \alpha\\ with effective sample size
\\\tilde{n} = n + \alpha \cdot r \cdot c\\.

**Jeffreys prior (\\\alpha = 0.5\\).** This is the standard
non-informative choice for categorical data. It corresponds to adding
half a pseudocount to each cell, which stabilises the chi-square
statistic for sparse tables without materially distorting the estimate
when tables are well-populated.

**Relationship to classical Cramér's V.** As \\n \to \infty\\, the
smoothed estimate converges to the classical estimator because the
pseudocounts \\\alpha\\ become negligible relative to \\n\\. For large
samples (such as the Titanic dataset with n = 2201) the difference is
therefore very small. The practical advantage of the Bayesian estimator
appears on small samples or sparse contingency tables.

**p-value.** The p-value is taken from the unsmoothed chi-square test
(via
[`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md)).
This is intentional: smoothing inflates the effective sample size and
would otherwise produce anti-conservative p-values. Users who want a
fully Bayesian decision criterion should use posterior credible
intervals (not yet implemented) rather than the p-value.

## References

Good, I. J. (1965). *The Estimation of Probabilities: An Essay on Modern
Bayesian Methods*. MIT Press.

Agresti, A. (2002). *Categorical Data Analysis* (2nd ed.). Wiley.
[doi:10.1002/0471249688](https://doi.org/10.1002/0471249688)

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

## See also

[`effect_size`](https://atinakosta.github.io/catgraph/reference/effect_size.md),
[`nmi_assoc`](https://atinakosta.github.io/catgraph/reference/nmi_assoc.md),
[`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md)

## Examples

``` r
set.seed(1)
x <- sample(c("A", "B", "C"), 120, replace = TRUE)
y <- sample(c("yes", "no"),    120, replace = TRUE)

# Classical
effect_size(x, y)$effect_size
#> [1] 0.170973

# Bayesian (Jeffreys prior)
bayesian_cramers_v(x, y)$effect_size
#> [1] 0.1669501

# Bayesian (Laplace prior)
bayesian_cramers_v(x, y, alpha = 1)$effect_size
#> [1] 0.1631126

# On a sparse table: Bayesian estimate is more stable
x_sparse <- sample(c("A","B","C","D"), 20, replace = TRUE)
y_sparse <- sample(c("P","Q","R","S"), 20, replace = TRUE)
effect_size(x_sparse, y_sparse)$effect_size
#> Warning: At least one expected cell frequency is < 5 for pair (x, y). Consider setting simulate_p = TRUE.
#> Warning: Sparse contingency table for pair (x, y): 3x4 = 12 cells, 20 obs, 100% cells with E < 5. Cramer's V and chi-square p-values may be unstable; consider collapsing categories, enabling bias correction, or simulate_p = TRUE.
#> [1] 0.5672383
bayesian_cramers_v(x_sparse, y_sparse)$effect_size
#> Warning: At least one expected cell frequency is < 5 for pair (x, y). Consider setting simulate_p = TRUE.
#> Warning: Sparse contingency table for pair (x, y): 3x4 = 12 cells, 20 obs, 100% cells with E < 5. Cramer's V and chi-square p-values may be unstable; consider collapsing categories, enabling bias correction, or simulate_p = TRUE.
#> [1] 0.4377301
```
