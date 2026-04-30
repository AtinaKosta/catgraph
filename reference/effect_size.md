# Compute phi or Cramer's V effect size for a pair of categorical variables

Dispatches to the phi coefficient for 2x2 tables and to Cramer's V for
larger tables. Optionally applies the bias correction of Bergsma (2013).
All calculations are based on the chi-square statistic and sample size
returned by
[`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md).

## Usage

``` r
effect_size(
  x,
  y,
  corrected = FALSE,
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

- corrected:

  Logical. If `TRUE`, returns the bias-corrected version of Cramer's V
  (Bergsma, 2013) for n x m tables, or the bias-corrected phi for 2x2
  tables. If `FALSE` (default), returns the classical estimators.

- correct:

  Logical. Yates' continuity correction for chi-square. Passed to
  [`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md).
  Default `FALSE`.

- simulate_p:

  Logical. Monte Carlo p-value simulation. Passed to
  [`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md).
  Default `FALSE`.

- B:

  Integer. Monte Carlo resamples. Default `2000L`.

- x_name, y_name:

  Optional character. Variable names used in warning messages to
  identify which pair triggered a sparsity warning. If `NULL` (default),
  the names are deparsed from the call. Primarily used internally by
  [`build_graph`](https://atinakosta.github.io/catgraph/reference/build_graph.md);
  users do not normally supply these.

## Value

A named list with:

- effect_size:

  The phi or Cramer's V value (numeric, in \[0, 1\] for classical; may
  be 0 when bias-corrected formula yields a negative value, which is
  clamped to 0).

- metric:

  Character: `"phi"` or `"cramers_v"`.

- corrected:

  Logical indicating whether bias correction was applied.

- type:

  Observed contingency-table type: `"2x2"`, `"RxC"`, or `"degenerate"`.

- statistic:

  The chi-square statistic.

- p_value:

  The p-value from the chi-square test.

- df:

  Degrees of freedom (NA if Monte Carlo simulation was used).

- n:

  Number of pairwise-complete observations.

## Details

**Phi coefficient (2x2 tables):**

\$\$\phi = \sqrt{\chi^2 / n}\$\$

This is equivalent to the Pearson correlation coefficient computed on
two binary variables coded 0/1 (Agresti, 2002, p. 60).

**Bias-corrected phi (Bergsma, 2013):**

\$\$\tilde{\phi} = \max\\\left(0,\\ \phi^2 -
\frac{1}{n-1}\right)^{1/2}\$\$

**Classical Cramer's V (n x m tables):**

\$\$V = \sqrt{\frac{\chi^2 / n}{\min(r-1,\\ c-1)}}\$\$

where \\r\\ and \\c\\ are the number of rows and columns in the
contingency table (Cramer, 1946).

**Bias-corrected Cramer's V (Bergsma, 2013):**

Let \\\phi^2 = \chi^2/n\\. The bias-corrected estimator is computed as

\$\$\phi^2\_{corr} = \max\\\left(0,\\ \phi^2 -
\frac{(r-1)(c-1)}{n-1}\right)\$\$

with corrected effective dimensions

\$\$r\_{corr} = r - \frac{(r-1)^2}{n-1}, \qquad c\_{corr} = c -
\frac{(c-1)^2}{n-1}\$\$

and

\$\$\tilde{V} = \sqrt{\frac{\phi^2\_{corr}}{\min(r\_{corr}-1,\\
c\_{corr}-1)}}\$\$

whenever the denominator is positive. If the corrected denominator is
non-positive, the function returns 0 with a warning.

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

Agresti, A. (2002). *Categorical Data Analysis* (2nd ed.). John Wiley &
Sons. [doi:10.1002/0471249688](https://doi.org/10.1002/0471249688)

Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
*Journal of the Korean Statistical Society*, 42(3), 323–328.
[doi:10.1016/j.jkss.2012.10.002](https://doi.org/10.1016/j.jkss.2012.10.002)

Cramer, H. (1946). *Mathematical Methods of Statistics*. Princeton
University Press.

## See also

[`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md),
[`detect_type`](https://atinakosta.github.io/catgraph/reference/detect_type.md),
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md)

## Examples

``` r
set.seed(1)
x <- sample(c("A", "B"), 80, replace = TRUE)
y <- sample(c("yes", "no"), 80, replace = TRUE)
effect_size(x, y)
#> $effect_size
#> [1] 0.1001252
#> 
#> $metric
#> [1] "phi"
#> 
#> $corrected
#> [1] FALSE
#> 
#> $type
#> [1] "2x2"
#> 
#> $statistic
#> [1] 0.802005
#> 
#> $p_value
#> [1] 0.3704946
#> 
#> $df
#> [1] 1
#> 
#> $n
#> [1] 80
#> 
effect_size(x, y, corrected = TRUE)
#> $effect_size
#> [1] 0
#> 
#> $metric
#> [1] "phi"
#> 
#> $corrected
#> [1] TRUE
#> 
#> $type
#> [1] "2x2"
#> 
#> $statistic
#> [1] 0.802005
#> 
#> $p_value
#> [1] 0.3704946
#> 
#> $df
#> [1] 1
#> 
#> $n
#> [1] 80
#> 

# Multinomial example
z <- sample(c("low", "mid", "high"), 80, replace = TRUE)
effect_size(x, z)
#> $effect_size
#> [1] 0.09053383
#> 
#> $metric
#> [1] "cramers_v"
#> 
#> $corrected
#> [1] FALSE
#> 
#> $type
#> [1] "RxC"
#> 
#> $statistic
#> [1] 0.6557099
#> 
#> $p_value
#> [1] 0.7204675
#> 
#> $df
#> [1] 2
#> 
#> $n
#> [1] 80
#> 
effect_size(x, z, corrected = TRUE)
#> $effect_size
#> [1] 0
#> 
#> $metric
#> [1] "cramers_v"
#> 
#> $corrected
#> [1] TRUE
#> 
#> $type
#> [1] "RxC"
#> 
#> $statistic
#> [1] 0.6557099
#> 
#> $p_value
#> [1] 0.7204675
#> 
#> $df
#> [1] 2
#> 
#> $n
#> [1] 80
#> 
```
