# Dense pairwise similarity matrix of categorical variables

Computes the full \\p \times p\\ matrix of pairwise effect sizes for
categorical variables, including pairs with zero association.

## Usage

``` r
assoc_similarity(
  data,
  corrected = FALSE,
  correct = FALSE,
  simulate_p = FALSE,
  B = 2000L,
  what = c("effect_size", "p_value", "n", "all")
)
```

## Arguments

- data:

  A data frame of categorical variables (same requirements as
  [`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md)).

- corrected:

  Logical. If `TRUE`, use the bias-corrected estimator (Bergsma, 2013).
  Default `FALSE`.

- correct:

  Logical. Yates' continuity correction for chi-square. Default `FALSE`.

- simulate_p:

  Logical. Monte Carlo simulation for p-values (affects only the p-value
  matrix, not the effect-size matrix). Default `FALSE`.

- B:

  Integer. Monte Carlo resamples. Default `2000L`.

- what:

  Character. What to return: `"effect_size"` (default), `"p_value"`,
  `"n"` (pairwise-complete observation count), or `"all"` (a list of
  matrices).

## Value

A symmetric numeric matrix (or a list of matrices when `what = "all"`).
Diagonal is `NA`. Row and column names are the variable names.

## Details

This function is the **correct input for heatmap-style visualisation**
and any analysis that requires a dense similarity matrix. The igraph
object returned by
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
is the correct input for *topology* (centrality, clustering, density,
community detection): it represents zero-association pairs as absent
edges and therefore would give misleading heatmaps.

This function duplicates the computation done by
[`build_graph`](https://atinakosta.github.io/catgraph/reference/build_graph.md)
but does not collapse the result into a graph, so all pairs are
represented. In v0.3.0 and earlier, the same output was extracted from
the graph via
[`assoc_matrix()`](https://atinakosta.github.io/catgraph/reference/assoc_matrix.md),
but because the graph forced zero-weight pairs to `.Machine$double.eps`,
the resulting matrix silently conflated "zero association" with
"near-zero association". From 0.4.0 onwards, use `assoc_similarity()`
when you want the *full* dense matrix and
[`assoc_matrix()`](https://atinakosta.github.io/catgraph/reference/assoc_matrix.md)
(the graph extractor) when you want the matrix of *actual* edges.

## See also

[`build_graph`](https://atinakosta.github.io/catgraph/reference/build_graph.md),
[`assoc_matrix`](https://atinakosta.github.io/catgraph/reference/assoc_matrix.md),
[`plot_heatmap`](https://atinakosta.github.io/catgraph/reference/plot_heatmap.md)

## Examples

``` r
df <- expand_table(Titanic)
S <- assoc_similarity(df)
round(S, 3)
#>          Class   Sex   Age Survived
#> Class       NA 0.399 0.232    0.294
#> Sex      0.399    NA 0.111    0.456
#> Age      0.232 0.111    NA    0.098
#> Survived 0.294 0.456 0.098       NA

# All three components at once
out <- assoc_similarity(df, what = "all")
str(out, max.level = 1)
#> List of 3
#>  $ effect_size: num [1:4, 1:4] NA 0.399 0.232 0.294 0.399 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ p_value    : num [1:4, 1:4] NA 1.56e-75 1.69e-25 5.00e-41 1.56e-75 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ n          : num [1:4, 1:4] NA 2201 2201 2201 2201 ...
#>   ..- attr(*, "dimnames")=List of 2
```
