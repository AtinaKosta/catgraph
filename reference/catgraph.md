# Construct a categorical association network

The primary user-facing constructor for catgraph. It computes pairwise
effect sizes (phi or Cramer's V) for all categorical variable pairs,
stores the resulting weighted igraph network, and preserves processed
data and metadata for downstream analysis. Use this function for
standard workflows; use
[`build_graph`](https://atinakosta.github.io/catgraph/reference/build_graph.md)
only when a raw igraph object is required.

## Usage

``` r
catgraph(
  data,
  corrected = FALSE,
  correct = FALSE,
  simulate_p = FALSE,
  B = 2000L
)

# S3 method for class 'catgraph'
print(x, ...)

# S3 method for class 'catgraph'
summary(object, top = 10L, ...)
```

## Arguments

- data:

  A data frame or tibble whose columns represent categorical variables.
  Factor, character, and logical columns are supported. Numeric columns
  are coerced to character with a message.

- corrected:

  Logical. If `FALSE` (default), classical phi and Cramer's V are
  computed. If `TRUE`, the bias-corrected estimators of Bergsma (2013)
  are used.

- correct:

  Logical. Yates' continuity correction for the chi-square test. Default
  `FALSE`.

- simulate_p:

  Logical. Monte Carlo p-value simulation. Default `FALSE`.

- B:

  Integer. Monte Carlo resamples when `simulate_p = TRUE`. Default
  `2000L`.

- x:

  A `catgraph` object.

- ...:

  Ignored.

- object:

  A `catgraph` object.

- top:

  Integer. Number of strongest edges to display. Use `Inf` for all
  edges. Default `10L`.

## Value

An S3 object of class `catgraph` containing:

- `graph`:

  An undirected weighted igraph object. True zero associations are
  absent edges, not near-zero edges.

- `data`:

  The **processed** data frame actually used for estimation (after
  non-categorical coercion and constant-column removal). Downstream
  functions such as
  [`catgraph_ci`](https://atinakosta.github.io/catgraph/reference/catgraph_ci.md)
  resample from this object. Changed from `raw_data` in v0.4.0 to fix an
  internal-consistency bug.

- `raw_data`:

  The original input data frame, for reference.

- `corrected`:

  Logical flag indicating which estimator is active.

- `n_vars`:

  Number of variables (graph vertices).

- `n_pairs_total`:

  Number of variable pairs evaluated.

- `n_pairs`:

  Number of retained graph edges (pairs with non-zero effect size).

- `call`:

  The matched call.

## Details

**Scope.** A `catgraph` is a *pairwise association network*, not a
conditional-independence graphical model. Edges encode bivariate
dependence between two variables and do not imply that the two variables
remain dependent after controlling for the remaining variables.
Interpret centrality, community, and bridge measures accordingly. See
the package vignette for a full discussion.

All variable pairs with non-zero effect size are retained by default (no
thresholding at construction time). To remove weak or non-significant
edges, pass the object to
[`prune_edges`](https://atinakosta.github.io/catgraph/reference/prune_edges.md).

## Methods (by generic)

- `print(catgraph)`: Print a concise summary of a `catgraph` object.

- `summary(catgraph)`: Summarise a `catgraph` object, listing edges
  sorted by effect size.

## References

Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
*Journal of the Korean Statistical Society*, 42(3), 323–328.
[doi:10.1016/j.jkss.2012.10.002](https://doi.org/10.1016/j.jkss.2012.10.002)

## See also

[`prune_edges`](https://atinakosta.github.io/catgraph/reference/prune_edges.md),
[`detect_clusters`](https://atinakosta.github.io/catgraph/reference/detect_clusters.md),
[`plot.catgraph`](https://atinakosta.github.io/catgraph/reference/plot.catgraph.md),
[`assoc_matrix`](https://atinakosta.github.io/catgraph/reference/assoc_matrix.md),
[`assoc_similarity`](https://atinakosta.github.io/catgraph/reference/assoc_similarity.md)

## Examples

``` r
df <- expand_table(Titanic)
cg <- catgraph(df)
cg
#> catgraph object (pairwise association network)
#>   Variables : 4 
#>   Edges     : 6 
#>   Estimator : classical 
#>   Weights   : min = 0.0976  median = 0.2630  max = 0.4556
#>   Metric mix: cramers_v = 3, phi = 3 
#>   Note      : edges encode pairwise marginal association, not
#>               conditional independence. Edge weights use phi
#>               (2x2) and Cramer's V (RxC); both lie on [0, 1],
#>               but are not strictly exchangeable across table
#>               dimensions. Interpret mixed-metric graphs with care.
#>               See vignette 'Methodological caveats', item 2.
summary(cg)
#> catgraph summary
#>   Variables       : 4 
#>   Pairs evaluated : 6 
#>   Edges retained  : 6 
#> 
#>   Estimator       : classical 
#> 
#>   Top 6 edges by effect size:
#> 
#>    var1     var2 effect_size    metric    p_value    n type
#> 1   Sex Survived     0.45560       phi 2.302e-101 2201  2x2
#> 2 Class      Sex     0.39872 cramers_v  1.557e-75 2201  RxC
#> 3 Class Survived     0.29412 cramers_v  5.000e-41 2201  RxC
#> 4 Class      Age     0.23195 cramers_v  1.695e-25 2201  RxC
#> 5   Sex      Age     0.11101       phi  1.907e-07 2201  2x2
#> 6   Age Survived     0.09758       phi  4.701e-06 2201  2x2

cg_bc <- catgraph(df, corrected = TRUE)
```
