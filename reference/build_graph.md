# Build the underlying igraph association network

Computes pairwise effect sizes (phi or Cramer's V) for all pairs of
categorical columns in a data frame and returns the underlying igraph
object. This is the lower-level computational engine used by
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md).
Most users should call
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
unless they specifically need direct access to the raw igraph
representation.

## Usage

``` r
build_graph(
  data,
  corrected = FALSE,
  correct = FALSE,
  simulate_p = FALSE,
  B = 2000L
)
```

## Arguments

- data:

  A data frame or tibble. All columns are treated as categorical.
  Non-factor, non-character, non-logical columns are coerced to
  character with a message. Columns with only one unique observed value
  (after pairwise deletion) are dropped with a warning.

- corrected:

  Logical. Whether to use bias-corrected Cramer's V / phi (Bergsma,
  2013). Default `FALSE`.

- correct:

  Logical. Yates' continuity correction for chi-square. Default `FALSE`.

- simulate_p:

  Logical. Use Monte Carlo simulation for p-values. Default `FALSE`.

- B:

  Integer. Number of Monte Carlo resamples. Default `2000L`.

## Value

An igraph undirected graph. Pairs with a true zero effect size (no
association whatsoever) are represented as *absent edges* rather than
near-zero edges, so the graph is sparse rather than structurally
complete. The attribute `"processed_data"` on the returned graph holds
the data frame actually used for estimation (after coercion and
constant-column removal), which downstream functions such as
[`catgraph_ci`](https://atinakosta.github.io/catgraph/reference/catgraph_ci.md)
use when resampling. The graph attribute `"pair_results"` stores the
full pairwise results table before zero-weight edges are omitted.

For ordinary package use, prefer
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
which wraps this graph together with processed data, metadata, and S3
methods.

Vertex and edge attributes:

- Vertices:

  One per column in `data`, with vertex attribute `name` set to the
  column name. Isolated vertices are preserved even if all their pairs
  have zero effect size.

- Edge attribute `weight`:

  The phi or Cramer's V value.

- Edge attribute `metric`:

  `"phi"` or `"cramers_v"`.

- Edge attribute `corrected`:

  Whether bias correction was applied.

- Edge attribute `p_value`:

  Chi-square p-value.

- Edge attribute `statistic`:

  Chi-square statistic.

- Edge attribute `df`:

  Degrees of freedom.

- Edge attribute `n`:

  Pairwise-complete observation count for that pair. Values can differ
  across edges when missingness is present (pairwise deletion).

- Edge attribute `type`:

  `"2x2"` or `"RxC"`.

- Edge attribute `estimable`:

  Logical indicating whether the pairwise effect size was estimable
  before zero-weight omission.

## Details

Computes effect sizes (phi or Cramer's V) for all pairs of categorical
columns in a data frame and returns an igraph object whose edge weights
correspond to those effect sizes. This is the main computational engine
used by the top-level
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
constructor.

**Scope.** The returned graph represents *pairwise marginal* association
strength. It is not a conditional-independence graphical model: an edge
between `A` and `B` does not imply that the variables remain dependent
after controlling for the other variables in the data. See the package
vignette section "Scope and interpretation" for details.

All variable pairs with non-zero effect size are included by default.
Use
[`prune_edges`](https://atinakosta.github.io/catgraph/reference/prune_edges.md)
to remove edges below a weight or adjusted-p threshold after
construction.

**Note on zero-weight pairs.** In earlier versions of the package (\<=
0.3.0), zero-weight pairs were stored as edges with weight
`.Machine$double.eps` to guarantee a fully connected graph. This made
the graph structurally complete and silently inflated density-based
measures. From 0.4.0 onwards, zero-weight pairs are absent edges; a
dense similarity matrix suitable for heatmaps is available separately
via
[`assoc_similarity`](https://atinakosta.github.io/catgraph/reference/assoc_similarity.md)
or
[`assoc_matrix`](https://atinakosta.github.io/catgraph/reference/assoc_matrix.md).

## References

Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
*Journal of the Korean Statistical Society*, 42(3), 323–328.
[doi:10.1016/j.jkss.2012.10.002](https://doi.org/10.1016/j.jkss.2012.10.002)

Csardi, G., & Nepusz, T. (2006). The igraph software package for complex
network research. *InterJournal, Complex Systems*, 1695.
<https://igraph.org>

## See also

[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`prune_edges`](https://atinakosta.github.io/catgraph/reference/prune_edges.md),
[`effect_size`](https://atinakosta.github.io/catgraph/reference/effect_size.md),
[`assoc_similarity`](https://atinakosta.github.io/catgraph/reference/assoc_similarity.md)

## Examples

``` r
data(HairEyeColor)
df <- expand_table(HairEyeColor)
g  <- build_graph(df[, c("Hair", "Eye")])
igraph::E(g)$weight
#> [1] 0.2790446
```
