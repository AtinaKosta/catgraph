# Add bootstrap confidence intervals to all edges of a catgraph

Calls
[`bootstrap_ci`](https://atinakosta.github.io/catgraph/reference/bootstrap_ci.md)
for every edge in a `catgraph` object and stores the lower and upper
bounds as additional edge attributes (`ci_lower`, `ci_upper`, `ci_conf`,
`ci_type`).

## Usage

``` r
catgraph_ci(
  x,
  R = 1000L,
  conf = 0.95,
  type = c("percentile", "bca"),
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  A `catgraph` object.

- R:

  Integer. Number of bootstrap resamples per pair. Default `1000L`.

- conf:

  Numeric. Confidence level. Default `0.95`.

- type:

  Character. `"percentile"` or `"bca"`. Default `"percentile"`.

- seed:

  Integer or `NULL`. Base seed; each pair uses `seed + i` to ensure
  per-pair reproducibility without global state. Default `NULL`.

- verbose:

  Logical. Print a progress counter. Default `TRUE`.

## Value

The input `catgraph` object with three new edge attributes: `ci_lower`,
`ci_upper`, `ci_conf`, `ci_type`.

## Details

For large graphs (many variable pairs) this function can be slow because
it runs `R` resamples per edge. Consider lowering `R` or running on a
pruned graph (via
[`prune_edges`](https://atinakosta.github.io/catgraph/reference/prune_edges.md))
to reduce computation.

## See also

[`bootstrap_ci`](https://atinakosta.github.io/catgraph/reference/bootstrap_ci.md),
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md)

## Examples

``` r
df <- as.data.frame(Titanic)
df_exp <- df[rep(seq_len(nrow(df)), df$Freq), -5]
cg <- catgraph(df_exp)
# \donttest{
cg <- catgraph_ci(cg, R = 500, seed = 1)
#>   Bootstrap CI: edge 1 / 6  Bootstrap CI: edge 2 / 6  Bootstrap CI: edge 3 / 6  Bootstrap CI: edge 4 / 6  Bootstrap CI: edge 5 / 6  Bootstrap CI: edge 6 / 6
igraph::E(cg$graph)$ci_lower
#> [1] 0.37158251 0.20622301 0.25051349 0.06191884 0.41857071 0.05197030
igraph::E(cg$graph)$ci_upper
#> [1] 0.4286221 0.2609237 0.3356203 0.1651256 0.4913456 0.1415683
# }
```
