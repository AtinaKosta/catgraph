# Plot a jointbalance diagnostic

Two-panel layout: the left panel is a bar chart of adjusted marginal
p-values per variable (the "Table 1" view); the right panel is the
modality-difference graph for one selected group pair (from
[`plot_modality_difference`](https://atinakosta.github.io/catgraph/reference/plot_modality_difference.md)).
For \\k \> 2\\ groups the default pair is the one with the smallest
Bonferroni-adjusted omnibus p-value.

## Usage

``` r
# S3 method for class 'jointbalance'
plot(x, pair = NULL, layout = c("side_by_side", "marginal_only"), ...)
```

## Arguments

- x:

  A `jointbalance` object.

- pair:

  Optional length-2 character vector specifying which group pair to
  visualise. Default `NULL` picks the most significant pair.

- layout:

  Character. One of `"side_by_side"` (default) or `"marginal_only"`. Use
  `"marginal_only"` when no pairwise omnibus rejected and there is no
  edge-wise result to plot.

- ...:

  Further arguments passed to
  [`plot_modality_difference`](https://atinakosta.github.io/catgraph/reference/plot_modality_difference.md).

## Value

Invisibly returns `x`.
