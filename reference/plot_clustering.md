# Plot weighted clustering coefficients for a catgraph

Produces a bar chart (single method) or a grouped/faceted comparison
chart (multiple methods) of clustering coefficients.

## Usage

``` r
plot_clustering(
  x,
  method = "barrat",
  normalize = TRUE,
  title = NULL,
  engine = c("ggplot2", "base")
)
```

## Arguments

- x:

  A `catgraph` object.

- method:

  Character. One method name or `"all"`. Default `"barrat"`.

- normalize:

  Logical. Default `TRUE`.

- title:

  Character. Plot title. Default `NULL`.

- engine:

  Character. `"ggplot2"` (default) or `"base"`.

## Value

For `engine = "ggplot2"`: a `ggplot` object. For `engine = "base"`:
`NULL`, invisibly.

## See also

[`clustering_coef`](https://atinakosta.github.io/catgraph/reference/clustering_coef.md),
[`compare_clustering`](https://atinakosta.github.io/catgraph/reference/compare_clustering.md)

## Examples

``` r
df <- expand_table(Titanic)
cg <- catgraph(df)
plot_clustering(cg)

plot_clustering(cg, method = "all")

```
