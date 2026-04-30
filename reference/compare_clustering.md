# Compare all weighted clustering coefficients side by side

A convenience wrapper around
[`clustering_coef`](https://atinakosta.github.io/catgraph/reference/clustering_coef.md)
that always computes all five methods and returns them in a single wide
data frame, making it easy to identify variables that score consistently
high (robust cluster hubs) versus those that vary across methods
(structurally ambiguous variables).

## Usage

``` r
compare_clustering(x, normalize = TRUE)
```

## Arguments

- x:

  A `catgraph` object.

- normalize:

  Logical. Default `TRUE`.

## Value

A `data.frame` with columns `variable`, `watts_strogatz`, `barrat`,
`onnela`, `zhang`, `redundancy`, and `mean_cc` (the row mean across all
five methods, excluding `NA`s). Sorted by `mean_cc` descending.

## See also

[`clustering_coef`](https://atinakosta.github.io/catgraph/reference/clustering_coef.md),
[`plot_clustering`](https://atinakosta.github.io/catgraph/reference/plot_clustering.md)

## Examples

``` r
df <- expand_table(Titanic)
cg <- catgraph(df)
compare_clustering(cg)
#>   variable watts_strogatz barrat    onnela     zhang redundancy   mean_cc
#> 1    Class              1    0.5 0.2861019 0.2530913  1.0000000 0.6078387
#> 2      Sex              1    0.5 0.2795694 0.2482151  0.8088828 0.5673334
#> 3      Age              1    0.5 0.2106029 0.3691475  0.7145339 0.5588569
#> 4 Survived              1    0.5 0.2688767 0.3138760  0.6472822 0.5460070
```
