# Extract pairwise association weights as a matrix or tidy data frame

Returns the pairwise effect size values from a `catgraph` object as a
symmetric matrix or a long-format tidy data frame suitable for further
analysis or export.

## Usage

``` r
assoc_matrix(
  x,
  format = c("matrix", "tidy"),
  include_p = TRUE,
  include_n = TRUE
)
```

## Arguments

- x:

  A `catgraph` object.

- format:

  Character. Output format: `"matrix"` (default) returns a symmetric
  numeric matrix with `NA` on the diagonal; `"tidy"` returns a data
  frame with one row per pair.

- include_p:

  Logical. When `format = "tidy"`, whether to include the p-value
  column. Default `TRUE`.

- include_n:

  Logical. When `format = "tidy"`, whether to include the pairwise
  observation count. Default `TRUE`.

## Value

- `format = "matrix"`:

  A symmetric numeric matrix of dimension p x p (p = number of
  variables). Diagonal is `NA`. Row and column names are the variable
  names.

- `format = "tidy"`:

  A data frame with columns `var1`, `var2`, `effect_size`, `metric`,
  `type`, and optionally `p_value` and `n`.

## See also

[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`prune_edges`](https://atinakosta.github.io/catgraph/reference/prune_edges.md)

## Examples

``` r
df <- as.data.frame(Titanic)
df_exp <- df[rep(seq_len(nrow(df)), df$Freq), -5]
cg <- catgraph(df_exp)
assoc_matrix(cg)
#>              Class       Sex        Age   Survived
#> Class           NA 0.3987227 0.23194779 0.29412010
#> Sex      0.3987227        NA 0.11101269 0.45560478
#> Age      0.2319478 0.1110127         NA 0.09757511
#> Survived 0.2941201 0.4556048 0.09757511         NA
assoc_matrix(cg, format = "tidy")
#>    var1     var2 effect_size    metric type       p_value    n
#> 5   Sex Survived  0.45560478       phi  2x2 2.302151e-101 2201
#> 1 Class      Sex  0.39872269 cramers_v  RxC  1.556637e-75 2201
#> 3 Class Survived  0.29412010 cramers_v  RxC  4.999928e-41 2201
#> 2 Class      Age  0.23194779 cramers_v  RxC  1.694884e-25 2201
#> 4   Sex      Age  0.11101269       phi  2x2  1.907432e-07 2201
#> 6   Age Survived  0.09757511       phi  2x2  4.700752e-06 2201
```
