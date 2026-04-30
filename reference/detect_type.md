# Detect contingency table type for a pair of categorical variables

For a given pair of columns, this function determines whether the
cross-tabulation yields a 2x2 or a larger (n x m) contingency table,
after removing rows with missing values in either column (pairwise
deletion). This drives the dispatch between the phi coefficient and
Cramer's V in
[`effect_size`](https://atinakosta.github.io/catgraph/reference/effect_size.md).

## Usage

``` r
detect_type(x, y)
```

## Arguments

- x:

  A factor, character, or logical vector.

- y:

  A factor, character, or logical vector of the same length as `x`.

## Value

A character string: `"2x2"` if both variables have exactly two levels
observed in the pairwise-complete data, or `"RxC"` otherwise.

## Details

Level counts are computed on the pairwise-complete subset, not on the
full vector. This means that a variable with three levels in the full
data may be classified as having two levels for a specific pair if one
level is entirely missing when the other variable is observed.

## See also

[`effect_size`](https://atinakosta.github.io/catgraph/reference/effect_size.md),
[`compute_assoc`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md)

## Examples

``` r
x <- c("a", "b", "a", "b", NA)
y <- c("yes", "no", "yes", "yes", "no")
detect_type(x, y)  # "2x2"
#> [1] "2x2"

z <- c("low", "mid", "high", "low", "mid")
detect_type(x[-5], z[-5])  # "RxC"
#> [1] "RxC"
```
