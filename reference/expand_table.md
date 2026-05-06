# Expand a contingency table or frequency data frame to observation-level format

Converts multi-dimensional contingency tables (`table`, `array`,
`ftable`), or data frames that contain a frequency/count column, into a
flat data frame where every row represents one observation. This is the
required input format for
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md).

## Usage

``` r
expand_table(tbl, freq_col = NULL, as_factor = TRUE, drop_zero = TRUE)
```

## Arguments

- tbl:

  A `table`, `array`, `ftable`, or `data.frame`. All standard R
  contingency table formats are supported, including multi-dimensional
  arrays such as `Titanic` (4-D) and `HairEyeColor` (3-D).

- freq_col:

  Character or integer. Only used when `tbl` is a `data.frame`. The name
  or column index of the frequency/count column. If `NULL` (default),
  the function looks for a column named `"Freq"`, `"freq"`, `"n"`,
  `"count"`, or `"Count"` in that order. An error is raised if none is
  found.

- as_factor:

  Logical. If `TRUE` (default), all resulting columns are coerced to
  factors, preserving the level order from the original object's
  `dimnames`. Set to `FALSE` to return character columns.

- drop_zero:

  Logical. If `TRUE` (default), rows corresponding to zero-count cells
  are silently dropped. Set to `FALSE` to keep them (those rows will
  appear zero times and therefore not affect results, but the factor
  levels will still be present).

## Value

A `data.frame` with one row per observation and one column per
categorical variable. The column names are taken from
[`dimnames()`](https://rdrr.io/r/base/dimnames.html) of the input
object, or from the non-frequency columns of the input data frame. Row
names are reset to `NULL`.

## Details

**Accepted input formats:**

- One row per observation:

  Already in the correct format — pass directly to
  [`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
  without calling `expand_table()`.

- `table` / `array`:

  The standard output of [`table()`](https://rdrr.io/r/base/table.html),
  [`xtabs()`](https://rdrr.io/r/stats/xtabs.html), or built-in datasets
  such as `Titanic`, `HairEyeColor`, and `UCBAdmissions`.

- `ftable`:

  Converted to a `table` internally before expansion.

- `data.frame` with frequency column:

  The output of `as.data.frame(some_table)`, which always contains a
  `Freq` column.

**Not accepted:**

- Raw numeric matrices:

  A plain matrix of counts without `dimnames` cannot be safely
  converted. Assign `dimnames` first.

- Numeric 0/1 columns:

  Columns coded as integers or doubles are not treated as categorical.
  Coerce with [`as.factor()`](https://rdrr.io/r/base/factor.html) or
  [`as.character()`](https://rdrr.io/r/base/character.html) first, or
  pass them through
  [`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
  which will coerce and warn automatically.

## References

R Core Team (2024). *R: A Language and Environment for Statistical
Computing*. R Foundation for Statistical Computing, Vienna, Austria.
<https://www.R-project.org/>

## See also

[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`assoc_matrix`](https://atinakosta.github.io/catgraph/reference/assoc_matrix.md)

## Examples

``` r
# Built-in 4-D table
df <- expand_table(Titanic)
str(df)
#> 'data.frame':    2201 obs. of  4 variables:
#>  $ Class   : Factor w/ 4 levels "1st","2nd","3rd",..: 3 3 3 3 3 3 3 3 3 3 ...
#>  $ Sex     : Factor w/ 2 levels "Female","Male": 2 2 2 2 2 2 2 2 2 2 ...
#>  $ Age     : Factor w/ 2 levels "Adult","Child": 2 2 2 2 2 2 2 2 2 2 ...
#>  $ Survived: Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 1 1 1 ...
nrow(df)  # 2201 passengers
#> [1] 2201

# Built-in 3-D table
df2 <- expand_table(HairEyeColor)
str(df2)
#> 'data.frame':    592 obs. of  3 variables:
#>  $ Hair: Factor w/ 4 levels "Black","Blond",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Eye : Factor w/ 4 levels "Blue","Brown",..: 2 2 2 2 2 2 2 2 2 2 ...
#>  $ Sex : Factor w/ 2 levels "Female","Male": 2 2 2 2 2 2 2 2 2 2 ...

# data.frame with Freq column (output of as.data.frame on a table)
tab_df <- as.data.frame(UCBAdmissions)
df3 <- expand_table(tab_df)
nrow(df3)  # 4526 applicants
#> [1] 4526

# Custom data frame with a count column
survey <- data.frame(
  gender   = c("M", "F", "M", "F"),
  smokes   = c("yes", "yes", "no", "no"),
  n        = c(23L, 15L, 48L, 61L)
)
df4 <- expand_table(survey, freq_col = "n")
nrow(df4)  # 147 observations
#> [1] 147

# Use directly with catgraph
cg <- catgraph(expand_table(Titanic))
cg
#> catgraph object (pairwise association network)
#>   Variables : 4 
#>   Edges     : 6 
#>   Method    : Cramer's V (classical) 
#>   Weights   : min = 0.0976  median = 0.2630  max = 0.4556
#>   Note      : edges encode pairwise marginal association, not
#>               conditional independence. All metrics lie on [0, 1].
#>               NMI / AMI weights are not exchangeable with Cramer's V
#>               weights across graph objects. See vignette
#>               'Methodological caveats'.
```
