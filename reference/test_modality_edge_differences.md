# Edge-wise post-hoc test for modality-network differences

After a significant global test with
[`test_modality_graph_equality`](https://atinakosta.github.io/catgraph/reference/test_modality_graph_equality.md),
this function identifies which specific edges contribute most strongly
to the observed difference. For each edge, an empirical two-sided
p-value is computed from label-permutation distributions of the absolute
edge-weight difference, and Benjamini-Hochberg FDR correction is applied
across all tested edges.

## Usage

``` r
test_modality_edge_differences(
  x,
  y,
  n_perm = 500L,
  edges = c("all", "union"),
  p_adjust = "BH",
  strata = NULL,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- x, y:

  Two `catmodgraph` objects on the same variable set, each carrying
  their `$data` component.

- n_perm:

  Integer, number of permutations. Default `500`.

- edges:

  Character, which edges to test. One of `"all"` (default, every
  upper-triangle pair of modalities) or `"union"` (only edges present in
  at least one of the two input graphs). Using `"union"` reduces the
  multiple-testing burden.

- p_adjust:

  Character, adjustment method for multiple testing, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default `"BH"`.

- strata:

  Optional stratification vector of length equal to the combined sample
  size. See
  [`test_modality_graph_equality`](https://atinakosta.github.io/catgraph/reference/test_modality_graph_equality.md).

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Logical. If `TRUE` (default), prints progress.

## Value

An object of class `catmodedgetest` with components:

- `edge_table`:

  Data frame with one row per tested edge: columns `from`, `to`,
  `weight_x`, `weight_y`, `obs_diff`, `p_empirical`, `p_adjusted`.

- `n_perm`:

  Integer number of permutations.

- `n_x`, `n_y`:

  Sample sizes.

- `edges`:

  Character, which edge-subset criterion was used.

- `p_adjust_method`:

  Adjustment method.

- `strata_used`:

  Logical.

## Details

For each edge, the permutation null is that the observed edge-weight
difference is no larger than expected under random reassignment of
sample labels. The test statistic is the observed edge-weight
difference; the null distribution is obtained by recomputing that
difference across `n_perm` random label reassignments of the combined
data. Edges not present in the pruned graphs are treated as having
weight 0.

The function is designed to be called *after* a significant global test.
Calling it without a significant omnibus inflates the family-wise Type I
error beyond the nominal FDR level, because the omnibus test acts as a
gatekeeper under the closed-testing principle.

## See also

[`test_modality_graph_equality`](https://atinakosta.github.io/catgraph/reference/test_modality_graph_equality.md)
for the global test.

## Examples

``` r
# \donttest{
data(survey_health)
df_f <- subset(survey_health, sex == "female")[, -1]
df_m <- subset(survey_health, sex == "male")[, -1]

mg_f <- build_modality_graph(df_f)
mg_m <- build_modality_graph(df_m)

edge_test <- test_modality_edge_differences(
  mg_f, mg_m, n_perm = 200, edges = "union",
  seed = 1, verbose = FALSE
)
head(edge_test$edge_table)
#>                       from                   to    weight_x     weight_y
#> 74       diet_quality=poor   exercise_freq=high 0.031381364 0.1950809523
#> 24           age_group=55+    diet_quality=poor 0.093745693 0.1986175040
#> 89       exercise_freq=low  health_insurance=no 0.207140819 0.0009701169
#> 102      exercise_freq=low health_insurance=yes 0.207140819 0.0009701169
#> 17         age_group=35-54   bmi_category=obese 0.006335803 0.1114963387
#> 90  exercise_freq=moderate  health_insurance=no 0.200021474 0.0156441558
#>       obs_diff p_empirical p_adjusted  abs_diff
#> 74  -0.1636996 0.004975124  0.4228856 0.1636996
#> 24  -0.1048718 0.004975124  0.4228856 0.1048718
#> 89   0.2061707 0.014925373  0.6343284 0.2061707
#> 102  0.2061707 0.014925373  0.6343284 0.2061707
#> 17  -0.1051605 0.019900498  0.6766169 0.1051605
#> 90   0.1843773 0.049751244  0.7611940 0.1843773
# }
```
