# node_centrality method for catmodgraph objects

Extends
[`node_centrality`](https://atinakosta.github.io/catgraph/reference/node_centrality.md)
to accept a `catmodgraph` object (modality-level network). Returns the
standard centrality measures - strength, betweenness, closeness,
eigenvector, PageRank - augmented with the gravity indices MGI+, MGI-,
dMGI, OS, and role from
[`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md).

## Usage

``` r
# S3 method for class 'catmodgraph'
node_centrality(x, normalize = TRUE, ...)
```

## Arguments

- x:

  A `catmodgraph` object as returned by
  [`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
  or
  [`prune_modality_edges`](https://atinakosta.github.io/catgraph/reference/prune_modality_edges.md).

- normalize:

  Logical. If `TRUE` (default), the five traditional centrality columns
  are normalised to \\\[0, 1\]\\. Gravity indices are never normalised
  here (use `mgi_plus_norm` for a normalised MGI+).

- ...:

  Additional arguments passed to
  [`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md).

## Value

A data frame with one row per modality node and columns: `node`,
`variable`, `modality`, `prevalence`, `degree`, `strength`,
`w_betweenness`, `w_closeness`, `w_eigenvector`, `w_pagerank`
(traditional, optionally normalised), plus `mgi_plus`, `mgi_minus`,
`delta_mgi`, `mgi_plus_norm`, `os`, `role` (gravity indices, always on
their natural scale). Rows are ordered by `delta_mgi` descending.

## See also

[`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md),
[`node_centrality`](https://atinakosta.github.io/catgraph/reference/node_centrality.md),
[`plot_gravity`](https://atinakosta.github.io/catgraph/reference/plot_gravity.md)

## Examples

``` r
data(survey_health)
mg  <- build_modality_graph(survey_health)
mg  <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
nc  <- node_centrality(mg)
print(nc[, c("node", "strength", "w_eigenvector", "delta_mgi", "os", "role")])
#>                        node   strength w_eigenvector delta_mgi     os      role
#> 1           lung_disease=no 0.51871365  6.396867e-01    1.9579 0.3290 attractor
#> 2      smoking_status=never 0.68459219  8.532382e-01    1.4575 0.9350 attractor
#> 3      health_insurance=yes 0.33155716  4.655303e-01    1.2913 0.1978 attractor
#> 4           age_group=18-34 1.00000000  1.000000e+00    1.1862 1.5947 attractor
#> 5       bmi_category=normal 0.47252704  3.454526e-01    0.9877 0.5081 attractor
#> 6         exercise_freq=low 0.74218627  6.499834e-01    0.5800 1.0977 attractor
#> 7    smoking_status=current 0.50985388  5.122413e-01    0.4179 1.4511 attractor
#> 8                sex=female 0.14810840  9.806502e-02    0.3783 0.1358 attractor
#> 9                  sex=male 0.14810840  9.806502e-02    0.3581 0.1435 attractor
#> 10   exercise_freq=moderate 0.16641723  4.437841e-02    0.3279 0.1971 attractor
#> 11        diet_quality=good 0.28239169  1.900231e-01    0.2424 0.3981 attractor
#> 12    smoking_status=former 0.33187246  5.011307e-01    0.0993 0.5588 attractor
#> 13 bmi_category=underweight 0.00000000  8.639119e-18    0.0000 0.0000   neutral
#> 14  bmi_category=overweight 0.00000000  8.639119e-18    0.0000 0.0000   neutral
#> 15       exercise_freq=high 0.55327029  5.640982e-01   -0.0007 0.9496 satellite
#> 16          age_group=35-54 0.24690951  5.194967e-02   -0.0233 0.4181 satellite
#> 17        diet_quality=fair 0.08497285  6.333688e-03   -0.1223 0.1311 satellite
#> 18        diet_quality=poor 0.31516366  1.752797e-01   -0.1259 0.5244 satellite
#> 19         lung_disease=yes 0.51871365  6.396867e-01   -0.3311 1.9454 satellite
#> 20      health_insurance=no 0.33155716  4.655303e-01   -0.3578 0.7139 satellite
#> 21            age_group=55+ 0.94189376  7.958089e-01   -0.5073 2.4026 satellite
#> 22       bmi_category=obese 0.50272921  3.178106e-01   -0.5081 1.1120 satellite
```
