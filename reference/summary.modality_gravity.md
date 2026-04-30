# Summary method for modality_gravity output

Returns a compact role-distribution table and per-variable breakdown,
and prints a readable synopsis to the console.

## Usage

``` r
# S3 method for class 'modality_gravity'
summary(object, ...)
```

## Arguments

- object:

  A data frame returned by
  [`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md).

- ...:

  Ignored.

## Value

Invisibly returns a list with elements:

- role_counts:

  Named integer vector of attractor/neutral/satellite counts.

- by_variable:

  Data frame with one row per variable showing the dominant role and
  mean `delta_mgi` across its modalities.

- top_attractor:

  Single-row data frame: modality with highest `delta_mgi`.

- top_satellite:

  Single-row data frame: modality with lowest `delta_mgi`.

- spearman_strength:

  Spearman correlation between `strength` and `delta_mgi`.

## See also

[`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md),
[`print.modality_gravity`](https://atinakosta.github.io/catgraph/reference/print.modality_gravity.md)

## Examples

``` r
data(survey_health)
mg   <- build_modality_graph(survey_health)
mg   <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
grav <- modality_gravity(mg)
summary(grav)
#> 
#> === Modality Gravity - Summary ===
#> 
#> Role distribution:
#> 
#> attractor   neutral satellite 
#>        12         2         8 
#> 
#> Dominant attractor : lung_disease=no  (dMGI = 1.958, prev = 0.855)
#> Strongest satellite: bmi_category=obese  (dMGI = -0.508, prev = 0.224)
#> 
#> Spearman rho (strength vs dMGI): 0.191
#> (Values far from +/-1.0 indicate MGI captures structure beyond connectivity)
#> 
#> Per-variable breakdown:
#>          variable n_modalities mean_delta dominant_role
#>      lung_disease            2      0.813     attractor
#>    smoking_status            3      0.658     attractor
#>  health_insurance            2      0.467     attractor
#>               sex            2      0.368     attractor
#>     exercise_freq            3      0.302     attractor
#>         age_group            3      0.219     satellite
#>      bmi_category            4      0.120       neutral
#>      diet_quality            3     -0.002     satellite
```
