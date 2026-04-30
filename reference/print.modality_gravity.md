# Print method for modality_gravity output

Displays a formatted, role-grouped summary of the data frame returned by
[`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md).
Attractors are shown first in descending `delta_mgi` order, then
neutrals, then satellites in ascending order (strongest satellite last).

## Usage

``` r
# S3 method for class 'modality_gravity'
print(x, digits = 3L, max_width = 45L, ...)
```

## Arguments

- x:

  A data frame returned by
  [`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md).

- digits:

  Integer. Number of decimal places for numeric columns. Default `3L`.

- max_width:

  Integer. Maximum characters for the node label column before
  truncation with `...`. Default `45L`.

- ...:

  Ignored. Present for S3 compatibility.

## Value

Invisibly returns `x`.

## See also

[`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md),
[`summary.modality_gravity`](https://atinakosta.github.io/catgraph/reference/summary.modality_gravity.md)

## Examples

``` r
data(survey_health)
mg   <- build_modality_graph(survey_health)
mg   <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
grav <- modality_gravity(mg)
print(grav)
#> 
#> === Modality Gravity Index ===
#> (22 nodes: 12 attractors, 2 neutral, 8 satellites)
#> 
#> ATTRACTORS  [dMGI > 0]
#>   Node                                               Prev      MGI+      MGI-      dMGI      OS
#> --------------------------------------------------------------------------------------------- 
#>   lung_disease=no                                   0.855     1.958     0.000     1.958   0.329
#>   smoking_status=never                              0.469     1.630     0.173     1.458   0.935
#>   health_insurance=yes                              0.783     1.291     0.000     1.291   0.198
#>   age_group=18-34                                   0.379     1.502     0.316     1.186   1.595
#>   bmi_category=normal                               0.434     0.988     0.000     0.988   0.508
#>   exercise_freq=low                                 0.352     0.940     0.360     0.580   1.098
#>   smoking_status=current                            0.257     0.579     0.161     0.418   1.451
#>   sex=female                                        0.514     0.378     0.000     0.378   0.136
#>   sex=male                                          0.486     0.358     0.000     0.358   0.144
#>   exercise_freq=moderate                            0.364     0.328     0.000     0.328   0.197
#>   diet_quality=good                                 0.337     0.362     0.119     0.242   0.398
#>   smoking_status=former                             0.274     0.280     0.181     0.099   0.559
#> 
#> NEUTRAL     [dMGI = 0]
#>   Node                                               Prev      MGI+      MGI-      dMGI      OS
#> --------------------------------------------------------------------------------------------- 
#>   bmi_category=underweight                          0.022     0.000     0.000     0.000   0.000
#>   bmi_category=overweight                           0.319     0.000     0.000     0.000   0.000
#> 
#> SATELLITES  [dMGI < 0]
#>   Node                                               Prev      MGI+      MGI-      dMGI      OS
#> --------------------------------------------------------------------------------------------- 
#>   bmi_category=obese                                0.224     0.000     0.508    -0.508   1.112
#>   age_group=55+                                     0.239     0.232     0.739    -0.507   2.403
#>   health_insurance=no                               0.217     0.000     0.358    -0.358   0.714
#>   lung_disease=yes                                  0.145     0.000     0.331    -0.331   1.945
#>   diet_quality=poor                                 0.312     0.160     0.286    -0.126   0.524
#>   diet_quality=fair                                 0.352     0.000     0.122    -0.122   0.131
#>   age_group=35-54                                   0.382     0.159     0.182    -0.023   0.418
#>   exercise_freq=high                                0.284     0.381     0.382    -0.001   0.950
#> 
#> Columns: Prev = prevalence  MGI+ = gravitational mass  MGI- = orbital drag  dMGI = net gravity  OS = orbital score
```
