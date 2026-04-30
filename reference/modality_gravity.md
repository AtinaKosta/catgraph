# Modality Gravity Index for catmodgraph objects

Computes four prevalence-weighted structural indices for every node in a
modality-level network: gravitational mass (MGI+), orbital drag (MGI-),
net gravity (dMGI), and orbital score (OS). These indices extend
standard centrality measures by incorporating the empirical prevalence
of each modality, allowing direction-sensitive identification of
attractor and satellite modalities.

## Usage

``` r
modality_gravity(x, weight_attr = "weight", min_prevalence = 0.01)
```

## Arguments

- x:

  A `catmodgraph` object as returned by
  [`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
  or
  [`prune_modality_edges`](https://atinakosta.github.io/catgraph/reference/prune_modality_edges.md).
  The object must contain a `$data` slot with the original data frame
  from which prevalences are computed.

- weight_attr:

  Character. Name of the edge attribute to use as association weight.
  Defaults to `"weight"` (absolute phi / Cramer's V stored by
  `build_modality_graph`).

- min_prevalence:

  Numeric in (0, 1). Modalities with prevalence below this value are
  retained in the output but their prevalence-ratio contributions are
  flagged. Default `0.01` (1%).

## Value

A data frame with one row per modality node and columns:

- node:

  Full node label in `variable=modality` format.

- variable:

  Variable name (left of `=`).

- modality:

  Modality label (right of `=`).

- prevalence:

  Empirical relative frequency computed from `x$data`, ignoring `NA`s.

- degree:

  Number of edges incident to the node.

- strength:

  Weighted degree: sum of incident edge weights.

- mgi_plus:

  MGI+: sum of \\w\_{ij} \cdot p_i/p_j\\ over neighbours \\j\\ with
  \\p_j \< p_i\\.

- mgi_minus:

  MGI-: sum of \\w\_{ij} \cdot p_i/p_j\\ over neighbours \\j\\ with
  \\p_j \> p_i\\.

- delta_mgi:

  Net gravity: `mgi_plus - mgi_minus`. Positive = net attractor;
  negative = net satellite.

- mgi_plus_norm:

  MGI+ normalised to \\\[0, 1\]\\ within the graph, for cross-graph
  comparisons.

- os:

  Orbital Score: sum of \\w\_{ij} \cdot p_j/p_i\\ over all neighbours.
  High OS indicates strong gravitational capture.

- role:

  Character: `"attractor"` if `delta_mgi > 0`, `"satellite"` if
  `delta_mgi < 0`, `"neutral"` if `delta_mgi == 0` (isolated nodes or
  perfectly balanced).

Rows are ordered by `delta_mgi` descending (strongest attractors first).

## Details

**Prevalence computation.** Prevalences are computed from `x$data` as
\\n\_{ij} / n_j\\ where \\n\_{ij}\\ is the number of non-missing
observations of variable \\j\\ taking value \\i\\, and \\n_j\\ is the
total non-missing count for that variable. Each variable is normalised
independently, so prevalences sum to 1 within each variable.

**Relationship to standard centrality.** When all nodes in the graph
have equal prevalence, \\p_i/p_j = 1\\ for every edge, and
`mgi_plus + mgi_minus` equals weighted degree (strength). `delta_mgi` is
then determined solely by whether a node's neighbours are uniformly more
or less prevalent, which is uninformative. The metric is most useful
when prevalences vary substantially across modalities.

**Edge weights.** The function takes the absolute value of edge weights
before computing ratios so that signed phi values (stored for the
modality graph) do not cancel positive and negative contributions. The
sign of the association is a separate dimension captured in `std_resid`;
`delta_mgi` captures prevalence-weighted directional dominance.

## See also

[`plot_gravity`](https://atinakosta.github.io/catgraph/reference/plot_gravity.md),
[`compare_gravity`](https://atinakosta.github.io/catgraph/reference/compare_gravity.md),
[`node_centrality`](https://atinakosta.github.io/catgraph/reference/node_centrality.md),
[`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)

## Examples

``` r
data(survey_health)
mg <- build_modality_graph(survey_health)
mg <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
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

# Top attractors
head(grav[grav$role == "attractor", ], 5)
#> 
#> === Modality Gravity Index ===
#> (5 nodes: 5 attractors, 0 neutral, 0 satellites)
#> 
#> ATTRACTORS  [dMGI > 0]
#>   Node                                               Prev      MGI+      MGI-      dMGI      OS
#> --------------------------------------------------------------------------------------------- 
#>   lung_disease=no                                   0.855     1.958     0.000     1.958   0.329
#>   smoking_status=never                              0.469     1.630     0.173     1.458   0.935
#>   health_insurance=yes                              0.783     1.291     0.000     1.291   0.198
#>   age_group=18-34                                   0.379     1.502     0.316     1.186   1.595
#>   bmi_category=normal                               0.434     0.988     0.000     0.988   0.508
#> 
#> Columns: Prev = prevalence  MGI+ = gravitational mass  MGI- = orbital drag  dMGI = net gravity  OS = orbital score

# Top satellites (highest orbital score)
head(grav[order(-grav$os), ], 5)
#> 
#> === Modality Gravity Index ===
#> (5 nodes: 2 attractors, 0 neutral, 3 satellites)
#> 
#> ATTRACTORS  [dMGI > 0]
#>   Node                                               Prev      MGI+      MGI-      dMGI      OS
#> --------------------------------------------------------------------------------------------- 
#>   age_group=18-34                                   0.379     1.502     0.316     1.186   1.595
#>   smoking_status=current                            0.257     0.579     0.161     0.418   1.451
#> 
#> SATELLITES  [dMGI < 0]
#>   Node                                               Prev      MGI+      MGI-      dMGI      OS
#> --------------------------------------------------------------------------------------------- 
#>   bmi_category=obese                                0.224     0.000     0.508    -0.508   1.112
#>   age_group=55+                                     0.239     0.232     0.739    -0.507   2.403
#>   lung_disease=yes                                  0.145     0.000     0.331    -0.331   1.945
#> 
#> Columns: Prev = prevalence  MGI+ = gravitational mass  MGI- = orbital drag  dMGI = net gravity  OS = orbital score
```
