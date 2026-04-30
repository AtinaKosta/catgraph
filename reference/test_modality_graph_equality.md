# Permutation test for equality of modality-graph structure

Performs an omnibus permutation test of whether two samples show the
same modality-level association structure.

## Usage

``` r
test_modality_graph_equality(
  x,
  y,
  n_perm = 500L,
  statistic = c("frobenius", "jaccard", "max"),
  test_type = c("unfiltered", "pipeline"),
  min_weight = 0.1,
  max_p = 0.05,
  strata = NULL,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- x, y:

  Two `catmodgraph` objects constructed from datasets with the same
  variable set (same column names). The underlying row data must be
  stored in each object's `$data` component, which is standard output
  from
  [`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md).

- n_perm:

  Integer. Number of permutations. Default `500`. For publication-grade
  p-values use `2000` or more.

- statistic:

  Character. Test statistic, one of `"frobenius"` (default),
  `"jaccard"`, or `"max"`. See Details.

- test_type:

  Character. How to compute the test statistic. One of `"unfiltered"`
  (default) or `"pipeline"`. See Details.

- min_weight:

  Numeric. Minimum edge-weight threshold used only when
  `test_type = "pipeline"`. Default `0.10`.

- max_p:

  Numeric. Maximum edge p-value used only when `test_type = "pipeline"`.
  Default `0.05`.

- strata:

  Optional vector of length equal to the combined sample size. If
  supplied, permutations are conducted *within* levels of `strata`,
  preserving the joint distribution of the stratification variable under
  the null. Use to remove the confounding of study membership with a
  known nuisance variable (e.g., fieldwork year, geographic region).

- seed:

  Optional integer seed for reproducibility.

- verbose:

  Logical. If `TRUE` (default), prints progress.

## Value

An object of class `catmodtest` with components:

- `statistic`:

  Character name of the test statistic used.

- `observed`:

  Numeric, the observed value of the test statistic on the two input
  graphs.

- `null_distribution`:

  Numeric vector of length `n_perm` with the test statistic computed
  under each label permutation.

- `p_value`:

  Numeric empirical p-value, `(sum(null >= obs) + 1) / (n_perm + 1)`.

- `n_perm`:

  Integer number of permutations.

- `n_x`, `n_y`:

  Sample sizes of the two inputs.

- `test_type`:

  Character, which pipeline mode was used.

- `strata_used`:

  Logical indicator of whether stratified permutation was applied.

## Details

The test evaluates whether the observed difference between two
modality-level association matrices is larger than expected under random
reassignment of sample labels. Rejection supports a difference in the
overall marginal association structure. The test is omnibus: it does not
identify which specific modality pairs drive the difference. Use
[`test_modality_edge_differences`](https://atinakosta.github.io/catgraph/reference/test_modality_edge_differences.md)
for edge-wise follow-up.

**Test statistics.** Let \\A\\ and \\B\\ denote the weighted adjacency
matrices of the two graphs on a common node set.

- `"frobenius"` (default) uses \\\\A - B\\\_F^2 = \sum\_{i,j} (A\_{ij} -
  B\_{ij})^2\\, sensitive to all edge-weight changes and dominated by
  high-weight edges.

- `"jaccard"` uses \\1 - \|E_A \cap E_B\|/\|E_A \cup E_B\|\\, the
  complement of edge-set agreement. Topological, ignores weight
  magnitudes.

- `"max"` uses \\\max\_{i,j} \|A\_{ij} - B\_{ij}\|\\, sensitive to any
  single sharp edge-weight change.

**Pipeline modes.**

- `"unfiltered"` (default): the test statistic is computed on the *full
  unpruned* phi matrix from each sample (i.e., every pair of modalities
  contributes its raw phi). This tests the joint distribution cleanly,
  with no confounding between "edges differ" and "different edges
  survived pruning."

- `"pipeline"`: each permutation re-runs
  [`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
  plus
  [`prune_modality_edges`](https://atinakosta.github.io/catgraph/reference/prune_modality_edges.md)
  and compares the pruned adjacency matrices. Matches the actual
  analysis a user ran but the resulting null mixes edge-weight and
  edge-set changes. Slower.

**Stratification.** When `strata` is supplied, permutations rearrange
study labels *within* strata, so the joint distribution of the
stratification variable is preserved under the null. This evaluates
sample-label differences conditional on the supplied strata and should
be interpreted as a stratified permutation analysis, not as causal
adjustment. **Assumptions.** Respondents are assumed i.i.d. within each
input. For clustered or repeated-measures data the test is
anticonservative.

## References

Anderson, M. J. (2001). A new method for non-parametric multivariate
analysis of variance. *Austral Ecology*, 26(1), 32-46.

van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski, J. J., Tio,
P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2022). Comparing
network structures on three aspects: A permutation test. *Psychological
Methods*. [doi:10.1037/met0000476](https://doi.org/10.1037/met0000476)

## See also

[`test_modality_edge_differences`](https://atinakosta.github.io/catgraph/reference/test_modality_edge_differences.md)
for post-hoc edge-wise testing;
[`compare_modality_graphs`](https://atinakosta.github.io/catgraph/reference/compare_modality_graphs.md)
for visual comparison.

## Examples

``` r
# Split survey_health by sex and test whether the joint structure
# differs. Using a small n_perm for the example; in practice use 2000+.
data(survey_health)
df_f <- subset(survey_health, sex == "female")[, -1]
df_m <- subset(survey_health, sex == "male")[, -1]

mg_f <- build_modality_graph(df_f)
mg_m <- build_modality_graph(df_m)

test_result <- test_modality_graph_equality(
  mg_f, mg_m, n_perm = 200, seed = 1, verbose = FALSE
)
print(test_result)
#> Permutation test of modality-graph equality
#>   Statistic       : frobenius 
#>   Test type       : unfiltered 
#>   Stratified      : FALSE 
#>   Sample sizes    : n_x = 206 , n_y = 195 
#>   Permutations    : 200 valid / 200 
#>   Observed stat   : 1.912 
#>   Null mean / sd  : 1.611 / 0.267 
#>   Empirical p     : 0.1493 
#> 
#>   FAIL TO REJECT null at alpha = 0.05 
```
