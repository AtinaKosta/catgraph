# Detect communities of co-associated modalities

Applies graph community detection to a `catmodgraph` object and writes
the resulting community membership onto the graph vertices. Communities
here are groups of modalities (factor levels) that tend to co-associate
across different variables. This is a category co-association analysis,
not a respondent-segmentation method; for the latter, see the poLCA or
FactoMineR packages.

## Usage

``` r
cluster_modalities(x, method = c("louvain", "walktrap"), signed = TRUE)
```

## Arguments

- x:

  A `catmodgraph` object.

- method:

  Character. Community detection method. One of `"louvain"` (default) or
  `"walktrap"`.

- signed:

  Logical. If `TRUE`, only edges with positive standardised Pearson
  residual (attraction: modalities co-occurring more than expected under
  independence) are used for community detection, so communities are
  defined solely by positive co-association. Repulsion edges (negative
  `std_resid`) are excluded from the clustering graph but are retained
  on the original graph for downstream plotting with
  `plot(mg, signed = TRUE)`. Only supported for `method = "louvain"`.
  Default `TRUE` (changed from v0.9.0 which defaulted to `FALSE`).

## Value

The input `catmodgraph` object with additional vertex attribute
`cluster`, and a new component `membership` giving the cluster label for
each modality node.

## Details

**Why `signed = TRUE` is now the default.** The `weight` edge attribute
stores the absolute phi coefficient so that edge thickness in plots
scales with association strength regardless of direction. When unsigned
Louvain uses these absolute weights, pairs of modalities with large
*negative* residuals (strong repulsion, e.g. `smoking_status=current`
and `lung_disease=no`) are pulled into the same community — the opposite
of what subject-matter knowledge expects. Restricting clustering to
positive-residual edges produces communities defined by genuine
co-occurrence surplus, which are substantively more interpretable.

The `std_resid` edge attribute is retained on all edges in the original
graph, so `plot(mg, signed = TRUE)` still shows both attraction (green)
and repulsion (red) edges even when communities were detected on the
positive-only subgraph.

**Signed Louvain is a pragmatic adaptation.** Modern igraph rejects
negative edge weights in `cluster_louvain()`, so `signed = TRUE` simply
drops repulsion edges rather than implementing a proper signed
algorithm. For principled signed community detection (Traag & Bruggeman,
2009), see the signnet package.

`signed = TRUE` requires the `std_resid` edge attribute, which is
attached by
[`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
and preserved through
[`prune_modality_edges`](https://atinakosta.github.io/catgraph/reference/prune_modality_edges.md).

## References

Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
Fast unfolding of communities in large networks. *Journal of Statistical
Mechanics*, 2008(10), P10008.
[doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)

Traag, V. A., & Bruggeman, J. (2009). Community detection in networks
with positive and negative links. *Physical Review E*, 80(3), 036115.
[doi:10.1103/PhysRevE.80.036115](https://doi.org/10.1103/PhysRevE.80.036115)

## Examples

``` r
df <- expand_table(Titanic)
mg <- build_modality_graph(df)

# Default: positive-only Louvain (recommended)
mg <- cluster_modalities(mg)
table(mg$membership)
#> 
#> 1 2 
#> 6 4 

# Legacy unsigned behaviour: abs(phi) drives community detection
mg_abs <- cluster_modalities(mg, signed = FALSE)
table(mg_abs$membership)
#> 
#> 1 2 
#> 7 3 
```
