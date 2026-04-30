# Detect variable communities in a catgraph using graph clustering algorithms

Applies a community detection algorithm to the weighted undirected graph
of a `catgraph` object. Assigns cluster membership as a vertex attribute
and returns the updated object. The Louvain algorithm (Blondel et al.,
2008) is the default as it is fast, handles weighted edges, and requires
no pre-specification of the number of clusters.

## Usage

``` r
detect_clusters(
  x,
  method = c("louvain", "walktrap", "fast_greedy", "label_prop", "edge_betweenness"),
  resolution = 1,
  steps = 4L
)
```

## Arguments

- x:

  A `catgraph` object.

- method:

  Character. Community detection algorithm to use. One of `"louvain"`
  (default), `"walktrap"`, `"fast_greedy"`, `"label_prop"`, or
  `"edge_betweenness"`. All methods use edge weights when available.

- resolution:

  Numeric. Resolution parameter for the Louvain method
  ([`igraph::cluster_louvain`](https://r.igraph.org/reference/cluster_louvain.html)
  `resolution` argument). Higher values favour smaller communities.
  Default `1`.

- steps:

  Integer. Number of random walk steps for the Walktrap method. Default
  `4L`.

## Value

The input `catgraph` object with two additions:

- `graph`:

  Vertex attribute `cluster` (integer membership vector) and
  `cluster_method` (character) are added to the graph.

- `clustering`:

  A named list with the raw igraph communities object (`communities`),
  the membership vector (`membership`), the number of communities
  (`n_clusters`), and the modularity score (`modularity`).

## Details

Cluster detection is performed on the graph at the time of the call. If
you want to cluster a pruned graph, call
[`prune_edges`](https://atinakosta.github.io/catgraph/reference/prune_edges.md)
first. This function clusters vertices of the *variable-level*
association graph. Since vertices represent variables rather than factor
levels or respondents, the detected communities should be read as blocks
of variables that co-vary pairwise. For community detection at the
modality level (factor-level co-association), see
[`cluster_modalities`](https://atinakosta.github.io/catgraph/reference/cluster_modalities.md).

The Louvain algorithm is non-deterministic; set a seed with
[`set.seed()`](https://rdrr.io/r/base/Random.html) before calling if
reproducibility is needed.

## References

Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
Fast unfolding of communities in large networks. *Journal of Statistical
Mechanics: Theory and Experiment*, 2008(10), P10008.
[doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)

Pons, P., & Latapy, M. (2005). Computing communities in large networks
using random walks. In *Computer and Information Sciences – ISCIS 2005*
(pp. 284–293). Springer.
[doi:10.1007/11569596_31](https://doi.org/10.1007/11569596_31)

## See also

[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`prune_edges`](https://atinakosta.github.io/catgraph/reference/prune_edges.md)

## Examples

``` r
df <- as.data.frame(Titanic)
df_exp <- df[rep(seq_len(nrow(df)), df$Freq), -5]
cg <- catgraph(df_exp)
set.seed(42)
cg <- detect_clusters(cg)
cg$clustering$n_clusters
#> [1] 1
cg$clustering$modularity
#> [1] 0
```
