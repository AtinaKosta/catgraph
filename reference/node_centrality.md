# Compute weighted centrality indices for all variables in a catgraph

Returns a data frame of five weighted centrality measures for every node
(variable) in a `catgraph` object. All measures use the phi / Cramer's V
edge weights stored in the graph, so they reflect association strength
rather than mere connectivity.

## Usage

``` r
node_centrality(x, ...)
```

## Arguments

- x:

  A `catgraph` or `catmodgraph` object.

- ...:

  Additional arguments passed to the method.

## Value

A `data.frame` with one row per variable and the following columns:

- `variable`:

  Variable name (character).

- `strength`:

  Weighted degree / strength centrality: the sum of all edge weights
  incident to a node. The primary measure of how strongly a variable is
  globally associated with all others (Barrat et al., 2004).

- `w_betweenness`:

  Weighted betweenness centrality: the number of weighted shortest paths
  passing through a node, where path lengths are \\1/w\_{ij}\\ so that
  stronger associations correspond to shorter distances (Brandes, 2001).

- `w_closeness`:

  Weighted closeness centrality: the inverse of the mean weighted
  distance to all other reachable nodes.

- `w_eigenvector`:

  Weighted eigenvector centrality: a node scores highly if it is
  connected to other high-scoring nodes (Bonacich, 1987).

- `w_pagerank`:

  Weighted PageRank: a random-walk based centrality robust to isolated
  subgraphs (Brin & Page, 1998). Uses the igraph default damping factor
  of 0.85.

The data frame is returned sorted by `strength` in descending order.

## Details

Edges with `NA` or non-positive weights are dropped from a local copy of
the graph before centrality is computed; the input `catgraph` object is
never modified. If the graph has no edges after this cleanup, a
zero-valued result is returned (with a uniform `1/p` PageRank, the
stationary distribution on an edgeless graph).

**Strength centrality (Barrat et al., 2004):**

\$\$s_i = \sum\_{j \in \mathcal{N}(i)} w\_{ij}\$\$

where \\\mathcal{N}(i)\\ is the set of neighbours of node \\i\\ and
\\w\_{ij}\\ is the edge weight (phi or Cramer's V). In a categorical
association graph this is the most interpretable measure: it is simply
the total association strength of a variable with all other variables.

**Weighted betweenness (Brandes, 2001):**

Path lengths are defined as \\1/w\_{ij}\\ so that stronger edges are
traversed preferentially. A variable with high weighted betweenness sits
on many shortest paths and therefore acts as a statistical bridge
between otherwise weakly associated variable groups.

**Weighted closeness:**

\$\$C_i = \frac{n-1}{\sum\_{j \neq i} d(i,j)}\$\$

where \\d(i,j)\\ is the weighted shortest-path distance. A high value
indicates a variable that reaches all other variables through short
(strong-association) paths.

**Eigenvector centrality (Bonacich, 1987):**

The leading eigenvector of the weighted adjacency matrix. A variable
scores highly when it is associated with other variables that are
themselves highly associated variables — a second-order hub measure.

**PageRank (Brin & Page, 1998):**

A random-walk measure that down-weights the contribution of nodes with
many weak connections. Robust when the graph is sparse or has near-zero
weight edges.

## References

Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A.
(2004). The architecture of complex weighted networks. *Proceedings of
the National Academy of Sciences*, 101(11), 3747–3752.
[doi:10.1073/pnas.0400087101](https://doi.org/10.1073/pnas.0400087101)

Bonacich, P. (1987). Power and centrality: A family of measures.
*American Journal of Sociology*, 92(5), 1170–1182.
[doi:10.1086/228631](https://doi.org/10.1086/228631)

Brandes, U. (2001). A faster algorithm for betweenness centrality.
*Journal of Mathematical Sociology*, 25(2), 163–177.
[doi:10.1080/0022250X.2001.9990249](https://doi.org/10.1080/0022250X.2001.9990249)

Brin, S., & Page, L. (1998). The anatomy of a large-scale hypertextual
web search engine. *Computer Networks*, 30(1-7), 107–117.
[doi:10.1016/S0169-7552(98)00110-X](https://doi.org/10.1016/S0169-7552%2898%2900110-X)

## See also

[`plot_centrality`](https://atinakosta.github.io/catgraph/reference/plot_centrality.md),
[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`clustering_coef`](https://atinakosta.github.io/catgraph/reference/clustering_coef.md)

## Examples

``` r
df <- expand_table(Titanic)
cg <- catgraph(df)
nc <- node_centrality(cg)
nc
#>   variable  strength w_betweenness w_closeness w_eigenvector w_pagerank
#> 1      Sex 1.0000000             0   0.8869209     1.0000000  1.0000000
#> 2    Class 0.9579945             1   1.0000000     0.9284092  0.9719876
#> 3 Survived 0.8777217             0   0.7680133     0.9201859  0.8882785
#> 4      Age 0.4563527             0   0.5423700     0.4929748  0.5186247

# Identify the most central variable
nc$variable[1]
#> [1] "Sex"

# Raw (unnormalised) values
node_centrality(cg, normalize = FALSE)
#>   variable  strength w_betweenness w_closeness w_eigenvector w_pagerank
#> 1      Sex 0.9653402             0  0.08678885     1.0000000  0.2959551
#> 2    Class 0.9247906             2  0.09785410     0.9284092  0.2876647
#> 3 Survived 0.8473000             0  0.07515326     0.9201859  0.2628906
#> 4      Age 0.4405356             0  0.05307313     0.4929748  0.1534896
```
