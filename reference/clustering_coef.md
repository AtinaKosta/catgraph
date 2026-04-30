# Compute weighted clustering coefficients for all variables in a catgraph

Returns a data frame of five local clustering coefficients for every
node in a `catgraph` object. Four are established weighted extensions of
the Watts-Strogatz coefficient; the fifth (redundancy) is a measure
specific to effect-size graphs that quantifies how much of a pairwise
association is explained by indirect paths through other variables.

## Usage

``` r
clustering_coef(x, method = "all", normalize = TRUE)
```

## Arguments

- x:

  A `catgraph` object.

- method:

  Character vector. Which coefficient(s) to compute. Any subset of
  `c("watts_strogatz", "barrat", "onnela", "zhang", "redundancy")` or
  `"all"` (default).

- normalize:

  Logical. If `TRUE` (default), each coefficient is normalised to \[0,
  1\] by dividing by its theoretical maximum. For the four classical
  coefficients the theoretical maximum is already 1. For redundancy the
  observed maximum is used.

## Value

A `data.frame` with one row per variable. Columns are `variable` plus
one column per requested method. The data frame is sorted by the mean
clustering coefficient across all requested methods in descending order.

## Details

All five coefficients measure the extent to which the neighbours of a
node are also connected to each other — the triangle-closing tendency —
but they differ in how they incorporate edge weights.

**Watts-Strogatz (unweighted; Watts & Strogatz, 1998):**

\$\$C_i^{WS} = \frac{t_i}{k_i(k_i-1)}\$\$

where \\t_i\\ is the number of triangles through node \\i\\ and \\k_i\\
is its degree. This is the classical unweighted coefficient included as
a baseline.

**Barrat et al. (2004):**

\$\$C_i^{B} = \frac{1}{s_i(k_i-1)} \sum\_{j,h} \frac{w\_{ij}+w\_{ih}}{2}
a\_{ij} a\_{ih} a\_{jh}\$\$

where \\s_i = \sum_j w\_{ij}\\ is the node strength and \\a\_{ij}\\ is
the binary adjacency. Weights the contribution of each triangle by the
average weight of the two edges incident to the focal node. Nodes that
close triangles through their strongest edges score highly.

**Onnela et al. (2005):**

\$\$C_i^{O} = \frac{1}{k_i(k_i-1)} \sum\_{j,h} (\hat{w}\_{ij}
\hat{w}\_{ih} \hat{w}\_{jh})^{1/3}\$\$

where \\\hat{w} = w / \max(w)\\ are normalised weights. Uses the
geometric mean of triangle edge weights, making it sensitive to weak
ties: a triangle with one very weak edge scores low.

**Zhang & Horvath (2005):**

\$\$C_i^{ZH} = \frac{\sum\_{j,h} w\_{ij} w\_{ih} w\_{jh}}
{\left(\sum\_{j \neq i} w\_{ij}\right)^2 - \sum\_{j \neq i}
w\_{ij}^2}\$\$

The numerator sums products of three edge weights around all triangles;
the denominator is the squared strength minus the sum of squared
weights. Originally proposed for co-expression networks (WGCNA).
Amplifies strong-edge triangles via cubic weighting.

**Redundancy (original):**

For each edge \\(i,j)\\ in the graph, the redundancy ratio is defined as
the direct weight divided by the maximum indirect path weight:

\$\$R\_{ij} = \frac{w\_{ij}}{\max\_{k \neq i,j} \min(w\_{ik},
w\_{kj})}\$\$

The node-level redundancy coefficient is the mean of \\R\_{ij}\\ over
all edges incident to node \\i\\. A value near 1 indicates that the
direct association between two variables is no stronger than what their
shared neighbours would predict — the edge may be mediated. A value well
above 1 indicates a genuine direct association that exceeds indirect
explanation. When no indirect path exists (k \< 3), the edge is assigned
\\R\_{ij} = \infty\\ and excluded from the node mean.

## References

Barrat, A., Barthelemy, M., Pastor-Satorras, R., & Vespignani, A.
(2004). The architecture of complex weighted networks. *Proceedings of
the National Academy of Sciences*, 101(11), 3747–3752.
[doi:10.1073/pnas.0400087101](https://doi.org/10.1073/pnas.0400087101)

Onnela, J.-P., Saramaki, J., Kertesz, J., & Kaski, K. (2005). Intensity
and coherence of motifs in weighted complex networks. *Physical Review
E*, 71(6), 065103.
[doi:10.1103/PhysRevE.71.065103](https://doi.org/10.1103/PhysRevE.71.065103)

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
'small-world' networks. *Nature*, 393(6684), 440–442.
[doi:10.1038/30918](https://doi.org/10.1038/30918)

Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
co-expression network analysis. *Statistical Applications in Genetics
and Molecular Biology*, 4(1), Article 17.
[doi:10.2202/1544-6115.1128](https://doi.org/10.2202/1544-6115.1128)

## See also

[`compare_clustering`](https://atinakosta.github.io/catgraph/reference/compare_clustering.md),
[`plot_clustering`](https://atinakosta.github.io/catgraph/reference/plot_clustering.md),
[`node_centrality`](https://atinakosta.github.io/catgraph/reference/node_centrality.md)

## Examples

``` r
df <- expand_table(Titanic)
cg <- catgraph(df)

cc <- clustering_coef(cg)
cc
#>   variable watts_strogatz barrat    onnela     zhang redundancy
#> 1    Class              1    0.5 0.2861019 0.2530913  1.0000000
#> 2      Sex              1    0.5 0.2795694 0.2482151  0.8088828
#> 3      Age              1    0.5 0.2106029 0.3691475  0.7145339
#> 4 Survived              1    0.5 0.2688767 0.3138760  0.6472822

# Single method
clustering_coef(cg, method = "barrat")
#>   variable barrat
#> 1      Sex    0.5
#> 2      Age    0.5
#> 3    Class    0.5
#> 4 Survived    0.5

# Compare all methods
compare_clustering(cg)
#>   variable watts_strogatz barrat    onnela     zhang redundancy   mean_cc
#> 1    Class              1    0.5 0.2861019 0.2530913  1.0000000 0.6078387
#> 2      Sex              1    0.5 0.2795694 0.2482151  0.8088828 0.5673334
#> 3      Age              1    0.5 0.2106029 0.3691475  0.7145339 0.5588569
#> 4 Survived              1    0.5 0.2688767 0.3138760  0.6472822 0.5460070
```
