# Plot gravity indices alongside traditional centrality for a catmodgraph

Produces a 2 x 3 panel figure showing six structural measures for the
same modality network on a shared layout: strength, betweenness, and
eigenvector centrality (top row) alongside MGI+, OS, and dMGI (bottom
row). Node size encodes the magnitude of each measure; colour encodes
community membership (top row) or attractor/satellite role (dMGI panel).

## Usage

``` r
plot_gravity(
  x,
  gravity = NULL,
  layout_fn = igraph::layout_with_fr,
  community_attr = "cluster",
  palette = NULL,
  seed = 1L,
  title = "Modality network: traditional centrality vs. gravity indices",
  node_size_range = c(4, 22),
  show_labels = FALSE,
  bars = FALSE,
  bars_n = 12L,
  attractor_col = "#1D9E75",
  satellite_col = "#D85A30",
  ...
)
```

## Arguments

- x:

  A `catmodgraph` object.

- gravity:

  A data frame returned by
  [`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md).
  If `NULL` (default), `modality_gravity(x)` is called internally.

- layout_fn:

  An igraph layout function. Defaults to
  [`igraph::layout_with_fr`](https://r.igraph.org/reference/layout_with_fr.html).

- community_attr:

  Character. Vertex attribute name for community membership (set by
  [`cluster_modalities`](https://atinakosta.github.io/catgraph/reference/cluster_modalities.md)).
  Defaults to `"cluster"`.

- palette:

  Character vector of colours for communities. If `NULL`,
  [`grDevices::hcl.colors`](https://rdrr.io/r/grDevices/palettes.html)
  with palette `"Dark 3"` is used.

- seed:

  Integer. Random seed for reproducible layouts. Default `1L`.

- title:

  Character. Overall figure title. Default
  `"Modality network: traditional centrality vs. gravity indices"`.

- node_size_range:

  Numeric vector of length 2. Minimum and maximum node sizes. Default
  `c(4, 22)`.

- show_labels:

  Logical. Whether to draw node labels. Default `FALSE`.

- bars:

  Logical. If `TRUE`, a second figure is produced showing a 2 x 3 bar
  chart panel with the top `bars_n` modalities per measure, coloured by
  community. Default `FALSE`.

- bars_n:

  Integer. Number of top modalities to show per bar chart panel when
  `bars = TRUE`. Default `12L`.

- attractor_col:

  Colour for attractor nodes in the dMGI panel. Default `"#1D9E75"`.

- satellite_col:

  Colour for satellite nodes in the dMGI panel. Default `"#D85A30"`.

- ...:

  Additional arguments passed to
  [`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md)
  when `gravity = NULL`.

## Value

Invisibly returns the `gravity` data frame used for plotting. The
primary effect is the figure drawn on the current graphics device.

## See also

[`modality_gravity`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md),
[`compare_gravity`](https://atinakosta.github.io/catgraph/reference/compare_gravity.md)

## Examples

``` r
data(survey_health)
mg  <- build_modality_graph(survey_health)
mg  <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
mg  <- cluster_modalities(mg, method = "louvain")
plot_gravity(mg)

```
