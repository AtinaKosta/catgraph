# Summarise modality communities in a catmodgraph

Summarises the structure of communities detected in a `catmodgraph`
object. The output reports which modalities belong to each community,
which variables contribute to each community, and the internal edge
cohesion. This is descriptive output for interpreting the community
structure of a modality co-association graph.

## Usage

``` r
summarise_modality_communities(x)
```

## Arguments

- x:

  A `catmodgraph` object that has already been processed by
  [`cluster_modalities()`](https://atinakosta.github.io/catgraph/reference/cluster_modalities.md).

## Value

A list of class `"catmodcommunity"` with components:

- `community_summary`:

  Data frame with one row per modality community, reporting size,
  variables represented, and mean internal edge weight.

- `community_members`:

  Data frame listing all modality nodes and their assigned community.

- `variable_composition`:

  Table of variable counts by community, useful for seeing which
  variables span which communities.

## Details

Communities are interpreted as **groups of modalities (factor levels)
that co-associate across different variables**, not as respondent
segments or latent classes. For respondent-level segmentation, use the
poLCA or FactoMineR packages.

## Examples

``` r
df <- expand_table(Titanic)
mg <- build_modality_graph(df)
mg <- cluster_modalities(mg)
comm <- summarise_modality_communities(mg)
comm
#> catmodcommunity object
#>   Modality communities: 2 
#>   Community sizes     : 6, 4 
```
