# Construct a bipartite respondent-modality graph

Builds a two-mode (bipartite) graph where one partition of vertices
represents respondents (rows of the data) and the other represents
modalities (factor levels). An edge connects respondent \\i\\ to
modality \\V{=}l\\ if row \\i\\ has value \\l\\ on variable \\V\\.

## Usage

``` r
bipartite_modality_graph(data, remove_na = TRUE, row_prefix = "r")
```

## Arguments

- data:

  A data frame of categorical variables.

- remove_na:

  Logical. If `TRUE` (default), rows with any missing value are dropped
  before graph construction.

- row_prefix:

  Character. Prefix used to name respondent-side vertices. Default
  `"r"`, so respondents are labelled `"r1", "r2", ...`.

## Value

An object of class `"catbipartite"`, a list with:

- `graph`:

  An `igraph` undirected graph with vertex attribute `type` (FALSE =
  respondent, TRUE = modality).

- `modalities`:

  Data frame describing modality-side vertices, with columns `node`,
  `variable`, `modality`.

- `n_rows`:

  Number of respondent-side vertices.

- `n_modalities`:

  Number of modality-side vertices.

- `data`:

  The processed categorical data used to build the graph (after NA
  removal).

## Details

Bipartite (two-mode) graphs preserve the full respondent-to-modality
incidence structure. They are the unprojected counterpart of the
modality co-association graph returned by
[`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md):
projecting a `catbipartite` onto its modality partition yields a
respondents-in-common modality network, which is a **count-based
co-occurrence network** distinct from the chi-square / phi-filtered
network produced by `build_modality_graph`.

The bipartite graph is the correct tool when the scientific object of
interest is the raw incidence relationship between units and categories
— for example, affiliation networks (persons to events), ecological
surveys (sites to species), or survey data seen as respondents endorsing
response levels.

The resulting graph satisfies the igraph convention for bipartite
graphs: `igraph::bipartite_projection(g)` can be called directly to
obtain either the respondent-side or the modality-side projection.

**Scope.** A `catbipartite` is a raw incidence graph; it applies no
statistical filter. Edges reflect endorsement, not association. For
statistically-filtered category co-association, see
[`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md).

## See also

[`build_modality_graph`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
for the chi-square / phi-filtered modality co-association graph.

## Examples

``` r
df <- expand_table(Titanic)
bg <- bipartite_modality_graph(df)
bg
#> catbipartite object (respondent-modality incidence graph)
#>   Respondent nodes : 2201 
#>   Modality nodes   : 10 
#>   Edges            : 8804 
#>   Variables        : Class, Sex, Age, Survived 
summary(bg)
#> $n_rows
#> [1] 2201
#> 
#> $n_modalities
#> [1] 10
#> 
#> $n_edges
#> [1] 8804
#> 
#> $variables
#> [1] "Class"    "Sex"      "Age"      "Survived"
#> 
#> $modalities_per_variable
#> 
#>      Age    Class      Sex Survived 
#>        2        4        2        2 
#> 

# Project onto modalities: edges weighted by shared respondents
proj <- igraph::bipartite_projection(bg$graph, which = "true")
```
