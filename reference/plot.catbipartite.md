# Plot a bipartite respondent-modality graph

Visualises a `catbipartite` object with respondents on one side and
modalities on the other. Respondent nodes are plotted small and
unlabelled; modality nodes are plotted large and coloured by originating
variable.

## Usage

``` r
# S3 method for class 'catbipartite'
plot(
  x,
  show_respondents = TRUE,
  max_respondents = 500L,
  vertex_size_mod = 18,
  vertex_size_row = 2,
  ...
)
```

## Arguments

- x:

  A `catbipartite` object.

- show_respondents:

  Logical. If `TRUE` (default), respondent-side vertices are drawn as
  small dots. If `FALSE`, only the modality partition is plotted (useful
  when the number of respondents would make the plot unreadable).

- max_respondents:

  Integer or `NULL`. If the number of respondents exceeds this value and
  `show_respondents = TRUE`, a random sample of `max_respondents` is
  drawn. Set to `NULL` to disable. Default `500`.

- vertex_size_mod:

  Numeric. Vertex size for modality nodes. Default `18`.

- vertex_size_row:

  Numeric. Vertex size for respondent nodes. Default `2`.

- ...:

  Further arguments passed to `plot.igraph()`.

## Value

Invisibly returns the input object.

## Examples

``` r
df <- expand_table(Titanic)
bg <- bipartite_modality_graph(df)
plot(bg)

plot(bg, show_respondents = FALSE)

```
