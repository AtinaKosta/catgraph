# Prune edges from a catgraph by effect size or adjusted p-value

Removes edges whose effect size or (adjusted) p-value does not meet a
specified threshold, returning a new `catgraph` object with a sparser
graph. Multiple-testing adjustment is applied across all edges by
default.

## Usage

``` r
prune_edges(
  x,
  min_weight = 0,
  max_p = 1,
  p_adjust = c("BH", "holm", "bonferroni", "none"),
  remove_isolates = FALSE
)
```

## Arguments

- x:

  A `catgraph` object.

- min_weight:

  Numeric, non-negative. Edges with effect size strictly below this
  value are removed. Since phi and Cramer's V lie in \[0, 1\],
  meaningful thresholds are in that range; values \>= 1 remove all
  edges. Effect-size pruning is the **primary** filter and is always a
  safer choice than p-value pruning, because the package is centred on
  effect sizes. Default `0` (no filtering).

- max_p:

  Numeric in \[0, 1\]. Edges with *adjusted* p-value strictly above this
  value are removed. Adjustment method is controlled by `p_adjust`.
  Default `1` (no filtering).

- p_adjust:

  Character. Multiple-testing correction applied across all edges
  (`choose(p, 2)` simultaneous tests). One of:

  `"BH"`

  :   Benjamini-Hochberg false discovery rate (default). Recommended for
      exploratory work.

  `"holm"`

  :   Holm-Bonferroni step-down; strong family-wise error rate control.

  `"bonferroni"`

  :   Bonferroni; conservative FWER control.

  `"none"`

  :   Raw p-values (unadjusted). Not recommended when many variables are
      analysed; retained for reproducing pre-0.4.0 behaviour.

- remove_isolates:

  Logical. If `TRUE`, vertices with degree 0 after pruning are also
  removed. Default `FALSE`.

## Value

A `catgraph` object with the filtered graph. The graph gains two new
edge attributes: `p_value_adj` (adjusted p-values) and `p_adjust_method`
(the method string).

## Details

Pruning uses the edge attribute `weight` (the active effect size) and
the *adjusted* `p_value`. Both thresholds apply simultaneously; an edge
is retained only when **both** conditions are met.

Edges with `NA` weights or p-values (from degenerate pairs) are always
removed.

**Chained calls and multiplicity scoping.** Multiple-testing correction
is applied across the edges present in the graph at the time of the
call. When `prune_edges()` is called on a graph that has *already* been
pruned with a non-`"none"` p-value adjustment, the second call
re-adjusts on the surviving subset, not on the original `choose(p, 2)`
tests. This is anti-conservative and the function emits a warning in
this case. To change the adjustment method mid-analysis, rebuild the
`catgraph` with
[`catgraph()`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
and prune once; do not chain two adjusted prunes. A single
`prune_edges()` call that specifies both `min_weight` and `max_p` is
always safe because the BH / Holm denominators are computed before any
filtering.

Conventional Cohen (1988) thresholds for phi and Cramer's V: small
\\\approx 0.1\\, medium \\\approx 0.3\\, large \\\geq 0.5\\.

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. *JRSS-B*,
57(1), 289–300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics*, 6(2), 65–70.

Cohen, J. (1988). *Statistical Power Analysis for the Behavioral
Sciences* (2nd ed.). Lawrence Erlbaum Associates.

## See also

[`catgraph`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`detect_clusters`](https://atinakosta.github.io/catgraph/reference/detect_clusters.md)

## Examples

``` r
df <- expand_table(Titanic)
cg <- catgraph(df)

# Default: BH-adjusted p-values, effect-size floor of 0.1
cg_pruned <- prune_edges(cg, min_weight = 0.1, max_p = 0.05)
cg_pruned
#> catgraph object (pairwise association network)
#>   Variables : 4 
#>   Edges     : 5 
#>   Method    : Cramer's V (classical) 
#>   Weights   : min = 0.1110  median = 0.2941  max = 0.4556
#>   Note      : edges encode pairwise marginal association, not
#>               conditional independence. All metrics lie on [0, 1].
#>               NMI / AMI weights are not exchangeable with Cramer's V
#>               weights across graph objects. See vignette
#>               'Methodological caveats'.

# Stricter: Holm adjustment
prune_edges(cg, min_weight = 0.1, max_p = 0.05, p_adjust = "holm")
#> catgraph object (pairwise association network)
#>   Variables : 4 
#>   Edges     : 5 
#>   Method    : Cramer's V (classical) 
#>   Weights   : min = 0.1110  median = 0.2941  max = 0.4556
#>   Note      : edges encode pairwise marginal association, not
#>               conditional independence. All metrics lie on [0, 1].
#>               NMI / AMI weights are not exchangeable with Cramer's V
#>               weights across graph objects. See vignette
#>               'Methodological caveats'.

# Pre-0.4.0 behaviour (raw p-values)
prune_edges(cg, min_weight = 0.1, max_p = 0.05, p_adjust = "none")
#> catgraph object (pairwise association network)
#>   Variables : 4 
#>   Edges     : 5 
#>   Method    : Cramer's V (classical) 
#>   Weights   : min = 0.1110  median = 0.2941  max = 0.4556
#>   Note      : edges encode pairwise marginal association, not
#>               conditional independence. All metrics lie on [0, 1].
#>               NMI / AMI weights are not exchangeable with Cramer's V
#>               weights across graph objects. See vignette
#>               'Methodological caveats'.
```
