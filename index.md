# catgraph

`catgraph` provides **network-based exploratory analysis of categorical
data** at two complementary levels.

1.  **Variable-level association network.** Variables are nodes. Edges
    are weighted by the phi coefficient (2x2 tables) or Cramer’s V
    (larger tables), with optional bias correction (Bergsma, 2013). The
    workflow supports structural exploration, edge pruning with
    multiple-testing adjustment, bootstrap confidence intervals, and
    descriptive network summaries such as centrality and community
    structure.

2.  **Modality-level co-association network.** Modalities (factor
    levels) are nodes. Cross-variable edges are weighted by absolute phi
    coefficients, with signed standardised Pearson residuals stored
    separately to indicate whether co-occurrence is above or below
    independence expectation. The workflow supports edge pruning, signed
    edge visualisation, and community detection over modalities.

    The modality layer sits in the tradition of Multiple Correspondence
    Analysis and two-mode affiliation networks. It operates on pairwise
    associations and does not model higher-order interactions or
    conditional dependencies. It is a descriptive **category
    co-association map**, not a respondent-segmentation tool. For
    respondent segmentation use
    [`poLCA`](https://cran.r-project.org/package=poLCA) or
    [`FactoMineR::HCPC`](https://cran.r-project.org/package=FactoMineR).

3.  **Modality gravity indices.** A novel extension to standard graph
    centrality that incorporates the empirical prevalence of each
    modality. The **Modality Gravity Index (MGI)** and **Orbital Score
    (OS)** identify which modalities act as gravitational attractors
    (dominant, pulling rarer modalities toward them) and which are
    satellites (minority modalities orbiting more prevalent ones). This
    addresses a fundamental limitation of standard centrality indices,
    which treat all nodes as exchangeable regardless of their empirical
    frequency.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("AtinaKosta/catgraph")
```

## What’s new in 0.10.0

- **[`cluster_modalities()`](https://atinakosta.github.io/catgraph/reference/cluster_modalities.md)
  now defaults to `signed = TRUE`** (breaking change from 0.9.0).
  Communities are defined by positive co-association only — edges where
  modalities co-occur *less* than expected under independence (negative
  standardised Pearson residual) are excluded from clustering. This
  produces substantively more interpretable communities: for example,
  `smoking_status=current` and `lung_disease=no` are no longer pulled
  into the same community by their large absolute phi weight despite
  being a repulsion pair. Use `signed = FALSE` to restore the previous
  behaviour.

- **[`build_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
  now stores `phi_signed`** as an additional edge attribute alongside
  `weight` (absolute phi), `p_value`, and `std_resid`. All downstream
  functions are unaffected.

## Quick start

``` r

library(catgraph)
data(survey_health)

# Variable-level: which categorical variables show pairwise association?
cg <- catgraph(survey_health, corrected = TRUE)
cg_p <- prune_edges(cg, min_weight = 0.05, max_p = 0.05, p_adjust = "BH")
plot(cg_p)

# Modality-level: which category levels tend to co-occur across variables?
mg <- build_modality_graph(survey_health)
mg <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
mg <- cluster_modalities(mg)
plot(mg, color_by = "cluster", signed = TRUE)

# Gravity indices: which modalities are attractors vs satellites?
grav <- modality_gravity(mg)
print(grav)                       # role-grouped: ATTRACTORS / SATELLITES
summary(grav)                     # role counts, Spearman rho diagnostic
plot_gravity(mg)                  # 6-panel: traditional centrality vs MGI
plot_gravity_scatter(grav, mg)    # eigenvector vs dMGI contradiction plot

# Compare gravity profiles across subgroups
mg_f <- build_conditional_modality_graph(survey_health, given = list(sex = "female"))
mg_m <- build_conditional_modality_graph(survey_health, given = list(sex = "male"))
mg_f <- prune_modality_edges(mg_f, min_weight = 0.10, max_p = 0.05)
mg_m <- prune_modality_edges(mg_m, min_weight = 0.10, max_p = 0.05)
compare_gravity(list(female = mg_f, male = mg_m))

# Formal test: do two groups differ in overall association structure?
test_modality_graph_equality(mg_f, mg_m, n_perm = 500)
```

See
[`vignette("introduction", package = "catgraph")`](https://atinakosta.github.io/catgraph/articles/introduction.md)
for the full variable-level and modality-level workflow, and
[`vignette("comparison", package = "catgraph")`](https://atinakosta.github.io/catgraph/articles/comparison.md)
for a worked comparison of `catgraph`’s modality layer with MCA,
bipartite affiliation networks, and naive co-occurrence projection.

## Scope — use `catgraph` for

| Question | Tool |
|----|----|
| Which categorical variables co-vary pairwise? | [`catgraph()`](https://atinakosta.github.io/catgraph/reference/catgraph.md) |
| Which category levels bundle together across variables? | [`build_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md) |
| Which modalities are structurally dominant vs peripheral? | [`modality_gravity()`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md) |
| How do category-level association patterns differ across groups? | `compare_*_graphs()` |
| What’s the unprojected respondent-modality incidence like? | [`bipartite_modality_graph()`](https://atinakosta.github.io/catgraph/reference/bipartite_modality_graph.md) |

## Scope — do not use `catgraph` for

- Causal inference or claims about direct effects.
- Conditional-independence structure — use a graphical model instead
  (`bnlearn`, `gRim`).
- Respondent segmentation or latent-class analysis — use `poLCA` or
  [`FactoMineR::HCPC()`](https://rdrr.io/pkg/FactoMineR/man/HCPC.html).
- Estimating conditional relationships between variables — edges are
  marginal and do not control for other variables.

## Interpretation

All graphs produced by `catgraph` represent **marginal association
structure**. Edge weights quantify pairwise dependence and should be
interpreted descriptively. Differences between graphs indicate changes
in association patterns, not necessarily changes in underlying causal
mechanisms.

Gravity indices (MGI, OS) incorporate empirical modality prevalence. A
positive dMGI indicates a modality that exerts net gravitational pull
over less prevalent neighbours; a negative dMGI indicates a satellite
modality being pulled toward more prevalent ones.

## License

GPL (\>= 3)
