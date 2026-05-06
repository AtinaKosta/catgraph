# Changelog

## catgraph 0.11.0

### New feature: multiple association paradigms

This release introduces three additional association metrics for edge
weighting, expanding catgraph from a single-paradigm tool to a
multi-paradigm framework covering frequentist, information-theoretic,
and Bayesian traditions. All four metrics are available at both the
variable level
([`catgraph()`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`build_graph()`](https://atinakosta.github.io/catgraph/reference/build_graph.md),
[`assoc_similarity()`](https://atinakosta.github.io/catgraph/reference/assoc_similarity.md))
and the modality level
([`build_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)).

#### New `method` argument

[`catgraph()`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`build_graph()`](https://atinakosta.github.io/catgraph/reference/build_graph.md),
[`assoc_similarity()`](https://atinakosta.github.io/catgraph/reference/assoc_similarity.md),
and
[`build_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
now accept a `method` argument (default `"cramers_v"`) controlling which
association measure is used as edge weights:

- `"cramers_v"` — classical phi / Cramér’s V (previous default, fully
  backward-compatible).
- `"cramers_v_corrected"` — bias-corrected Cramér’s V via Bergsma
  (2013). Previously accessed through `corrected = TRUE`; that argument
  is retained for backward compatibility and now resolves internally to
  `method = "cramers_v_corrected"`.
- `"nmi"` — Normalised Mutual Information. Symmetric, bounded \[0, 1\],
  information-theoretic. Sensitive to marginal entropy imbalances;
  recommended when variables have unequal numbers of categories or
  skewed marginal distributions.
- `"ami"` — Adjusted Mutual Information. Subtracts expected MI under
  random permutation (Vinh et al., 2010), correcting for the upward bias
  of NMI on sparse contingency tables. Recommended when any variable has
  rare categories (expected cell counts \< 5).
- `"bayesian_cramers_v"` — Dirichlet-smoothed Cramér’s V. Applies a
  symmetric Dirichlet(alpha) prior to cell counts before computing the
  association, stabilising edge weights when tables are sparse. Default
  `alpha = 0.5` (Jeffreys prior). Converges to classical Cramér’s V as n
  → ∞.

#### New `alpha` argument

[`catgraph()`](https://atinakosta.github.io/catgraph/reference/catgraph.md),
[`build_graph()`](https://atinakosta.github.io/catgraph/reference/build_graph.md),
[`assoc_similarity()`](https://atinakosta.github.io/catgraph/reference/assoc_similarity.md),
and
[`build_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
accept `alpha` (default `0.5`), the Dirichlet prior concentration for
`method = "bayesian_cramers_v"`. Ignored for all other methods.
Validated at call time.

#### New exported functions

- `nmi_assoc(x, y, adjusted = FALSE)` — computes NMI or AMI for a pair
  of categorical vectors. Returns the same list structure as
  [`effect_size()`](https://atinakosta.github.io/catgraph/reference/effect_size.md)
  for consistency.
- `bayesian_cramers_v(x, y, alpha = 0.5)` — computes Dirichlet- smoothed
  Cramér’s V for a pair of categorical vectors. p-value is taken from
  the unsmoothed chi-square test to avoid anti-conservative inference.

#### `catgraph` object: new fields

- `$method` — character string recording the association paradigm used.
- `$alpha` — the Dirichlet prior used (`NA_real_` when method is not
  `"bayesian_cramers_v"`).
- `$corrected` — retained for backward compatibility; now derived from
  `$method` rather than stored independently.

#### `catmodgraph` object: new fields

- `$method` — as above.
- `$alpha` — as above.

#### `print()` and `summary()` changes

[`print.catgraph()`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
and
[`summary.catgraph()`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
now display a human-readable `Method` line and, when applicable, an
`Alpha` line showing the Dirichlet prior. The old `Estimator` line and
`Metric mix` line have been replaced.

### Backward compatibility

All existing code using `catgraph(df)`,
`catgraph(df, corrected = TRUE)`, `build_graph(df)`, or
`assoc_similarity(df)` continues to work without any changes. The
`corrected` argument is retained and silently resolved to
`method = "cramers_v_corrected"`.

### References

Bergsma, W. (2013). A bias-correction for Cramér’s V and Tschuprow’s T.
*Journal of the Korean Statistical Society*, 42(3), 323–328.

Cover, T. M., & Thomas, J. A. (2006). *Elements of Information Theory*
(2nd ed.). Wiley.

Good, I. J. (1965). *The Estimation of Probabilities*. MIT Press.

Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic
measures for clusterings comparison. *Journal of Machine Learning
Research*, 11, 2837–2854.

## catgraph 0.10.0

### Breaking change in default behaviour

- [`cluster_modalities()`](https://atinakosta.github.io/catgraph/reference/cluster_modalities.md)
  now defaults to `signed = TRUE` (was `FALSE`). Communities are now
  defined by **positive co-association only** (edges with negative
  standardised Pearson residual — repulsion — are excluded from the
  clustering graph). This produces substantively more interpretable
  communities: modalities that co-occur *less* than expected
  (e.g. `smoking_status=current` and `lung_disease=no`) are no longer
  pulled into the same community by their large absolute phi weight.
  Users who need the previous unsigned behaviour can pass
  `signed = FALSE` explicitly.

### New edge attribute in `build_modality_graph()`

- [`build_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
  now stores `phi_signed` as an additional edge attribute alongside the
  existing `weight` (absolute phi), `p_value`, and `std_resid`. The
  `weight` attribute is unchanged so all downstream functions (gravity
  indices, centrality, pruning, plotting) are unaffected. `phi_signed`
  is available for users who want the directed association value
  directly from the graph object.

### Functions unaffected by this change

The following functions operate on `weight` (absolute phi) or
`std_resid` directly and are unaffected by the clustering default
change:
[`modality_gravity()`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md),
[`node_centrality()`](https://atinakosta.github.io/catgraph/reference/node_centrality.md),
[`plot_gravity()`](https://atinakosta.github.io/catgraph/reference/plot_gravity.md),
[`prune_modality_edges()`](https://atinakosta.github.io/catgraph/reference/prune_modality_edges.md),
[`plot.catmodgraph()`](https://atinakosta.github.io/catgraph/reference/plot.catmodgraph.md),
[`compare_modality_graphs()`](https://atinakosta.github.io/catgraph/reference/compare_modality_graphs.md),
[`test_modality_graph_equality()`](https://atinakosta.github.io/catgraph/reference/test_modality_graph_equality.md),
[`test_modality_edge_differences()`](https://atinakosta.github.io/catgraph/reference/test_modality_edge_differences.md),
[`build_conditional_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_conditional_modality_graph.md),
[`joint_balance()`](https://atinakosta.github.io/catgraph/reference/joint_balance.md).

## catgraph 0.9.0

### New functions

- [`modality_gravity()`](https://atinakosta.github.io/catgraph/reference/modality_gravity.md)
  — computes Modality Gravity Index (MGI+, MGI-, dMGI) and Orbital Score
  (OS) for every node in a catmodgraph. Identifies attractor and
  satellite modalities based on prevalence-weighted edge contributions.

- [`plot_gravity()`](https://atinakosta.github.io/catgraph/reference/plot_gravity.md)
  — 2x3 panel figure comparing traditional centrality (strength,
  betweenness, eigenvector) against gravity indices (MGI+, OS, dMGI) on
  a shared network layout. Set `bars = TRUE` for a bar chart version.

- [`plot_gravity_scatter()`](https://atinakosta.github.io/catgraph/reference/plot_gravity_scatter.md)
  — scatter plot of eigenvector centrality vs dMGI, with contradiction
  cases labelled. Primary diagnostic for demonstrating that MGI captures
  structure beyond standard centrality.

- [`compare_gravity()`](https://atinakosta.github.io/catgraph/reference/compare_gravity.md)
  — compares dMGI and OS profiles across two conditional catmodgraph
  objects (e.g. female vs male subgroup), returning a ranked difference
  table and optional bar chart.

### S3 methods

- [`print.modality_gravity()`](https://atinakosta.github.io/catgraph/reference/print.modality_gravity.md)
  — role-grouped formatted console output (ATTRACTORS / NEUTRAL /
  SATELLITES sections).

- [`summary.modality_gravity()`](https://atinakosta.github.io/catgraph/reference/summary.modality_gravity.md)
  — role counts, per-variable breakdown, and Spearman rho(strength,
  dMGI) diagnostic.

- [`node_centrality.catmodgraph()`](https://atinakosta.github.io/catgraph/reference/node_centrality.catmodgraph.md)
  — extends node_centrality() to accept catmodgraph objects, returning
  traditional centrality measures augmented with all gravity indices in
  a single table.

### Other changes

- [`node_centrality()`](https://atinakosta.github.io/catgraph/reference/node_centrality.md)
  converted to S3 generic to support dispatch to both catgraph and
  catmodgraph objects.

- 69 new tests added in `tests/testthat/test-modality-gravity.R`.

## catgraph 0.8.0

### New features

- **[`plot_modality_difference()`](https://atinakosta.github.io/catgraph/reference/plot_modality_difference.md)**:
  new exported function. Renders the output of
  [`test_modality_edge_differences()`](https://atinakosta.github.io/catgraph/reference/test_modality_edge_differences.md)
  as a single graph with a shared layout. Edge colour encodes the sign
  of the difference `weight_x - weight_y`, width scales with
  \|difference\|, and opacity scales with `-log10(p_adjusted)` (clamped
  at 1e-4). Complements the existing `plot.catmodedgetest()` bar chart:
  use the bar chart for a ranked-list read, this function for a
  network-structural read.

- **[`joint_balance()`](https://atinakosta.github.io/catgraph/reference/joint_balance.md)**:
  new high-level user-facing entry point for cross-group
  categorical-balance diagnostics. Combines (1) per-variable marginal
  chi-square tests (BH-adjusted across variables), (2) pairwise omnibus
  tests from
  [`test_modality_graph_equality()`](https://atinakosta.github.io/catgraph/reference/test_modality_graph_equality.md)
  (Bonferroni-adjusted across group pairs), and (3) edge-wise post-hoc
  from
  [`test_modality_edge_differences()`](https://atinakosta.github.io/catgraph/reference/test_modality_edge_differences.md)
  (BH-adjusted across edges, only for pairs rejecting the omnibus).
  Returns an S3 object of class `jointbalance` with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. The
  plot method is a two-panel layout: marginal bar chart + modality-
  difference graph for the most-rejecting pair.

- **[`build_conditional_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_conditional_modality_graph.md)**:
  new exported function. Thin wrapper over
  [`build_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)
  that subsets the data by one or more conditioning levels (e.g.,
  `list(sex = "female")`) and stores the conditioning specification on
  the returned object. The returned object is a plain `catmodgraph` so
  that all existing methods (plot, prune, cluster, compare) work
  unchanged. Intended for exploring conditional joint categorical
  structure in sociodemographic strata, survey waves, or study sites.
  Does NOT support conditioning on derived cluster labels — the
  conditioning must be an observed variable.

## catgraph 0.7.0

This release adds formal inference for cross-group structural
differences in modality networks. Two permutation-based tests are
introduced, complementing the visual comparison tools shipped in 0.6.0.

### New features

- **[`test_modality_graph_equality()`](https://atinakosta.github.io/catgraph/reference/test_modality_graph_equality.md)**:
  new exported function. Tests the null hypothesis that two modality
  networks are drawn from the same joint distribution, using label
  permutation. Supports three test statistics (Frobenius, Jaccard,
  maximum edge difference), two pipeline modes (unfiltered full weight
  matrix or the user’s pruned analysis), and stratified permutation via
  a `strata` argument. Returns an S3 object of class `catmodtest` with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

- **[`test_modality_edge_differences()`](https://atinakosta.github.io/catgraph/reference/test_modality_edge_differences.md)**:
  new exported function for edge-wise post-hoc testing after a
  significant global test. Computes empirical p-values per edge from the
  same permutation procedure, with Benjamini-Hochberg FDR correction
  applied across edges. Supports restricting the edge set to the union
  of the two input graphs to reduce testing burden. Returns an S3 object
  of class `catmodedgetest` with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

### Documentation

- New vignette section planned in `comparison.Rmd` demonstrating the
  two-step inferential workflow (omnibus test then edge-wise post-hoc)
  on the bundled `survey_health` dataset.
- Simulation study script added to `tests_manual/` documenting Type I
  error and power characteristics.

### Dependencies

- Adds [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html) to
  Imports (via `@importFrom`).

## catgraph 0.6.0

This release refocuses the modality layer as a **category co-association
network** in the tradition of Multiple Correspondence Analysis and
two-mode affiliation networks, removing functionality that over-claimed
the package as a respondent-segmentation tool. It is a **breaking
change** for users of the 0.5.x respondent-profile module.

### Scope refocus (breaking)

- **Removed:** `assign_users_to_profiles()`,
  `profile_modality_clusters()`, the `catmodprofile` S3 class, the
  `$user_cluster` and `$user_scores` list components, and all associated
  S3 methods. These implemented a post-hoc respondent-to-community
  assignment rule that mis-framed the modality layer as a segmentation
  method. The modality layer is now documented honestly as a descriptive
  category co-association map; respondent-level segmentation should use
  `poLCA` or
  [`FactoMineR::HCPC()`](https://rdrr.io/pkg/FactoMineR/man/HCPC.html)
  instead.
- **Renamed:** `profile_modality_clusters()` -\>
  [`summarise_modality_communities()`](https://atinakosta.github.io/catgraph/reference/summarise_modality_communities.md).
  The function now returns an object of class `catmodcommunity` (was
  `catmodprofile`), with list components `community_summary`,
  `community_members`, and `variable_composition`.
- **Documentation rewrite:** the modality-graph section of the vignette
  has been rewritten to frame the layer as a category co-association
  tool rather than a segmentation tool, with an explicit “when to use
  what” decision table and a dedicated caveats subsection for modality
  graphs.

### New features

- **`plot.catmodgraph(signed = TRUE)`**: new argument to colour edges by
  sign of the stored standardised Pearson residual (green = positive /
  attraction, red = negative / repulsion). Edge transparency scales with
  \|residual\|. Default `signed = FALSE` preserves the previous unsigned
  colouring.
- **`plot.catmodgraph(remove_isolates = TRUE)`**: new argument, default
  `TRUE`. Hides degree-0 vertices from the plot without modifying the
  input object. Singletons arising from pruning no longer clutter the
  visualization by default.
- **[`bipartite_modality_graph()`](https://atinakosta.github.io/catgraph/reference/bipartite_modality_graph.md)**:
  new exported function constructing the two-mode (respondent ×
  modality) incidence graph with an `igraph` bipartite `type` attribute.
  Returns a new S3 class `catbipartite` with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
  Complements the unipartite modality graph by preserving raw incidence
  structure.
- **[`compare_catgraphs()`](https://atinakosta.github.io/catgraph/reference/compare_catgraphs.md)**:
  new exported function for multi-panel visual comparison of
  variable-level association networks across groups, populations, or
  time points. Supports pooled / individual / overlay / none pruning
  strategies, with a shared layout anchored on the union graph.
- **[`compare_modality_graphs()`](https://atinakosta.github.io/catgraph/reference/compare_modality_graphs.md)**:
  new exported function for multi-panel comparison of modality-level
  networks. Supports `restrict = "common"` or `"union"` modality-set
  handling and signed-edge visualisation.
- **New vignette `comparison.Rmd`**: worked comparison of `catgraph`’s
  modality layer against Multiple Correspondence Analysis, bipartite
  affiliation networks, and naive co-occurrence projection on the
  bundled `survey_health` dataset. Includes a decision-table for when to
  use which method and an explicit list of where `catgraph` should not
  be overclaimed.
- **`cluster_modalities(signed = TRUE)`**: new argument for signed-aware
  community detection. Drops negative-residual (repulsion) edges before
  running Louvain, so communities are defined by positive co-association
  only. The `std_resid` attribute is retained on the original graph for
  visualization. Documented as a pragmatic adaptation; for principled
  signed-network community detection (Traag & Bruggeman 2009), see the
  `signnet` package.

### Tests

Test suite updated to reflect the removed functions. Tests referencing
`assign_users_to_profiles`, `$user_cluster`, or `$user_scores` have been
deleted; new tests cover
[`summarise_modality_communities()`](https://atinakosta.github.io/catgraph/reference/summarise_modality_communities.md)
and its S3 methods.

## catgraph 0.5.0

This release refocuses the package on its network-native intent and
ships a working modality-level association network module. It is a
**breaking change** for users of the 0.4.x respondent-segmentation
module.

### Scope refocus

The package is now organised around two network modules:

1.  **Variable-level association network** — shipped and stable since
    0.4.0.
2.  **Modality-level association network** — new in 0.5.0. Nodes are
    modalities (factor levels across all variables), edges are
    cross-variable modality associations (phi coefficient, with the
    signed standardised Pearson residual stored as an interpretive edge
    attribute). Louvain or Walktrap communities on that graph serve as
    modality clusters, and respondents are assigned to profiles via an
    argmax rule (count-based default, weighted-score option).

### New features

- **[`build_modality_graph()`](https://atinakosta.github.io/catgraph/reference/build_modality_graph.md)**:
  constructs a modality-level graph from a categorical data frame.
  Same-variable modality pairs are excluded by construction. Returns an
  object of new S3 class `catmodgraph`.
- **[`prune_modality_edges()`](https://atinakosta.github.io/catgraph/reference/prune_modality_edges.md)**:
  modality-level analogue of
  [`prune_edges()`](https://atinakosta.github.io/catgraph/reference/prune_edges.md),
  removing weak or non-significant edges by phi and p-value thresholds.
  Also supports isolate removal and keeps the modality metadata,
  indicator matrix, and membership aligned with the surviving vertices.
- **[`cluster_modalities()`](https://atinakosta.github.io/catgraph/reference/cluster_modalities.md)**:
  Louvain (default) or Walktrap community detection on a `catmodgraph`.
  Writes membership onto vertices.
- **`assign_users_to_profiles()`**: assigns each respondent to a
  modality cluster using either a count rule or a node-strength-weighted
  rule.
- **`profile_modality_clusters()`**: summarises modality communities
  (size, variable composition, internal cohesion, user-cluster counts).
  Returns an object of new S3 class `catmodprofile`.
- **[`plot.catmodgraph()`](https://atinakosta.github.io/catgraph/reference/plot.catmodgraph.md)**:
  dedicated plot method for modality graphs, with colouring by
  originating variable or by detected cluster.
- New S3 classes: `catmodgraph` and `catmodprofile`, with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods.
- New vignette section “Modality-level network analysis” documenting the
  full modality workflow on the expanded Titanic dataset.

### Breaking changes

- **Removed:** `segment_individuals()`, `profile_clusters()`,
  `plot_profiles()`, the `catsegment` and `catprofile` S3 classes, and
  all their S3 methods. The Gower-distance + PAM / hclust / k-modes row
  clustering they implemented was distance-based, not network-based, and
  therefore did not match the package’s network-first intent. The
  modality-network module described above replaces this functionality.
- **Dropped dependencies:** `cluster` and `klaR` are no longer required.

### Cleanup

- Corrected the `@source` tag in the `survey_health` dataset
  documentation (previously referenced a `data-raw/` script that was not
  in the tarball).
- Fixed an escaped-backtick rendering bug in the vignette’s
  `chained-caveat` chunk.
- Fixed a duplicated YAML block at the top of the vignette that rendered
  as visible text in the HTML output.

### Scope change (breaking for documentation only)

- **Package title changed** from *“Undirected Graphical Models for
  Categorical Data Using Effect Size Metrics”* to *“Effect-Size Weighted
  Association Networks for Categorical Data”*.
- **DESCRIPTION and vignette** now state explicitly that a `catgraph` is
  a pairwise *marginal* association network, not a conditional-
  independence graphical model. Edges encode bivariate dependence only.
- **New vignette sections** “Scope and interpretation” and
  “Methodological caveats” document the five caveats reviewers should
  know about: marginal-not-conditional structure, mixed phi/Cramer’s V
  metrics, pairwise-deletion edge-specific sample sizes, multiple
  testing, and unsigned weights.

### Bug fixes

- **[`build_graph()`](https://atinakosta.github.io/catgraph/reference/build_graph.md)
  no longer forces zero weights to `.Machine$double.eps`.** In v0.3.0
  and earlier, pairs with true zero association were stored as edges
  with weight `2.22e-16`, which structurally guaranteed a complete graph
  and silently inflated density, centrality, and community-detection
  measures. True zero weights are now absent edges. For dense
  similarity-matrix output (e.g. heatmaps), use the new
  [`assoc_similarity()`](https://atinakosta.github.io/catgraph/reference/assoc_similarity.md)
  function.
- **[`catgraph()`](https://atinakosta.github.io/catgraph/reference/catgraph.md)
  now stores the processed data, not the raw input.** In v0.3.0,
  `cg$data` held the user’s original data while
  [`build_graph()`](https://atinakosta.github.io/catgraph/reference/build_graph.md)
  internally coerced columns and dropped constant columns. This meant
  [`catgraph_ci()`](https://atinakosta.github.io/catgraph/reference/catgraph_ci.md)
  resampled from a dataset that differed from the one used for the
  original point estimate. `cg$data` now holds the same processed
  object; the original is kept in `cg$raw_data` for reference.
- **[`plot.catgraph()`](https://atinakosta.github.io/catgraph/reference/plot.catgraph.md)
  now honours the `layout` argument.** The igraph branch hard-coded
  Fruchterman-Reingold; `kk`, `circle`, `grid`, `graphopt`, `nicely`,
  and `random` are now respected.

### New features

- **[`assoc_similarity()`](https://atinakosta.github.io/catgraph/reference/assoc_similarity.md)**:
  a new exported function returning the dense `p x p` similarity matrix
  of all pairs (including zero-weight pairs), suitable for heatmaps.
  Separates the similarity-matrix object from the graph-topology object.
- **`prune_edges(..., p_adjust = ...)`**: multiple-testing correction
  across all `choose(p, 2)` simultaneous tests. Default is Benjamini-
  Hochberg (BH); `holm`, `bonferroni`, and `none` are also available.
  Two new edge attributes (`p_value_adj`, `p_adjust_method`) are added.
- **[`compute_assoc()`](https://atinakosta.github.io/catgraph/reference/compute_assoc.md)
  severe-sparsity warning**: a second-tier warning fires when \>20% of
  cells have expected count \< 5, or when a table of minimum dimension
  \>= 3 has observations-to-cells ratio below 5.

### New tests

- `test-build_graph.R` — zero weights are absent edges; `cg$data` is
  processed, not raw.
- `test-prune_edges.R` — multiplicity-correction behaviour.
- `test-layout.R` — layout argument is respected.
- `test-assoc_similarity.R` — dense similarity helper.

## catgraph 0.3.0

Starting version for testing.
