# catgraph

`catgraph` provides **weighted association networks for categorical data**
at two complementary levels.

1. **Variable-level association network.** Variables are nodes. Edges are
   weighted by one of four association paradigms — frequentist
   (classical and bias-corrected Cramér's V), information-theoretic
   (Normalised Mutual Information and Adjusted Mutual Information), or
   Bayesian (Dirichlet-smoothed Cramér's V) — selected via the `method`
   argument. The workflow supports structural exploration, edge pruning with
   multiple-testing adjustment, bootstrap confidence intervals, and
   descriptive network summaries such as centrality and community structure.

2. **Modality-level co-association network.** Modalities (factor levels)
   are nodes. Cross-variable edges are weighted by the association measure
   chosen at the variable level, with signed standardised Pearson residuals
   stored separately to indicate whether co-occurrence is above or below
   independence expectation. The workflow supports edge pruning, signed edge
   visualisation, and community detection over modalities.

   The modality layer sits in the tradition of Multiple Correspondence
   Analysis and two-mode affiliation networks.
   It operates on pairwise associations and does not model higher-order
   interactions or conditional dependencies. It is a descriptive
   **category co-association map**, not a respondent-segmentation tool.
   For respondent segmentation use
   [`poLCA`](https://cran.r-project.org/package=poLCA) or
   [`FactoMineR::HCPC`](https://cran.r-project.org/package=FactoMineR).

3. **Modality gravity indices.** A novel extension to standard graph
   centrality that incorporates the empirical prevalence of each modality.
   The **Modality Gravity Index (MGI)** and **Orbital Score (OS)** identify
   which modalities act as gravitational attractors (dominant, pulling rarer
   modalities toward them) and which are satellites (minority modalities
   orbiting more prevalent ones). This addresses a fundamental limitation of
   standard centrality indices, which treat all nodes as exchangeable
   regardless of their empirical frequency.

## Installation

```r
# install.packages("remotes")
remotes::install_github("AtinaKosta/catgraph")
```

## What's new in 0.11.0

- **Four association paradigms** are now supported as edge weights via the
  new `method` argument in `catgraph()`, `build_graph()`,
  `assoc_similarity()`, and `build_modality_graph()`:

  | `method` | Paradigm | Description |
  |---|---|---|
  | `"cramers_v"` | Frequentist | Classical phi / Cramér's V (default) |
  | `"cramers_v_corrected"` | Frequentist | Bias-corrected via Bergsma (2013) |
  | `"nmi"` | Information-theoretic | Normalised Mutual Information |
  | `"ami"` | Information-theoretic | Adjusted MI, corrects for chance |
  | `"bayesian_cramers_v"` | Bayesian | Dirichlet-smoothed Cramér's V |

- **Two new exported functions**: `nmi_assoc()` and `bayesian_cramers_v()`
  for computing information-theoretic and Bayesian association on arbitrary
  vector pairs.

- **`alpha` argument** controls the Dirichlet prior concentration for
  `method = "bayesian_cramers_v"` (default `0.5`, Jeffreys prior).

- **Backward compatible**: all existing code using `catgraph(df)`,
  `catgraph(df, corrected = TRUE)`, or `build_graph(df)` works without
  any changes.

## Quick start

```r
library(catgraph)
data(survey_health)

# --- Variable-level: choose your association paradigm ---

# Frequentist (default)
cg <- catgraph(survey_health, method = "cramers_v_corrected")

# Information-theoretic
cg_nmi <- catgraph(survey_health, method = "nmi")
cg_ami <- catgraph(survey_health, method = "ami")

# Bayesian
cg_bay <- catgraph(survey_health, method = "bayesian_cramers_v", alpha = 0.5)

# Prune and plot
cg_p <- prune_edges(cg, min_weight = 0.05, max_p = 0.05, p_adjust = "BH")
plot(cg_p)

# --- Modality-level: same method flows through ---
mg <- build_modality_graph(survey_health, method = "cramers_v_corrected")
mg <- prune_modality_edges(mg, min_weight = 0.10, max_p = 0.05)
mg <- cluster_modalities(mg)
plot(mg, color_by = "cluster", signed = TRUE)

# --- Gravity indices: which modalities are attractors vs satellites? ---
grav <- modality_gravity(mg)
print(grav)                       # role-grouped: ATTRACTORS / SATELLITES
summary(grav)                     # role counts, Spearman rho diagnostic
plot_gravity(mg)                  # 6-panel: traditional centrality vs MGI
plot_gravity_scatter(grav, mg)    # eigenvector vs dMGI contradiction plot

# --- Compare gravity profiles across subgroups ---
mg_f <- build_conditional_modality_graph(survey_health, given = list(sex = "female"))
mg_m <- build_conditional_modality_graph(survey_health, given = list(sex = "male"))
mg_f <- prune_modality_edges(mg_f, min_weight = 0.10, max_p = 0.05)
mg_m <- prune_modality_edges(mg_m, min_weight = 0.10, max_p = 0.05)
compare_gravity(list(female = mg_f, male = mg_m))

# --- Formal test: do two groups differ in overall association structure? ---
test_modality_graph_equality(mg_f, mg_m, n_perm = 500)
```

See `vignette("introduction", package = "catgraph")` for the full
variable-level and modality-level workflow, and
`vignette("comparison", package = "catgraph")` for a worked comparison of
`catgraph`'s modality layer with MCA, bipartite affiliation networks, and
naive co-occurrence projection.

## Scope — use `catgraph` for

| Question | Tool |
|---|---|
| Which categorical variables co-vary pairwise? | `catgraph()` |
| Which association paradigm is most stable? | `catgraph(method = ...)` |
| Which category levels bundle together across variables? | `build_modality_graph()` |
| Which modalities are structurally dominant vs peripheral? | `modality_gravity()` |
| How do category-level association patterns differ across groups? | `compare_*_graphs()` |
| What's the unprojected respondent-modality incidence like? | `bipartite_modality_graph()` |

## Scope — do not use `catgraph` for

- Causal inference or claims about direct effects.
- Conditional-independence structure — use a graphical model instead
  (`bnlearn`, `gRim`).
- Respondent segmentation or latent-class analysis — use `poLCA` or
  `FactoMineR::HCPC()`.
- Estimating conditional relationships between variables — edges are
  marginal and do not control for other variables.
- Mixed continuous/categorical data — use `linkspotter` or `corrgrapher`
  instead.

## Interpretation

All graphs produced by `catgraph` represent **marginal association
structure**. Edge weights quantify pairwise dependence under the chosen
paradigm and should be interpreted descriptively. High rank agreement
across paradigms (Spearman ρ > 0.90) indicates robust structural
conclusions; divergence reveals cases where the choice of metric matters
substantively.

Gravity indices (MGI, OS) incorporate empirical modality prevalence.
A positive dMGI indicates a modality that exerts net gravitational pull
over less prevalent neighbours; a negative dMGI indicates a satellite
modality being pulled toward more prevalent ones.

## License

GPL (>= 3)