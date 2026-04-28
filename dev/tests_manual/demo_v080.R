# ============================================================================
# catgraph v0.8.0 — end-to-end demonstration
#
# Constructs a synthetic two-group dataset with:
#   - IDENTICAL per-variable marginals across the two groups
#   - DIFFERENT joint structure (group A has a designed-in A<->B association
#     that is absent in group B)
#
# Expected behaviour from a correctly-working joint_balance():
#   - Marginal tests: NONE should reject at alpha = 0.05 (after BH)
#   - Pairwise omnibus: SHOULD reject at alpha = 0.05 (after Bonferroni)
#   - Edge-wise post-hoc: should flag the designed-in A<->B modality pairs
#
# This is a sanity check, not a calibration study. One run, one seed.
# ============================================================================

library(catgraph)
library(igraph)

set.seed(20260424)

# ---- 1. Construct the dataset ---------------------------------------------
#
# Variables:
#   V1 ("A"/"B"):   binary, 50/50 in both groups
#   V2 ("P"/"Q"):   binary, 50/50 in both groups
#   V3 ("M"/"N"):   binary, 50/50 in both groups (noise variable, unrelated)
#   group ("g1"/"g2"): grouping variable
#
# In group g1: V1 and V2 are strongly positively associated.
#   P(V2 = P | V1 = A) = 0.85;  P(V2 = P | V1 = B) = 0.15
# In group g2: V1 and V2 are independent.
#   P(V2 = P | V1 = A) = 0.50;  P(V2 = P | V1 = B) = 0.50
#
# By construction:
#   - Marginal P(V1 = A) = 0.5 in both groups
#   - Marginal P(V2 = P) = 0.5 in both groups  (symmetry of the conditional probs)
#   - V3 is independent of everything in both groups
# So every marginal chi-square should be null; only the JOINT structure differs.

make_group <- function(n, p_v2_given_v1_A, p_v2_given_v1_B) {
  V1 <- sample(c("A", "B"), n, replace = TRUE, prob = c(0.5, 0.5))
  V2 <- character(n)
  for (i in seq_len(n)) {
    if (V1[i] == "A") {
      V2[i] <- sample(c("P", "Q"), 1L,
                      prob = c(p_v2_given_v1_A, 1 - p_v2_given_v1_A))
    } else {
      V2[i] <- sample(c("P", "Q"), 1L,
                      prob = c(p_v2_given_v1_B, 1 - p_v2_given_v1_B))
    }
  }
  V3 <- sample(c("M", "N"), n, replace = TRUE, prob = c(0.5, 0.5))
  data.frame(V1 = V1, V2 = V2, V3 = V3, stringsAsFactors = FALSE)
}

n_per_group <- 400

df_g1 <- make_group(n_per_group, p_v2_given_v1_A = 0.85, p_v2_given_v1_B = 0.15)
df_g1$group <- "g1"

df_g2 <- make_group(n_per_group, p_v2_given_v1_A = 0.50, p_v2_given_v1_B = 0.50)
df_g2$group <- "g2"

dat <- rbind(df_g1, df_g2)
dat[] <- lapply(dat, factor)

cat("Dataset dimensions:", nrow(dat), "rows,", ncol(dat), "columns\n\n")
cat("Marginal distributions by group (should look identical):\n")
print(prop.table(table(dat$V1, dat$group), margin = 2))
print(prop.table(table(dat$V2, dat$group), margin = 2))
print(prop.table(table(dat$V3, dat$group), margin = 2))

cat("\nJoint V1-V2 distribution within g1 (should show strong association):\n")
print(prop.table(table(dat$V1[dat$group == "g1"],
                       dat$V2[dat$group == "g1"])))
cat("\nJoint V1-V2 distribution within g2 (should look ~independent):\n")
print(prop.table(table(dat$V1[dat$group == "g2"],
                       dat$V2[dat$group == "g2"])))


# ---- 2. Variable-level catgraph on pooled data ----------------------------
#
# This is the overall, marginalised view. We expect to see an edge between
# V1 and V2 (the pooled data still shows the association, because g1's signal
# dominates), but it will be weaker than the g1-only edge.

cat("\n=== Variable-level catgraph (pooled) ===\n")
cg <- catgraph(dat[, c("V1", "V2", "V3")])
print(cg)
summary(cg)

cat("\nPruned (BH at 0.05, min effect 0.1):\n")
cg_pruned <- prune_edges(cg, min_weight = 0.10, max_p = 0.05)
print(cg_pruned)


# ---- 3. Modality graphs per group -----------------------------------------

cat("\n=== Per-group modality graphs ===\n")
mg_g1 <- build_modality_graph(subset(dat, group == "g1")[, c("V1", "V2", "V3")])
mg_g2 <- build_modality_graph(subset(dat, group == "g2")[, c("V1", "V2", "V3")])

cat("g1 modality graph:\n"); print(mg_g1)
cat("\ng2 modality graph:\n"); print(mg_g2)

cat("\nEdge weights by graph (V1=A <-> V2=P should be strong in g1, weak in g2):\n")
print(igraph::as_data_frame(mg_g1$graph, what = "edges")[, c("from", "to", "weight")])
print(igraph::as_data_frame(mg_g2$graph, what = "edges")[, c("from", "to", "weight")])


# ---- 4. Side-by-side comparison -------------------------------------------

cat("\n=== compare_modality_graphs() side-by-side panel ===\n")
cat("Rendering... (two panels, shared layout)\n")
compare_modality_graphs(
  list(g1 = mg_g1, g2 = mg_g2),
  restrict = "common",
  signed   = TRUE
)


# ---- 5. joint_balance() full diagnostic -----------------------------------

cat("\n=== joint_balance(): full cross-group diagnostic ===\n")
jb <- joint_balance(
  dat,
  group       = "group",
  n_perm      = 500,
  n_perm_edge = 500,
  seed        = 1,
  verbose     = TRUE
)
print(jb)

cat("\n--- Expected outcomes given how we built the data ---\n")
cat("  Marginal rejections  : 0 of 3  (variables are designed to be marginally balanced)\n")
cat("  Omnibus rejections   : 1 of 1  (the joint structure differs by construction)\n")
cat("  Edge-wise rejections : edges involving V1=*/V2=* should dominate\n")


# ---- 6. Inspect the edgewise post-hoc result ------------------------------

if (length(jb$pairwise_edgewise) > 0L) {
  et <- jb$pairwise_edgewise[[1L]]
  cat("\n=== Edge-wise post-hoc table (top rows) ===\n")
  print(head(et$edge_table, 10), digits = 3)
}


# ---- 7. Single-graph difference plot --------------------------------------

if (length(jb$pairwise_edgewise) > 0L) {
  cat("\n=== plot_modality_difference() graph view ===\n")
  cat("Rendering... (edges: colour = sign of g1-g2 difference,\n")
  cat("                     width = |difference|, alpha = -log10 p_adj)\n")
  plot_modality_difference(
    jb$pairwise_edgewise[[1L]],
    reference    = list(mg_g1, mg_g2),
    group_labels = c("stronger in g1", "stronger in g2")
  )
}


# ---- 8. Full diagnostic plot ----------------------------------------------

cat("\n=== plot(jb): two-panel joint-balance diagnostic ===\n")
cat("Left:  -log10 adjusted marginal p-values (all should be short / grey)\n")
cat("Right: modality-difference graph for the g1 vs g2 pair\n")
plot(jb)


# ---- 9. Conditional modality graph demo ----------------------------------

cat("\n=== build_conditional_modality_graph(): conditioning on V3 = M ===\n")
mg_cond <- build_conditional_modality_graph(
  dat,
  given = list(V3 = "M")
)
print(mg_cond)
cat("\nConditioning metadata:\n")
str(mg_cond$conditioning)

cat("\nNote: V3 is dropped from the graph by construction (it was the\n")
cat("conditioning variable). The remaining graph shows V1-V2 associations\n")
cat("among the V3 = M subset.\n")


# ---- 10. Summary ----------------------------------------------------------

cat("\n\n=== Demo complete ===\n")
cat("All step-function outputs generated. Visually inspect the plots to confirm:\n")
cat(" - compare_modality_graphs panel: g1 shows a V1=A <-> V2=P edge, g2 does not\n")
cat(" - plot_modality_difference: one dominant coloured edge between V1 and V2\n")
cat(" - plot(jb): marginal bar chart is uniformly low, joint graph shows the difference\n")