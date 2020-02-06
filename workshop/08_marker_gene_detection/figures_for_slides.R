library(here)
source(here("workshop/08_marker_gene_detection/code.R"))

# Illustration of concept behind using pairwise tests to find marker genes -----

plotExpression(sce.pbmc, features = rownames(sce.pbmc)[1], x = "cluster", colour_by = "cluster")
plotExpression(sce.pbmc, features = rownames(sce.pbmc)[2], x = "cluster", colour_by = "cluster")
i <- which(rownames(sce.pbmc) == "CD3E")
plotExpression(sce.pbmc, features = rownames(sce.pbmc)[i], x = "cluster", colour_by = "cluster")

cd3e_ttest <- pairwiseTTests(
  x = logcounts(sce.pbmc)[i, , drop = FALSE],
  groups = sce.pbmc$cluster)

val <- cbind(
  DataFrame(
    comparison = paste0("cluster 1 vs. ", cd3e_ttest$pairs$second[which(cd3e_ttest$pairs$first == 1)])),
  do.call(rbind, cd3e_ttest$statistics[which(cd3e_ttest$pairs$first == 1)]))
val$logFC <- round(val$logFC, 2)
val$p.value <- signif(val$p.value, 2)
val$FDR <- NULL
val

scran:::combine_simes(val$p.value, logp = FALSE)

val <- cbind(
  DataFrame(
    comparison = paste0("cluster 2 vs. ", cd3e_ttest$pairs$second[which(cd3e_ttest$pairs$first == 2)])),
  do.call(rbind, cd3e_ttest$statistics[which(cd3e_ttest$pairs$first == 2)]))
val$logFC <- round(val$logFC, 2)
val$p.value <- signif(val$p.value, 2)
val$FDR <- NULL
val
scran:::combine_simes(val$p.value, logp = FALSE)

# TODO: Illustrate the idea of findMarkers() standard application by stepping through:
#
# 1. pairwiseTTest().
#   Run a t-test for each gene comparing cluster 1 and 2.
#   Then, repeat for clusters 1 and 3.
#   Then, repeat for cluster 1 and 4.
#   ...
#   Then, repeat for cluster 1 and 18.
#
# 2. combineMarkers() with default arguments, i.e. **one way** of aggregating all these pairwise comparisons
#   Construct ranking
#       The first-ranked gene in each pairwise comparison is given a Top value of 1.
#       The second-ranked gene in each pairwise comparison is given a Top value of 2 (unless it's already been given a Top value of 1).
#       The third-ranked gene in each pairwise comparison is given a Top value of 3 (unless it's already been given a Top value of 1 or 2).
#         ...
#         The last-ranked gene in each pairwise comparison is given a Top value of 'last' (unless it's already been given a Top value of 1, 2, ..., last - 1)
#     Then, compute the Simes-adjusted P-value for each gene
#     Then, compute the FDR value for each gene.

