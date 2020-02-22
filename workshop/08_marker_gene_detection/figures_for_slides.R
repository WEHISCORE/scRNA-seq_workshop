library(here)
source(here("workshop/08_marker_gene_detection/code.R"))

library(ggplot2)
library(cowplot)

# Illustration of looking at expression of each gene by cluster ----------------

# NOTE: These plots were copied to clipboard at 1500 x 800 resolution.

plotExpression(
  sce.pbmc,
  features = rownames(sce.pbmc)[1],
  x = "cluster",
  colour_by = "cluster") +
  theme_cowplot(font_size = 20)
plotExpression(
  sce.pbmc,
  features = rownames(sce.pbmc)[2],
  x = "cluster",
  colour_by = "cluster") +
  theme_cowplot(font_size = 20)
plotExpression(
  sce.pbmc,
  features = "PTPRC",
  x = "cluster",
  colour_by = "cluster") +
  theme_cowplot(font_size = 20)
plotExpression(
  sce.pbmc,
  features = "PF4",
  x = "cluster",
  colour_by = "cluster") +
  theme_cowplot(font_size = 20)
plotExpression(
  sce.pbmc,
  features = "CD3E",
  x = "cluster",
  colour_by = "cluster") +
  theme_cowplot(font_size = 20)

# Using pairwise t-tests -------------------------------------------------------

cd3e_ttest <- pairwiseTTests(
  x = logcounts(sce.pbmc)["CD3E", , drop = FALSE],
  groups = sce.pbmc$cluster)

val <- cbind(
  DataFrame(
    comparison = paste0("cluster 1 vs. ", cd3e_ttest$pairs$second[which(cd3e_ttest$pairs$first == 1)])),
  do.call(rbind, cd3e_ttest$statistics[which(cd3e_ttest$pairs$first == 1)]))
val$logFC <- round(val$logFC, 2)
val$p.value <- signif(val$p.value, 2)
val$FDR <- NULL
val

val <- cbind(
  DataFrame(
    comparison = paste0("cluster 2 vs. ", cd3e_ttest$pairs$second[which(cd3e_ttest$pairs$first == 2)])),
  do.call(rbind, cd3e_ttest$statistics[which(cd3e_ttest$pairs$first == 2)]))
val$logFC <- round(val$logFC, 2)
val$p.value <- signif(val$p.value, 2)
val$FDR <- NULL
val

val <- cbind(
  DataFrame(
    comparison = paste0("cluster 18 vs. ", cd3e_ttest$pairs$second[which(cd3e_ttest$pairs$first == 18)])),
  do.call(rbind, cd3e_ttest$statistics[which(cd3e_ttest$pairs$first == 2)]))
val$logFC <- round(val$logFC, 2)
val$p.value <- signif(val$p.value, 2)
val$FDR <- NULL
val

pf4_ttest <- pairwiseTTests(
  x = logcounts(sce.pbmc)["PF4", , drop = FALSE],
  groups = sce.pbmc$cluster)

pf4_ttest$statistics[1]
pf4_ttest$statistics[length(pf4_ttest$statistics)]

cd3e_markers <- combineMarkers(cd3e_ttest$statistics, cd3e_ttest$pairs)
combineMarkers(cd3e_ttest$statistics, cd3e_ttest$pairs, pval.type = "any")

# NOTE: The idea of findMarkers() standard application by stepping through:
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

# Ignoring blocking factors can lead to 'false' marker genes -------------------

m.out <- findMarkers(sce.416b, groups=sce.416b$cluster,
                     block=sce.416b$block, test.type="t", direction="up")
demo <- m.out[["1"]]

m.out2 <- m.out <- findMarkers(sce.416b, groups=sce.416b$cluster,
                               test.type="t", direction="up")

demo2 <- m.out2[["1"]]

# Find genes that are 'DE' when not accounting for blocking
i <- which(demo2$FDR < 0.05 & demo[rownames(demo2), ]$FDR > 0.05)

p1 <- plotExpression(sce.416b, features=rownames(demo2)[i[1]],
                     x="cluster", colour_by="cluster") +
  theme_cowplot(font_size = 20) +
  ggtitle("Ignoring experimental blocks")
p2 <- plotExpression(sce.416b, features=rownames(demo2)[i[1]],
                     x="cluster", colour_by="cluster", other_fields = "block") +
  facet_grid(~ block) +
  theme_cowplot(font_size = 20) +
  ggtitle("Accounting for experimental blocks")

plot_grid(p1, p2, ncol = 1)

# Invalidity of P-values -------------------------------------------------------

library(scran)
set.seed(0)
y <- matrix(rnorm(100000), ncol=200)
clusters <- kmeans(t(y), centers=2)$cluster
out <- findMarkers(y, clusters)
hist(out[[1]]$p.value, col="grey80", xlab="p-value", main = "findMarkers() P-values in the absence of any real clusters")
