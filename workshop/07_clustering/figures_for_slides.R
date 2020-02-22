library(here)
source(here("workshop/07_clustering/code.R"))

library(ggplot2)
library(cowplot)

# Examples of KNN / SNN graphs and cluster detection ---------------------------

sce.pbmc.subset <- sce.pbmc[, 1:20]

knng <- buildKNNGraph(sce.pbmc.subset, k=1, use.dimred = "PCA")

plot(knng, main = "KNN: k=1")
plot(buildKNNGraph(sce.pbmc.subset, k=2, use.dimred = "PCA"), main = "KNN: k=2")
plot(buildKNNGraph(sce.pbmc.subset, k=3, use.dimred = "PCA"), main = "KNN: k=3")
plot(buildKNNGraph(sce.pbmc.subset, k=19, use.dimred = "PCA"), main = "KNN: k=19 (20 - 1)")

snng <- buildSNNGraph(sce.pbmc.subset, k=1, use.dimred = "PCA")

plot(knng, main = "KNN: k=1")
plot(snng, main = "SNN: k=1")

plot(snng, main = "SNN: k = 1", edge.width = igraph::E(snng)$weight)

g3 <- buildSNNGraph(sce.pbmc.subset, k=3, use.dimred = "PCA")
plot(g3, main = "SNN: k = 3")
plot(g3, main = "SNN: k = 3", edge.width = igraph::E(g3)$weight)

clusters <- igraph::cluster_walktrap(g3)

plot(g3, main = "SNN: k = 3", edge.width = igraph::E(g3)$weight, vertex.color = clusters$membership)

# Naive comparison of scran-style and Seurat-style clustering ------------------

plot_grid(
  plotReducedDim(sce.pbmc, "TSNE", colour_by="cluster") +
    ggtitle("scran-style") +
    theme_cowplot(font_size = 16),
  plotReducedDim(sce.pbmc, "TSNE", colour_by="cluster2") +
    theme_cowplot(font_size = 16) +
    ggtitle("Seurat-style"),
  ncol = 2)
