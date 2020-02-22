library(here)
source(here("workshop/09_cell_annotation/code.R"))

library(ggplot2)
library(cowplot)

# Comparing to clustering ------------------------------------------------------

cluster_colours <- setNames(
  scater:::.get_palette("tableau20"),
  levels(sce.pbmc$cluster))
cluster_colours <- cluster_colours[!is.na(names(cluster_colours))]

plotTSNE(sce.pbmc, other_fields = "cluster") +
  geom_point(aes(fill = cluster, colour = cluster)) +
  scale_fill_manual(values = cluster_colours) +
  scale_colour_manual(values = cluster_colours) +
  guides(fill = FALSE, colour = FALSE) +
  theme_cowplot(font_size = 20)

plotScoreHeatmap(pred, clusters = sce.pbmc$cluster, order.by.clusters = TRUE, annotation_colors = list(Clusters = cluster_colours),
                 annotation_legend = FALSE, legend = FALSE)


plotTSNE(sce.pbmc, colour_by = "cluster", text_by = I(abbreviate(pred$labels, minlength = 3)), text_size = 10) +
  scale_fill_manual(values = cluster_colours) +
  scale_colour_manual(values = cluster_colours) +
  guides(fill = FALSE, colour = FALSE) +
  theme_cowplot(font_size = 20)
