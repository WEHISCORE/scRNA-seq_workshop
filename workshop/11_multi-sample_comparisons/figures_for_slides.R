library(here)
source(here("workshop/11_multi-sample_comparisons/code.R"))

library(ggplot2)
library(cowplot)

# Illustrative dataset: chimeric mouse embryos ---------------------------------


# TODO: Move this to figures_for_slides.R and spruce it up by faceting and perhaps adding pre-annoated cell labels

plot_grid(
  plotTSNE(merged, colour_by="tomato", text_by="cluster") +
    theme_cowplot(font_size = 20),
  plotTSNE(merged, colour_by=data.frame(pool=factor(merged$pool))) +
    theme_cowplot(font_size = 20),
  ncol = 2)

cluster_colours <- setNames(
  Polychrome::alphabet.colors(),
  sort(unique(merged$cluster)))
celltype.mapped_colours <- setNames(
  Polychrome::palette36.colors(),
  sort(unique(merged$celltype.mapped)))

plot_grid(
  plotTSNE(merged, colour_by="cluster", text_by="cluster", text_size=10) +
    theme_cowplot(font_size = 20) +
    scale_fill_manual(values = cluster_colours) +
    guides(fill = FALSE) +
    ggtitle("Coloured by cluster"),
  plotTSNE(merged, colour_by="celltype.mapped") +
    theme_cowplot(font_size = 20) +
    scale_fill_manual(values = celltype.mapped_colours) +
    guides(fill = FALSE) +
    ggtitle("Coloured by provided cell type labels"),
  ncol = 2)
