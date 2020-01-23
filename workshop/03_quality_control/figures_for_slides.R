library(here)
source(here("workshop/03_quality_control/code.R"))

library(ggplot2)
library(cowplot)

# colData summary --------------------------------------------------------------

plot_grid(
  ggplot(as.data.frame(colData(sce.416b))) +
    geom_bar(
      aes(x = block, fill = phenotype),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_brewer(palette = "Dark2") +
    theme_cowplot(font_size = 18) +
    theme(legend.position = "bottom", legend.direction = "vertical"),
  ggplot(as.data.frame(colData(sce.416b))) +
    geom_bar(
      aes(x = block, fill = spike.in.addition),
      position = position_fill(reverse = TRUE)) +
    coord_flip() +
    ylab("Frequency") +
    scale_fill_brewer(palette = "Set1") +
    theme_cowplot(font_size = 18) +
    theme(legend.position = "bottom", legend.direction = "vertical"),
  ggplot(as.data.frame(colData(sce.416b))) +
    geom_bar(aes(x = block, fill = block)) +
    coord_flip() +
    ylab("Number of wells") +
    scale_fill_brewer(palette = "Accent") +
    theme_cowplot(font_size = 18) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      # NOTE: Need extra margin on RHS.
      plot.margin = unit(c(5.5, 16.5, 5.5, 5.5), "pt")),
  align = "h",
  ncol = 3)

# Adaptive thresholds ----------------------------------------------------------

library(gganimate)

anim <- plotColData(sce.416b, y = "sum") +
  scale_y_log10() +
  theme_cowplot(font_size = 16) +
  ggtitle("{next_layer}") +
  stat_summary(fun.y = "median", geom = "point", col = "red", size = 4) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        ymin = median(x) - 0.5 * mad(x),
        y = median(x),
        ymax = median(x) + 0.5 * mad(x))
    },
    geom = "errorbar",
    col = "red",
    size = 1,
    width = 0.1) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        ymin = median(x) - 3 * mad(x),
        y = median(x),
        ymax = median(x) + 3 * mad(x))
    },
    geom = "errorbar",
    col = "red",
    size = 1,
    width = 0.1) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        ymin = median(x) - 3 * mad(x),
        y = median(x),
        ymax = median(x))
    },
    geom = "errorbar",
    col = "red",
    size = 1,
    width = 0.1) +
  geom_hline(
    yintercept = attr(qc.lib2, "thresholds")["lower"],
    lty = 2,
    col = "red") +
  transition_layers(
    layer_length = c(0, 40, 40, 40, 40, 40, 40),
    transition_length = c(0.001, 10, 10, 10, 10, 10, 50),
    from_blank = FALSE,
    layer_names = c(
      "Raw data (optionally log-transformed)",
      "Raw data (optionally log-transformed)",
      "Compute median",
      "Compute MAD",
      "Compute median ± 3 X MAD",
      'Apply type = "lower"',
      "Threshold"),
    keep_layers = c(Inf, Inf, 0, 0, 0, 0, Inf)) +
  enter_fade() +
  exit_fade()

anim_save(here("workshop/03_quality_control/isOutlier.gif"), anim)


# Adaptive outlier detection using a subset ------------------------------------

p1 <- plotColData(sce.grun, x="donor", y="altexps_ERCC_percent",
            colour_by=data.frame(discard=discard.ercc)) +
  ggtitle("Using all samples to compute thresholds") +
  theme_cowplot()
p2 <- plotColData(sce.grun, x="donor", y="altexps_ERCC_percent",
            colour_by=data.frame(discard=discard.ercc2)) +
  ggtitle("Using a subset of samples to compute thresholds") +
  theme_cowplot()
multiplot(p1, p2, cols=2)
