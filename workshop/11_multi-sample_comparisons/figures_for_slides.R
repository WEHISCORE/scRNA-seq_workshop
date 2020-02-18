library(here)
source(here("workshop/11_multi-sample_comparisons/code.R"))

library(ggplot2)
library(cowplot)

# Illustrative dataset: chimeric mouse embryos ---------------------------------

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

# Discarding ambient DEGs ------------------------------------------------------

library(MouseGastrulationData)
raw.tal1 <- Tal1ChimeraData(type="raw")

library(Matrix)
ambient <- vector("list", length(raw.tal1))
for (i in seq_along(raw.tal1)) {
  curmat <- counts(raw.tal1[[i]])
  is.empty <- colSums(curmat) < 100
  ambient[[i]] <- rowSums(curmat[,is.empty])
}

ambient <- do.call(cbind, ambient)
colnames(ambient) <- names(raw.tal1)
rownames(ambient) <- uniquifyFeatureNames(
  rowData(raw.tal1[[1]])$ENSEMBL, rowData(raw.tal1[[1]])$SYMBOL)
head(ambient)

is.hbb <- grep("^Hb[ab]-", rownames(summed.neural))
neural.hb <- colSums(counts(summed.neural)[is.hbb,])
ambient.hb <- colSums(ambient[is.hbb,])
scaled.ambient <- t(t(ambient) * neural.hb/ambient.hb)
head(scaled.ambient)

ratio <- rowMeans(scaled.ambient) / rowMeans(counts(summed.neural))
non.ambient <- ratio < 0.1
summary(non.ambient)

okay.genes <- names(non.ambient)[which(non.ambient)]
res.neural2 <- res.neural[rownames(res.neural) %in% okay.genes,]
summary(decideTests(res.neural2))

topTags(res.neural2)

# Subtracting ambient counts ---------------------------------------------------

subtracted <- counts(summed.neural) - scaled.ambient
subtracted <- round(subtracted)
subtracted[subtracted < 0] <- 0
subtracted[is.hbb,]

# Re-using keep.neural to simplify comparison.
y.ambient <- DGEList(ambient)
y.ambient <- y.ambient[keep.neural,]
y.ambient <- calcNormFactors(y.ambient)
y.ambient <- estimateDisp(y.ambient, design)
fit.ambient <- glmQLFit(y.ambient, design, robust=TRUE)
res.ambient <- glmQLFTest(fit.ambient, coef=ncol(design))
summary(decideTests(res.ambient))

topTags(res.ambient, n=10)

