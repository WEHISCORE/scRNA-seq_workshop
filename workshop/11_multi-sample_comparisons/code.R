# Illustrative dataset: chimeric mouse embryos ---------------------------------

library(MouseGastrulationData)
sce.chimera <- WTChimeraData(samples=5:10)

# feature-annotation
library(scater)
rownames(sce.chimera) <- uniquifyFeatureNames(
  rowData(sce.chimera)$ENSEMBL, rowData(sce.chimera)$SYMBOL)

# quality-control (already performed by authors)
# We also remove cells labelled as stripped nuclei or doublets
drop <- sce.chimera$celltype.mapped %in% c("stripped", "Doublet")
sce.chimera <- sce.chimera[,!drop]

# normalization (using pre-computed size factors)
sce.chimera <- logNormCounts(sce.chimera)

# variance-modelling
# Retaining all genes with a positive biological component, to
# preserve as much signal as possible in a very heterogeneous
# dataset
library(scran)
dec.chimera <- modelGeneVar(sce.chimera, block=sce.chimera$sample)
chosen.hvgs <- dec.chimera$bio > 0

# merging (hierarchical merge, first merge together reps with the # same genotype, then merge samples across different genotypes
library(batchelor)
set.seed(01001001)
merged <- correctExperiments(sce.chimera,
                             batch=sce.chimera$sample,
                             subset.row=chosen.hvgs,
                             PARAM=FastMnnParam(
                               merge.order=list(
                                 list(1,3,5), # WT (3 replicates)
                                 list(2,4,6)  # td-Tomato (3 replicates)
                               )
                             )
)

# clustering
g <- buildSNNGraph(merged, use.dimred="corrected")
clusters <- igraph::cluster_louvain(g)
merged$cluster <- factor(clusters$membership)

# dimensionality-reduction
# Using an external algorithm to compute nearest
# neighbors for greater speed
merged <- runTSNE(merged, dimred="corrected",
                  external_neighbors=TRUE)
merged <- runUMAP(merged, dimred="corrected",
                  external_neighbors=TRUE)

# Creating pseudo-bulk samples -------------------------------------------------

# Using 'label' and 'sample' as our two factors; each
# column of the output corresponds to one unique
# combination of these two factors.
summed <- aggregateAcrossCells(merged,
                               id=DataFrame(
                                 label=merged$celltype.mapped,
                                 sample=merged$sample))

label <- "Mesenchyme"
current <- summed[,label==summed$celltype.mapped]

# Creating up a DGEList object for use in edgeR:
library(edgeR)
y <- DGEList(counts(current), samples=colData(current))

# Discard sample-labels with very small library sizes
discarded <- isOutlier(y$samples$lib.size, log=TRUE,
                       type="lower")
y <- y[,!discarded]

# Remove genes that are lowly expressed
keep <- filterByExpr(y, group=current$tomato)
y <- y[keep,]

# Compute normalization factors for the pseudo-bulk
# samples (these are different to the single-cell
# size factors)
y <- calcNormFactors(y)

# Construct design matrix, blocking on pool and
# including term for td-Tomato-status
design <- model.matrix(
  ~0 + factor(pool) + factor(tomato),
  y$samples)
# Tidy up design matrix column names
colnames(design) <- gsub(
  "factor|\\(|\\)|TRUE", "", colnames(design))

# Estimate and plot the negative binomial dispersions
y <- estimateDisp(y, design)
plotBCV(y)

# Estimate and plot the quasi-likelihood dispersions
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

# Test for DE due to td-Tomato-status
res <- glmQLFTest(fit, coef="tomato")
# Summarise the results
summary(decideTests(res))

# Putting it all together ------------------------------------------------------

de.results <- list()

for (i in unique(summed$celltype.mapped)) {
  current <- summed[,i==summed$celltype.mapped]
  y <- DGEList(counts(current), samples=colData(current))

  discarded <- isOutlier(
    colSums(counts(current)), log=TRUE, type="lower")
  y <- y[,!discarded]
  y <- y[filterByExpr(y, group=current$tomato),]

  y <- calcNormFactors(y)

  design <- try(
    model.matrix(
      ~0 + factor(pool) + factor(tomato), y$samples),
    silent=TRUE)
  if (is(design, "try-error") ||
      qr(design)$rank==nrow(design) ||
      qr(design)$rank < ncol(design)) {
    next
  }
  colnames(design) <- gsub(
    "factor|\\(|\\)|TRUE", "", colnames(design))

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- glmQLFTest(fit, coef="tomato")
  de.results[[i]] <- res
}

# Tabulate DEGs per label
summaries <- lapply(
  de.results,
  FUN=function(x) summary(decideTests(x))[,1])
sum.tab <- do.call(rbind, summaries)

# Find genes consistently DE across labels
degs <- lapply(
  de.results,
  FUN=function(x) rownames(topTags(x, p.value=0.05)))
common.degs <- sort(table(unlist(degs)),
                    decreasing=TRUE)

# List labels that were skipped due to
# insufficient replicates or contrasts
setdiff(unique(summed$celltype.mapped),
        names(summaries))

# Differential abundance: overview ---------------------------------------------

abundances <- table(merged$celltype.mapped, merged$sample)
abundances <- unclass(abundances)

# Performing te DA analysis ----------------------------------------------------

# Attaching some column metadata.
extra.info <- colData(merged)[match(
  colnames(abundances), merged$sample),]
y.ab <- DGEList(abundances, samples=extra.info)

keep <- filterByExpr(y.ab, group=y.ab$samples$tomato)
y.ab <- y.ab[keep,]

# Construct design matrix, blocking on pool and
# including term for td-Tomato-status
design <- model.matrix(
  ~0 + factor(pool) + factor(tomato),
  y.ab$samples)
# Tidy up design matrix column names
colnames(design) <- gsub(
  "factor|\\(|\\)|TRUE", "", colnames(design))

# Estimate and plot the negative binomial dispersions
y.ab <- estimateDisp(y.ab, design, trend="none")
plotBCV(y.ab, cex=1)

# Estimate and plot the quasi-likelihood dispersions
fit.ab <- glmQLFit(y.ab, design, robust=TRUE,
                   abundance.trend=FALSE)
plotQLDisp(fit.ab, cex=1)

# Test for DA due to td-Tomato-status
res <- glmQLFTest(fit.ab, coef="tomato")
# Summarise the results
summary(decideTests(res))
# View the top DE genes
topTags(res)

# Handling composition effects: Assuming most labels don't change --------------

y.ab2 <- calcNormFactors(y.ab)
y.ab2 <- estimateDisp(y.ab2, design, trend="none")
fit.ab2 <- glmQLFit(y.ab2, design, robust=TRUE,
                    abundance.trend=FALSE)
res2 <- glmQLFTest(fit.ab2, coef=ncol(design))
topTags(res2)

# Handling composition effects: Removing the offending labels ------------------

offenders <- "ExE ectoderm"
y.ab3 <- y.ab[setdiff(
  rownames(y.ab), offenders),, keep.lib.sizes=FALSE]

y.ab3 <- estimateDisp(y.ab3, design, trend="none")
fit.ab3 <- glmQLFit(y.ab3, design, robust=TRUE,
                    abundance.trend=FALSE)
res3 <- glmQLFTest(fit.ab3, coef=ncol(design))
topTags(res3)

# Illustrative dataset: chimeric mouse embryos (slight return) -----------------

library(MouseGastrulationData)
sce.tal1 <- Tal1ChimeraData()
rownames(sce.tal1) <- uniquifyFeatureNames(
  rowData(sce.tal1)$ENSEMBL,
  rowData(sce.tal1)$SYMBOL)

# Using 'label' and 'sample' as our two factors; each
# column of the output corresponds to one unique
# combination of these two factors.
summed.tal1 <- aggregateAcrossCells(sce.tal1,
                                    ids=DataFrame(
                                      label=sce.tal1$celltype.mapped,
                                      sample=sce.tal1$sample))
# Focusing on cells labelled as 'neural crest'
summed.neural <- summed.tal1[
  ,summed.tal1$label=="Neural crest"]

# Standard edgeR analysis, as described earlier
y.neural <- DGEList(counts(summed.neural),
                    samples=colData(summed.neural))
keep.neural <- filterByExpr(y.neural,
                            group=y.neural$samples$tomato)
y.neural <- y.neural[keep.neural,]
y.neural <- calcNormFactors(y.neural)

# Samples were processed in 2 blocks
y.neural$samples$block <- c(1, 2, 1, 2)
# Construct design matrix, blocking on pool and
# including term for td-Tomato-status
design <- model.matrix(
  ~0 + factor(block) + factor(tomato),
  y.neural$samples)
# Tidy up design matrix column names
colnames(design) <- gsub(
  "factor|\\(|\\)|TRUE", "", colnames(design))

y.neural <- estimateDisp(y.neural, design)
fit.neural <- glmQLFit(y.neural, design, robust=TRUE)
res.neural <- glmQLFTest(fit.neural, coef="tomato")

summary(decideTests(res.neural))
topTags(res.neural, n=10)




