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





