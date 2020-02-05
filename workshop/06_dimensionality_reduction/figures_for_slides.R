# Effect of perplexity on t-SNE ------------------------------------------------

set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=5)
out5 <- plotReducedDim(sce.zeisel, dimred="TSNE",
                       colour_by="level1class") + ggtitle("perplexity = 5")

set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=20)
out20 <- plotReducedDim(sce.zeisel, dimred="TSNE",
                        colour_by="level1class") + ggtitle("perplexity = 20")

set.seed(100)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=80)
out80 <- plotReducedDim(sce.zeisel, dimred="TSNE",
                        colour_by="level1class") + ggtitle("perplexity = 80")

multiplot(out5, out20, out80, cols=3)

# Effect of n_neighbors on UMAP ------------------------------------------------

set.seed(100)
sce.zeisel <- runUMAP(sce.zeisel, dimred="PCA", n_neighbors=5)
out5 <- plotReducedDim(sce.zeisel, dimred="UMAP",
                       colour_by="level1class") + ggtitle("n_neighbors = 5")

set.seed(100)
sce.zeisel <- runUMAP(sce.zeisel, dimred="PCA", n_neighbors=15)
out15 <- plotReducedDim(sce.zeisel, dimred="UMAP",
                        colour_by="level1class") + ggtitle("n_neighbors = 15")

set.seed(100)
sce.zeisel <- runUMAP(sce.zeisel, dimred="PCA", n_neighbors=30)
out30 <- plotReducedDim(sce.zeisel, dimred="UMAP",
                        colour_by="level1class") + ggtitle("n_neighbors = 30")

multiplot(out5, out15, out30, cols=3)
