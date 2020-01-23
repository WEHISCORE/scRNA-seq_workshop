# Library size factors vs. convolution size factors ----------------------------

# Colouring points using the supplied cell-types
plot(lib.sf.zeisel, deconv.sf.zeisel, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(sce.zeisel$level1class)))
abline(a=0, b=1, col="red")
legend("topleft",
       legend = levels(factor(sce.zeisel$level1class)),
       col = seq_along(levels(factor(sce.zeisel$level1class))),
       pch = 16)
