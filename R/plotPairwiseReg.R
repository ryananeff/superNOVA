#' @title Plot multiple linear regressions for each subgroup for a pair of genes.
#' @description Visualize pairwise gene coexpression correlations for a set of subtypes for two particular genes.
#' @param geneA The first gene to plot.
#' @param geneB The second gene in coexpression with the first gene.
#' @datExpr A gene expression matrix (gene rows by sample columns).
#' @design_mat A model matrix (sample rows by group columns) with 0 or 1 denoting membership in the group.
#' @return A plot of geneA vs geneB expression with subgroup correlations highlighted.
#' @import grDevices graphics
#' @export
plotPairwiseReg <- function(geneA, geneB, datExpr,design_mat){

  x = datExpr[geneA,]
  y = datExpr[geneB,]

  cols = sample(grDevices::rainbow(length(colnames(design_mat))))

  plot(x,y,main="Regression models",xlab = geneA,ylab=geneB,col=cols_design)

  for(group in 1:ncol(design_mat)){
    x = datExpr[geneA,design_mat[,group]==1]
    y = datExpr[geneB,design_mat[,group]==1]
    abline(lm(y ~ x),col=cols[group])
  }
}
