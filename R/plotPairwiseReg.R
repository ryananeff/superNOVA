#' @title Plot multiple linear regressions for each subgroup for a pair of genes.
#' @description Visualize pairwise gene coexpression correlations for a set of subtypes for two particular genes.
#' @param geneA The first gene to plot.
#' @param geneB The second gene in coexpression with the first gene.
#' @param datExpr A gene expression matrix (gene rows by sample columns).
#' @param design_mat A model matrix (sample rows by group columns) with 0 or 1 denoting membership in the group.
#' @param legend_pos The position of the legend on the plot. Default = "topleft".
#' @param title The title of the plot. Default = "Regression models"
#' @return A plot of geneA vs geneB expression with subgroup correlations highlighted.
#' @import grDevices graphics
#' @export
plotPairwiseReg <- function(geneA, geneB, datExpr,design_mat,legend_pos="topleft",title="Regression models"){

  x = as.numeric(datExpr[geneA,])
  y = as.numeric(datExpr[geneB,])

  cols = grDevices::rainbow(length(colnames(design_mat)),alpha = 0.7)
  cols_translate = setNames(cols,colnames(design_mat))
  cols_design = rep(NA,length(design_mat))
  for(group in colnames(design_mat)){
    cols_design[design_mat[,group]==1] = cols_translate[group]
  }

  plot(x,y,main=title,xlab=geneA,ylab=geneB,col=cols_design,pch=16)

  abline(lm(y ~ x),col="black", lty=2) #this is for the global model

  for(group in 1:ncol(design_mat)){
    x = as.numeric(datExpr[geneA,design_mat[,group]==1])
    y = as.numeric(datExpr[geneB,design_mat[,group]==1])
    abline(lm(y ~ x),col=cols[group])
  }

  cols = c("#000000FF",cols)
  legend_names = c("global",colnames(design_mat))
  legend_types = c(2,rep(1,ncol(design_mat)))

  legend(legend_pos,legend=legend_names,col=cols,lty=legend_types)

}
