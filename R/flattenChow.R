#' @title Convert the output from superNOVA chowCor to an annotated edgelist, suitable for writing out to a .tsv file.
#' @description Takes in results from chowCor as a list of 2d matrices and converts it into a single data.frame, where the first two columns are an edgelist of all pairs compared.
#' @param chow_result A list of matrices as output by chowCor. Required
#' @return Returns a data.frame, sorted by superNOVA p-value, where each row is a feature pair, with columns "Gene1","Gene2","groupCor","groupCorPval","globalCor","globalCorP","pValDiff","Classes".
#' @keywords superNOVA
#' @export
flattenChow = function(chow_result){
  rownames = rownames(chow_result$corrs)
  colnames = colnames(chow_result$corrs)
  corrsflat = cbind(which(!is.na(chow_result$corrs),arr.ind = TRUE),na.omit(as.vector(chow_result$corrs)))
  corrsflat[,1] = sapply(corrsflat[,1],function(x){rownames[as.numeric(x)]})
  corrsflat[,2] = sapply(corrsflat[,2],function(x){colnames[as.numeric(x)]})
  columns = c("Gene1","Gene2","groupCor","groupCorPval","globalCor","globalCorP","pValDiff","qValDiff","Classes","group_order")
  output = data.frame(matrix(NA,nrow=nrow(corrsflat),ncol=length(columns)))
  output[,1:3] = corrsflat
  output[,4] = na.omit(as.vector(chow_result$corrsP))
  output[,5] = na.omit(as.vector(chow_result$globalCor))
  output[,6] = na.omit(as.vector(chow_result$globalCorP))
  output[,7] = as.numeric(na.omit(as.vector(chow_result$pvalues)))
  output[,8] = as.matrix(getQValue(output[,7])$qvalues)
  output[,9] = na.omit(as.vector(chow_result$classes))
  output[,10] = rep(chow_result$groups_compared)
  colnames(output) = columns
  output = output[output$Gene1!=output$Gene2,] #don't output self matching rows
  output = output[order(output$pValDiff,method = "radix"),] #order the rows with most significant first
  return(output)
}
