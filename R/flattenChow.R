#' @title Convert the output from superNOVA chowCor to an annotated edgelist, suitable for writing out to a .tsv file.
#' @description Takes in results from chowCor as a list of 2d matrices and converts it into a single data.frame, where the first two columns are an edgelist of all pairs compared.
#' @param chow_result A list of matrices as output by chowCor. Required
#' @param adjust_q Should the q-values be calculated? Otherwise, qValDiff will be the same as pValDiff
#' @param sort_output Should the final output be sorted from smallest to largest p-values?
#' @return Returns a data.frame, sorted by superNOVA p-value, where each row is a feature pair, with columns "Gene1","Gene2","groupCor","groupCorPval","globalCor","globalCorP","pValDiff","Classes".
#' @keywords superNOVA
#' @export
flattenChow = function(chow_result, method="BH", sort_output=T){
  rownames = rownames(chow_result$corrs)
  colnames = colnames(chow_result$corrs)
  corrsflat = cbind(which(!is.na(chow_result$corrs),arr.ind = TRUE),na.omit(as.vector(chow_result$corrs)))
  corrsflat[,1] = sapply(corrsflat[,1],function(x){rownames[as.numeric(x)]})
  corrsflat[,2] = sapply(corrsflat[,2],function(x){colnames[as.numeric(x)]})
  columns = c("Gene1","Gene2","groupCor","groupSlope","groupCorPval", #1-5
              "globalCor","globalSlope","globalCorP","pValDiff","pValDiff_adj", #6-10
              "Classes","group_order") #11-12
  output = data.frame(matrix(NA,nrow=nrow(corrsflat),ncol=length(columns)))
  output[,1:3] = corrsflat
  output[,4] = na.omit(as.vector(chow_result$slopes))
  output[,5] = na.omit(as.vector(chow_result$corrsP))
  if(rownames(chow_result$globalCor)[1]==colnames(chow_result$globalCor)[1]){
  output[,6] = as.vector(chow_result$globalCor[upper.tri(chow_result$globalCor,diag=TRUE)])
  output[,7] = as.vector(chow_result$globalSlope[upper.tri(chow_result$globalSlope,diag=TRUE)])
  output[,8] = as.vector(chow_result$globalCorP[upper.tri(chow_result$globalCorP,diag=TRUE)])
  } else {
  	output[,6] = as.vector(chow_result$globalCor)
  	output[,7] = as.vector(chow_result$globalSlope)
  	output[,8] = as.vector(chow_result$globalCorP)
  }
  output[,9] = as.numeric(na.omit(as.vector(chow_result$pvalues)))
  if (method=="qvalue"){
    output[,10] = as.matrix(getQValue(output[,9])$qvalues)
  }
  else if (method=="none"){
    output[,10] = output[,9]
  }
  else {
    output[,10] = p.adjust(output[,9],method=method)
  }
  output[,11] = na.omit(as.vector(chow_result$classes))
  output[,12] = rep(chow_result$groups_compared)
  colnames(output) = columns
  output = output[output$Gene1!=output$Gene2,] #don't output self matching rows
  output = output[order(output$pValDiff,method = "radix"),] #order the rows with most significant first
  return(output)
}
