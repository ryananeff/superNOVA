#' @title Plot ROC curves for superNOVA simulated data observed and analysed results for the two group case with a change in slope.
#' @description Provides a way to plot ROC curves for the two group case. This will not work for 3 or more simulated groups.
#' @param exp_vs_obs_table The expected versus observed table from analyse_sim_results().
#' @param method The name of the differential coexpression method on the plot. Default=superNOVA
#' @param samples The number of samples that were analyzed to generate the results.
#' @param output Should a file with the ROC curve be written?
#' @return The ROC curves for all 6 differentially coexpressed classes (excludes null classes).
#' @importFrom pROC plot.roc
#' @keywords superNOVA simulate plot slope
#' @export
plot_results_newslope <- function(exp_vs_obs_table,method="superNOVA",samples="unknown", output = TRUE){
    response = apply(X=exp_vs_obs_table, MARGIN=1, FUN=function(x){if(x["true_signif"]){x["true_class"]}else{"NonSig"}})
    predictor = as.numeric(exp_vs_obs_table$pvalue)

    order_factor = unique(response)
    response = factor(response, ordered=T,levels=order_factor)
    #predictor = factor(predictor, ordered=T,levels=order_factor)

    cbPalette = c("#000000", "darkgrey", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                  "#0072B2", "#CC79A7")

    if(output){
      pdf(paste0("roc_plot_",method,"_n_",samples,".pdf"))
      roc_p_0 = plot.roc(response,predictor,levels=c("NonSig","+/+"),smooth=F,col = cbPalette[3],main=paste0("ROC plot for ",method,", n=",samples))
      roc_p_n = plot.roc(response,predictor,levels=c("NonSig","-/-"),smooth=F,add=T,col = cbPalette[4])
      legend(0.2,0.5,legend=c("+/+","-/-"),lty=1,col = cbPalette[3:4],lwd=4)
      dev.off()
    }

    roc_p_0 = plot.roc(response,predictor,levels=c("NonSig","+/+"),smooth=F,col = cbPalette[3],main=paste0("ROC plot for ",method,", n=",samples))
    roc_p_n = plot.roc(response,predictor,levels=c("NonSig","-/-"),smooth=F,add=T,col = cbPalette[4])
    legend(0.2,0.5,legend=c("+/+","-/-"),lty=1,col = cbPalette[3:4],lwd=4)

    return(list(roc_pp = roc_p_0,
                roc_nn = roc_p_n))
}
