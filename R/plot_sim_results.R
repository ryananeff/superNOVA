#' @title Plot ROC curves for superNOVA simulated data observed and analysed results for the two group case.
#' @description Provides a way to plot ROC curves for the two group case. This will not work for 3 or more simulated groups.
#' @param exp_vs_obs_table The expected versus observed table from analyse_sim_results().
#' @param method The name of the differential coexpression method on the plot. Default=superNOVA
#' @param samples The number of samples that were analyzed to generate the results.
#' @param output Should a file with the ROC curve be written?
#' @return The ROC curves for all 6 differentially coexpressed classes (excludes null classes).
#' @keywords superNOVA simulate plot
#' @export
plot_sim_results <- function(exp_vs_obs_table,method="superNOVA",samples="unknown",output=T){

  response = apply(X=exp_vs_obs_table, MARGIN=1, FUN=function(x){if(x["true_signif"]){x["true_class"]}else{"NonSig"}})
  predictor = as.numeric(exp_vs_obs_table$pValDiff)

  order_factor = unique(response)
  response = factor(response, ordered=T,levels=order_factor)

  cbPalette = c("#000000", "darkgrey", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#CC79A7")
  if(output){
    grDevices::pdf(paste0("roc_plot_",method,"_n_",samples,".pdf"))
    roc_p_0 = pROC::plot.roc(response,predictor,levels=c("NonSig","+/0"),smooth=F,col = cbPalette[3],main=paste0("ROC plot for ",method,", n=",samples))
    roc_p_n = pROC::plot.roc(response,predictor,levels=c("NonSig","+/-"),smooth=F,add=T,col = cbPalette[4])
    roc_0_p = pROC::plot.roc(response,predictor,levels=c("NonSig","0/+"),smooth=F,add=T,col = cbPalette[5])
    roc_0_n = pROC::plot.roc(response,predictor,levels=c("NonSig","0/-"),smooth=F,add=T,col = cbPalette[6])
    roc_n_p = pROC::plot.roc(response,predictor,levels=c("NonSig","-/+"),smooth=F,add=T,col = cbPalette[7])
    roc_n_0 = pROC::plot.roc(response,predictor,levels=c("NonSig","-/0"),smooth=F,add=T,col = cbPalette[8])
    graphics::legend(0.2,0.5,legend=c("+/0","+/-","0/+","0/-","-/+","-/0"),lty=1,col = cbPalette[3:8],lwd=4)
    grDevices::dev.off()
  }

  roc_p_0 = pROC::plot.roc(response,predictor,levels=c("NonSig","+/0"),smooth=F,col = cbPalette[3],main=paste0("ROC plot for ",method,", n=",samples))
  roc_p_n = pROC::plot.roc(response,predictor,levels=c("NonSig","+/-"),smooth=F,add=T,col = cbPalette[4])
  roc_0_p = pROC::plot.roc(response,predictor,levels=c("NonSig","0/+"),smooth=F,add=T,col = cbPalette[5])
  roc_0_n = pROC::plot.roc(response,predictor,levels=c("NonSig","0/-"),smooth=F,add=T,col = cbPalette[6])
  roc_n_p = pROC::plot.roc(response,predictor,levels=c("NonSig","-/+"),smooth=F,add=T,col = cbPalette[7])
  roc_n_0 = pROC::plot.roc(response,predictor,levels=c("NonSig","-/0"),smooth=F,add=T,col = cbPalette[8])
  graphics::legend(0.5,0.6,legend=c("+/0","+/-","0/+","0/-","-/+","-/0"),lty=1,col = cbPalette[3:8],lwd=4,y.intersp=0.55)
  return(list(roc_p_0=roc_p_0,
              roc_p_n=roc_p_n,
              roc_0_p=roc_0_p,
              roc_0_n=roc_0_n,
              roc_n_p=roc_n_p,
              roc_n_0=roc_n_0))
}
