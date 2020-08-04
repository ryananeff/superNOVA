#' @title Analyse superNOVA simulated data results for two group and three group simulated data.
#' @description Provides a way to merge expected and observed differential coexpressed tables and calculate the sensitivity/specificity/FPR/FNR.
#' @param ddcor_res The flattened chow_result file
#' @param true_pairs The true_pairs table from simulate_data() or simulate_data_3gr().
#' @param sig_cutoff The significance p-value (alpha) to use to call observed results as significant.
#' @param verbose Whether to print additional output or not.
#' @return Returns a list, including an expected versus observed table for the simulated results, and stats including TPR, TNR, FPR, FNR.
#' @keywords superNOVA simulate analyse
#' @export
analyse_sim_results <- function(ddcor_res,true_pairs,sig_cutoff=0.05,verbose=F){

  ddcor_res["is_signif"] = apply(ddcor_res,1,function(x){if(as.numeric(x["pValDiff_adj"])<=sig_cutoff){T}else{F}})
  ddcor_res[ddcor_res["is_signif"]==F,"Classes"] = "NonSig"

  ddcor_res_signif = ddcor_res[ddcor_res["is_signif"]==T,]
  true_pairs_signif = true_pairs[true_pairs["true_signif"]==TRUE,]

  ddcor_res_nonsignif = ddcor_res[(ddcor_res["is_signif"]==F),]
  true_pairs_nonsignif = true_pairs[true_pairs["true_signif"]==FALSE,]

  obs_pairs_signif = paste(ddcor_res_signif$Gene1, ddcor_res_signif$Gene2)
  act_pairs_signif = paste(true_pairs_signif$Gene1, true_pairs_signif$Gene2)
  obs_pairs_nonsignif = paste(ddcor_res_nonsignif$Gene1, ddcor_res_nonsignif$Gene2)
  act_pairs_nonsignif = paste(true_pairs_nonsignif$Gene1, true_pairs_nonsignif$Gene2)

  true_pos = length(intersect(obs_pairs_signif,act_pairs_signif))
  true_neg = length(intersect(obs_pairs_nonsignif,act_pairs_nonsignif))
  false_pos = length(intersect(obs_pairs_signif,act_pairs_nonsignif))
  false_neg = length(intersect(obs_pairs_nonsignif,act_pairs_signif))
  tpr = true_pos/length(act_pairs_signif) #sensitivity, recall
  tnr = true_neg/length(act_pairs_nonsignif) #specificity, selectivity
  fpr = false_pos/length(act_pairs_nonsignif) #type I error, fall-out
  fnr = false_neg/length(act_pairs_signif) #type 2 error, miss rate

  if(verbose){
    print(paste("true positive rate, sensitivity:", tpr))
    print(paste("true negative rate, specificity:", tnr))
    print(paste("false positive rate, Type I error:", fpr))
    print(paste("false negative rate, Type II error:", fnr))
  }

  rownames(ddcor_res) = paste0(ddcor_res$Gene1,"_",ddcor_res$Gene2)
  exp_vs_obs_table = cbind(true_pairs,ddcor_res[rownames(true_pairs),])
  return(list(exp_vs_obs_table=exp_vs_obs_table,
              tpr=tpr,
              tnr=tnr,
              fpr=fpr,
              fnr=fnr))

}
