#' @title Calculate modular differential connectivity (MDC)
#' @description Takes modules of genes (possibly overlapping) and calculates the change in correlation among those genes between two conditions. Also reports the genes with the strongest gain in connectivity (i.e., average difference in z-score of > 0) and the strongest loss of correlation between conditions for each module, if any pass the significance measure specified.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param compare Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.
#' @param genes A character vector specifying gene symbols, present as rows in the inputMat, corresponding to each module label in the labels argument.
#' @param labels A character vector specifying module label names, one for each gene symbol in the genes argument, with overlap allowed (i.e., each gene can be in more than one module).
#' @param corr_cutoff Cutoff specifying correlation values beyond which will be truncated to this value, to reduce the effect of outlier correlation values when using small sample sizes. Note that this does NOT affect the underlying correlation values, but does affect the z-score difference of correlation calculation in the dcTopPairs table. Default = 0.99
#' @param signType Coerce all correlation coefficients to be either positive (via "positive"), negative (via "negative"), or none (via "none") prior to calculating differential correlation. This could be used if, e.g., you think that going from a positive to a negative correlation is unlikely to occur biologically and is more likely to be due to noise, and you want to ignore these effects. Note that this does NOT affect the reported underlying correlation values, but does affect the z-score difference of correlation calculation. Default = "none", for no coercing.
#' @param oneSidedPVal If the dCorAvgType test is total_average, this option specifies whether a one-sided p-value should be reported, as opposed to a two-sided p-value. That is, if the average difference of z-scores is greater than zero, test whether the permutation average difference of z-scores are less than that average to get the p-value, and vice versa for the case that the average difference of z-scores is less than 0. Otherwise, test whether the absolute value of the average difference in z-scores is greater than the absolute values of the permutation average difference in z-scores. Default = FALSE.
#' @param dCorAvgMethod Character vector specifying the method for calculating the "average" differential correlation calculation that should be used. Options = "median", "mean".
#' @param corrType The correlation type of the analysis, limited to "pearson" or "spearman". Default = "pearson".
#' @param nPerms Number of permutations to generate in order to calculate the significance of the result.
#' @param gene_avg_signif The gene average differential correlation significance (adjusted for MHTC) required in order for the a gene to be reported as having a gain or loss in connectivity.
#' @param number_DC_genes The number of top differentially correlated genes with more correlation in each condition in each module to return in the data frame.
#' @param cl A parallel cluster object created by parallel::makeCluster(). If FALSE, defaults to single-core implementation. Default = FALSE.
#' @return A data frame with the module labels, the average change in difference in z-score between conditions (i.e., one measure of the modular average differential connectivity, or MeDC), and the empirical p-value for the significance of the change in correlation.
#' @examples
#' data(darmanis)
#' module_genes = list(mod1 = rownames(darmanis)[1:100],
#'  mod2 = rownames(darmanis)[90:190], mod3 = rownames(darmanis)[190:290])
#' modules = stack(module_genes)
#' modules$ind = as.character(modules$ind)
#' moduleDC_res = moduleDC(inputMat = darmanis, design = design_mat,
#'  compare = c("oligodendrocyte", "neuron"), genes = modules$values,
#'  labels = modules$ind)
#' @export
moduleDC <- function(inputMat, design, compare, genes, labels,
                     gene_avg_signif = 0.05, number_DC_genes = 3, adjust="bonferroni",
                     bySlope=FALSE){

  if(!length(genes) == length(labels)) stop("Genes and labels vectors must be the same length.")

  labels_names = unique(labels)

  mdc_vector = vector()
  mdc_signif = vector()
  mdc_signif_adj = vector()
  module_size = vector()
  goc_genes = vector()
  loc_genes = vector()

  for(i in 1:length(labels_names)){
    message(paste0("Calculating MDC for module #", i, ", which is called ", labels_names[i]))
      genes_tmp = genes[labels == labels_names[i]]
      module_size[i] = length(genes_tmp)
      inputMat_tmp = inputMat[match(genes_tmp,rownames(inputMat),nomatch=F), ]

      chow_res = chowCor(matA = inputMat_tmp, design_mat = design, compare = compare, corrType = corrType)
      supernova_res = flattenChow(chow_res)

      #Fisher's method of combining p-values is most appropriate
      # since the alternative is that SS_global > SS_subgroups
      # see: http://dx.doi.org/10.1093/biomet/asx076
      log2_pval_sum = sum(-2*log(supernova_res$pValDiff))
      n_pval = 2*sum(supernova_res$pValDiff)
      combined_p = 1-pchisq(log2_pval_sum,df=2*n_pval)
      adjust_p = p.adjust(combined_p,method=method,n=n_pval)

      mdc_vector[i] = log2_pval_sum
      mdc_signif[i] = combined_p
      mdc_signif_adj = adjust_p

      tmp = chow_res$pvalues
      tmp[lower.tri(tmp)] = t(tmp)[lower.tri(t(tmp))]
      pvalues_arr = tmp

      gene_level_p = apply(pvalues_arr,1,function(x){1-pchisq(sum(-2*log(x)),df=2*length(x))})
      #gene_avg = gene_level_p
      #gene_avg_sig = gene_avg[gene_avg < gene_avg_signif, ]

      group_num = match(compare,all_groups)

      if(bySlope){
        tmp = chow_res$slopes
        tmp[lower.tri(tmp)] = t(tmp)[lower.tri(t(tmp))]
        corrs_arr = tmp

        tmp = chow_res$globalSlope
        tmp[lower.tri(tmp)] = t(tmp)[lower.tri(t(tmp))]
        globalcorrs_arr = tmp
      } else{
        tmp = chow_res$corrs
        tmp[lower.tri(tmp)] = t(tmp)[lower.tri(t(tmp))]
        corrs_arr = tmp

        tmp = chow_res$globalCor
        tmp[lower.tri(tmp)] = t(tmp)[lower.tri(t(tmp))]
        globalcorrs_arr = tmp
      }

      if (length(compare)==1) {
        #compare one group against global model
        group1 = global_corrs_arr
        group2 = matrix(sapply(corrs_arr,function(x){as.numeric(unlist(strsplit(x,"/"))[group_num])},USE.NAMES = F),
                        nrow=nrow(corrs_arr),ncol=ncol(corrs_arr))
      } else if (length(compare)==2) {
        #compare two groups against each other
        group1 = matrix(sapply(corrs_arr,function(x){as.numeric(unlist(strsplit(x,"/"))[group_num[1]])},USE.NAMES = F),
                        nrow=nrow(corrs_arr),ncol=ncol(corrs_arr))
        group2 = matrix(sapply(corrs_arr,function(x){as.numeric(unlist(strsplit(x,"/"))[group_num[2]])},USE.NAMES = F),
                        nrow=nrow(corrs_arr),ncol=ncol(corrs_arr))
      } else {
        stop("'compare' may only be used with one or two subgroups.")
      }

      gene_avg_diff = apply(group2-group1,1,function(x){mean(x)})
      names(gene_avg_diff) = rownames(corrs_arr)
      gene_avg_goc = names(head(gene_avg_diff[order(gene_avg_diff,decreasing = T)],number_DC_genes))
      gene_avg_loc = names(head(gene_avg_diff[order(gene_avg_diff)],number_DC_genes))
      goc_genes[i] = paste(gene_avg_goc, collapse = ", ")
      loc_genes[i] = paste(gene_avg_loc, collapse = ", ")

  }

  res_df = data.frame(Module = labels_names, Size = module_size, MeDC = mdc_vector,
    pVal = mdc_signif, pValadj = mdc_signif_adj, gene_pVal = gene_level_p,
    Top_GOC = goc_genes, Top_LOC = loc_genes)

  return(res_df)

}
