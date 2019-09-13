#' @title Calculate modular differential connectivity (MDC)
#' @description Takes modules of genes (possibly overlapping) and calculates the change in correlation among those genes between two conditions. Also reports the genes with the strongest gain in connectivity (i.e., average difference in z-score of > 0) and the strongest loss of correlation between conditions for each module, if any pass the significance measure specified.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param inputMatB The second matrix (or data.frame) to differentially correlate with the first of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param compare Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.
#' @param genes A character vector specifying gene symbols, present as rows in the inputMat, corresponding to each module label in the labels argument.
#' @param labels A character vector specifying module label names, one for each gene symbol in the genes argument, with overlap allowed (i.e., each gene can be in more than one module).
#' @param gene_avg_signif Minimum average gene level significance to report in the module.
#' @param number_DC_genes The number of differentially GOC or LOC genes to output in the final result.
#' @param adjust The method to adjust p-values (adjusts for multiple modules). Same inputs as p.adjust().
#' @param bySlope Whether or not to use the slopes instead of the correlations when calling GOC or LOC in the output for the two group case.
#' @param corrType The correlation type out of Pearson, Spearman, and bicor for the underlying linear model.
#' @param parallel Run the pipeline as single core or multithreaded? True or False.
#' @param batchConfig When running in parallel, the path to the batchtools configuration file. Best to run locally.
#' @param split When running in parallel, the number of times to split the input.
#' @return A data frame with the module labels, the average change in difference in z-score between conditions (i.e., one measure of the modular average differential connectivity, or MeDC), and the empirical p-value for the significance of the change in correlation.
#' @examples
#' data(design_mat)
#' data(darmanis)
#' rownames(design_mat) = colnames(darmanis)
#' module_genes = list(mod1 = rownames(darmanis)[1:100],
#'  mod2 = rownames(darmanis)[90:190], mod3 = rownames(darmanis)[190:290])
#' modules = stack(module_genes)
#' modules$ind = as.character(modules$ind)
#' moduleDC_res = moduleDC(inputMat = darmanis, design = design_mat,
#'  compare = c("oligodendrocyte", "neuron"), genes = modules$values,
#'  labels = modules$ind)
#' @export
moduleDC <- function(inputMat=inputMat, inputMatB=NULL, design=design, compare=compare, genes=genes, labels=labels,
                     gene_avg_signif = 0.05, number_DC_genes = 3, adjust="bonferroni",
                     bySlope=FALSE,corrType="bicor",parallel=F,
                     batchConfig = system.file("config/batchConfig_Local.R",package="superNOVA"),
                     split){

  if(!length(genes) == length(labels)) stop("Genes and labels vectors must be the same length.")

  labels_names = unique(labels)

  mdc_vector = vector()
  mdc_signif = vector()
  mdc_signif_adj = vector()
  module_size = vector()
  df = vector()
  genes_sig = vector()
  goc_genes = vector()
  loc_genes = vector()
  gene_pval = vector()

  for(i in 1:length(labels_names)){
    message(paste0("[moduleDC] Calculating MDC for module #", i, " of ",length(labels_names),", which is called ", labels_names[i]))
      genes_tmp = genes[labels == labels_names[i]]
      module_size[i] = length(genes_tmp)
      inputMat_tmp = inputMat[match(genes_tmp,rownames(inputMat),nomatch=F), ,drop=FALSE]
      if(nrow(inputMat_tmp) <= 1){
        message('[moduleDC] too few genes in set, skipping...')
        next
      }
      if(!(parallel)){
        chow_res = chowCor(matA = inputMat_tmp, matB=inputMatB, design_mat = design, compare = compare, corrType = corrType)
        supernova_res = flattenChow(chow_res, adjust_q=F,sort_output=F)
      } else {
        outputfile = "tmp_moduleDC.tsv"
        combined = rbind(inputMat_tmp,inputMatB) #temporarily for now
        superNOVA::chowParallel(design=design_mat,
                              inputMat = combined,
                              compare = compare,
                              corrType = corrType, outputFile = outputfile,
                              batchConfig = batchConfig,
                              perBatch=split,coresPerJob = 1,batchDir = "tmp_batch/",sigOutput=F)
        supernova_res = read.table(outputfile,header=T)
      }

      #Fisher's method of combining p-values is most appropriate
      # since the alternative is that SS_global > SS_subgroups
      # see: http://dx.doi.org/10.1093/biomet/asx076
      log2_pval_sum = sum(-2*log(supernova_res$pValDiff))
      n_pval = 2*length(supernova_res$pValDiff)
      combined_p = 1-pchisq(log2_pval_sum,df=n_pval)
      adjust_p = p.adjust(combined_p,method=adjust,n=length(labels_names))

      df[i] = n_pval
      mdc_vector[i] = log2_pval_sum
      mdc_signif[i] = combined_p
      mdc_signif_adj[i] = adjust_p

      message("[moduleDC] Calculating gene-level stats")

      tmp = chow_res$pvalues
      tmp[lower.tri(tmp)] = t(tmp)[lower.tri(t(tmp))]
      pvalues_arr = tmp

      gene_level_p = apply(pvalues_arr,1,function(x){1-pchisq(sum(-2*log(x)),df=2*length(x))})
      gene_level_p = gene_level_p[order(gene_level_p)]
      #gene_avg = gene_level_p
      genes_sig[i] = sum(gene_level_p < gene_avg_signif)

      #all_groups = colnames(design)

      #group_num = match(compare,all_groups)

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

      if (length(compare)==2) { #TODO: how to represent GOC and LOC for 3 or more groups?
        #compare two groups against each other
        group1 = matrix(sapply(corrs_arr,function(x){as.numeric(unlist(strsplit(x,"/"))[1])},USE.NAMES = F),
                        nrow=nrow(corrs_arr),ncol=ncol(corrs_arr))
        group2 = matrix(sapply(corrs_arr,function(x){as.numeric(unlist(strsplit(x,"/"))[2])},USE.NAMES = F),
                        nrow=nrow(corrs_arr),ncol=ncol(corrs_arr))
        gene_avg_diff = apply(group2-group1,1,function(x){mean(x)})
        names(gene_avg_diff) = rownames(corrs_arr)
        gene_avg_goc = head(gene_avg_diff[order(gene_avg_diff,decreasing = T)],number_DC_genes)
        gene_avg_loc = head(gene_avg_diff[order(gene_avg_diff)],number_DC_genes)
        gene_pval[i] = paste(lapply(seq_along(gene_level_p),
                                    function(y, n, i) { paste(n[[i]], round(y[[i]],4),sep=":") },
                                    y=gene_level_p, n=names(gene_level_p)),collapse="; ")
        goc_genes[i] = paste(lapply(seq_along(gene_avg_goc),
                                    function(y, n, i) { paste(n[[i]], round(y[[i]],4),sep=":") },
                                    y=gene_avg_goc, n=names(gene_avg_goc)),collapse="; ")
        loc_genes[i] = paste(lapply(seq_along(gene_avg_loc),
                                    function(y, n, i) { paste(n[[i]], round(y[[i]],4),sep=":") },
                                    y=gene_avg_loc, n=names(gene_avg_loc)),collapse="; ")
      }
  }

  res_df = data.frame(Module = labels_names, Size = module_size, genesSignif=genes_sig, MeDC = mdc_vector, df=df,
    pVal = mdc_signif, pValadj = mdc_signif_adj, gene_pVal = gene_pval,
    Top_GOC = goc_genes, Top_LOC = loc_genes)

  return(res_df)

}
