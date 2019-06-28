#' @title Find groups of differentially coexpressed gene symbols.
#' @description Takes a table of differentially correlated genes with respect to one gene in the Gene2 column and returns the a list of vectors with unique, non-NA gene symbols for genes in each of the differentially correlated classes.
#' @param superNOVA_res The table of differential correlations outputted from superNOVA. Expected to have pValDiff or qValDiff columns as well as group_order, groupCor, globalCor, groupSlope, globalSlope, Gene1, Classes columns.
#' @param adjusted Logical indicating whether adjusted p-values from the differential correlation table (i.e., column "qValDiff", when adjusted = TRUE) or unadjusted p-values (i.e., column "pValDiff", when adjusted = FALSE) should be used to subset the table into significant and non-significant portions. Default = FALSE
#' @param pval_gene_thresh p-value threshold to call a gene as having significant differential correlation or not. Default = 0.05
#' @param classes Logical indicator specifying whether individual differential correlation gene classes should be extracted from the table or not. If not, compare MUST be set. Default = TRUE
#' @param compare Names of one or two groups to compare. If set to one group name, compare that group to the global model. If set to two group names, pairwise compare with each other. Default = NULL.
#' @param bySlope Whether to use the slope instead of the correlation rho when comparing one or two groups. This parameter has no effect when classes=TRUE.
#' @param unique_genes Logical, if TRUE indicates that unique gene symbols from each category compared to the other groups should be chosen prior to GO enrichment analysis.
#' @param geneNameCol Character vector specifying the name of the columns that are used to extract the gene symbols. Note that the default is c("Gene1", "Gene2"), but this only makes sense in the context of a full DGCA experiment. In the case of a splitSet, you may want to use "Gene1" to avoid counting the splitSet names in all of the categories.
#' @return A list of significantly differentially coexpressed genes.
#' @export
superNOVAFindSignificant <- function(superNOVA_res, pval_gene_thresh = 0.05, adjusted = FALSE,
                                     classes = TRUE, geneNameCol = c("Gene1", "Gene2"), unique_genes = FALSE,
                                     compare=NULL, bySlope=FALSE){

  #this gives the names of all groups compared
  all_groups = unlist(strsplit(superNOVA_res$group_order[1],"/"))

  if(adjusted){
    if(!"qValDiff" %in% colnames(superNOVA_res)){
      stop("If adjusted p-values are desired, then the input table must have an adjusted p-value column.")
    }
    superNOVA_res_sig = superNOVA_res[superNOVA_res$qValDiff < pval_gene_thresh, ]
  } else {
    superNOVA_res_sig = superNOVA_res[superNOVA_res$pValDiff < pval_gene_thresh, ]
  }

  if(classes){ #break down everything by classes
    superNOVA_res_sig_list = split(superNOVA_res,superNOVA_res$Classes) #this creates a list of data frames
    if(unique_genes){
      for(i in 1:length(superNOVA_res_sig_list)){
        superNOVA_res_sig_list_i_rem = superNOVA_res_sig_list[-i]
        superNOVA_res_sig_list[[i]] = setdiff(superNOVA_res_sig_list[[i]], superNOVA_res_sig_list_i_rem)
      }
    }
  } else if(!(is.null(compare))){ #if there are specific groups to compare
    if(anyNA(match(compare,all_groups))){
      #sanity check - ensure input groups to compare match result groups
      stop("input group(s) in 'compare' do not match result groups in 'superNOVA_res$group_order'")
    }
    group_num = match(compare,all_groups)
    if (length(compare)==1) {
      #compare one group against global model
      if(bySlope){
        group1 = superNOVA_res_sig$globalSlope
        group2 = sapply(superNOVA_res_sig$groupSlope,
                        FUN=function(x){as.numeric(unlist(strsplit(x,"/"))[group_num])},
                        USE.NAMES = F)
      } else {
        group1 = superNOVA_res_sig$globalCor
        group2 = sapply(superNOVA_res_sig$groupCor,
                        FUN=function(x){as.numeric(unlist(strsplit(x,"/"))[group_num])},
                        USE.NAMES = F)
      }
    } else if (length(compare)==2) {
      #compare two groups against each other
      if(bySlope){
        group1 = sapply(superNOVA_res_sig$groupSlope,
                        FUN=function(x){as.numeric(unlist(strsplit(x,"/"))[group_num[1]])},
                        USE.NAMES = F)
        group2 = sapply(superNOVA_res_sig$groupSlope,
                        FUN=function(x){as.numeric(unlist(strsplit(x,"/"))[group_num[2]])},
                        USE.NAMES = F)
      } else {
        group1 = sapply(superNOVA_res_sig$groupCor,
                        FUN=function(x){as.numeric(unlist(strsplit(x,"/"))[group_num[1]])},
                        USE.NAMES = F)
        group2 = sapply(superNOVA_res_sig$groupCor,
                        FUN=function(x){as.numeric(unlist(strsplit(x,"/"))[group_num[2]])},
                        USE.NAMES = F)
      }
    } else {
      stop("'compare' may only be used with one or two subgroups.")
    }
    superNOVA_res_sig_pos = unique(unlist(superNOVA_res_sig[group2>group1,geneNameCol]))
    superNOVA_res_sig_neg = unique(unlist(superNOVA_res_sig[group2<group1,geneNameCol]))

    if(unique_genes){
      superNOVA_res_sig_pos_copy = superNOVA_res_sig_pos
      superNOVA_res_sig_neg_copy = superNOVA_res_sig_neg
      superNOVA_res_sig_pos = setdiff(superNOVA_res_sig_pos, superNOVA_res_sig_neg_copy)
      superNOVA_res_sig_neg = setdiff(superNOVA_res_sig_neg, superNOVA_res_sig_pos_copy)
    }
    superNOVA_res_sig_list = list(
      significant_gain_of_correlation_genes = superNOVA_res_sig_pos,
      significant_loss_of_correlation_genes = superNOVA_res_sig_neg)
  } else {
    stop("either 'classes' must be TRUE or 'compare' is not NULL.")
  }
  #filter all of the lists to remove NAs and select only the unique genes
  superNOVA_res_sig_list = lapply(superNOVA_res_sig_list, function(x) x[!is.na(x)])
  superNOVA_res_sig_list = lapply(superNOVA_res_sig_list, function(x) unique(x))
  return(superNOVA_res_sig_list)
}
