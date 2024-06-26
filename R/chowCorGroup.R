#' @title Runs the superNOVA pipeline directly on one or two feature matrices.
#' @description Evaluates differential coexpression between two or more subgroups of samples in the data versus the global model.
#' @param design_mat A model.matrix type design matrix, denoting which samples are in which subgroups of the dataset. Required.
#' @param matA A features (rows) by samples (columns) data.frame of feature values, such as a gene expression matrix. Required.
#' @param matB A second features (rows) by samples (columns) data.frame of feature values, such as a gene expression matrix. Values in this matrix will be compared to matA if provided. Default=NULL.
#' @param compare Which groups of samples in the design matrix should we evaluate? Otherwise, all groups in the design matrix will be evaluated. Default=NULL.
#' @param corrType The base correlation metric to evaluate coexpression. One of "pearson","spearman", or "bicor". Default=pearson.
#' @return Returns a list of matrices, includes pvalues for the superNOVA test (pvalues) and their classes (classes), correlation values for each subgroup model (corrs) and their pvalues (corrsP), correlation values for the global model (globalCor) and its pvalues (globalCorP), and a flag to indicate whether a second matrix was used (secondMat).
#' @keywords superNOVA
#' @importFrom stats cov na.omit pf pt var
#' @importFrom progress progress_bar
#' @export
chowCorGroup = function(design_mat,matA,matB=NULL,compare=NULL,corrType="pearson"){
  
  if(!is.null(matB)){secondMat=TRUE} else {secondMat=FALSE; matB = matA}
  if(!is.null(compare)){
    if(length(compare)<2){stop("number of groups to compare must be 2 or larger")}
    if(length(intersect(compare,colnames(design_mat)))<length(compare)){stop("group(s) to compare not in design matrix")}
    design_mat = design_mat[,compare]
  }
  design_mat = design_mat[rowSums(design_mat)>0,] #drop rows no longer relevant
  tryCatch({
    matA = matA[,rownames(design_mat),drop=FALSE]
    matB = matB[,rownames(design_mat),drop=FALSE]
  },error=function(e){
    stop("Row names in design matrix do not match column names in input.")
  })
  
  
  if(secondMat){ #if we are comparing two matrices
    
    results_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(results_mat) = rownames(matA) #gene1 names
    colnames(results_mat) = rownames(matB) #gene2 names
    
    ftest_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(ftest_mat) = rownames(matA) #gene1 names
    colnames(ftest_mat) = rownames(matB) #gene2 names
    
    corrs_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(corrs_mat) = rownames(matA) #gene1 names
    colnames(corrs_mat) = rownames(matB) #gene2 names
    
    slopes_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(slopes_mat) = rownames(matA) #gene1 names
    colnames(slopes_mat) = rownames(matB) #gene2 names
    
    pvals_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(pvals_mat) = rownames(matA) #gene1 names
    colnames(pvals_mat) = rownames(matB) #gene2 names
    
    global_corrs_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(global_corrs_mat) = rownames(matA) #gene1 names
    colnames(global_corrs_mat) = rownames(matB) #gene2 names
    
    global_slopes_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(global_slopes_mat) = rownames(matA) #gene1 names
    colnames(global_slopes_mat) = rownames(matB) #gene2 names
    
    global_pvals_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(global_pvals_mat) = rownames(matA) #gene1 names
    colnames(global_pvals_mat) = rownames(matB) #gene2 names
    
    classes_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(classes_mat) = rownames(matA) #gene1 names
    colnames(classes_mat) = rownames(matB) #gene2 names
    
    varsA_group = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matA))
    varsB_group = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
    covs_group = array(NA,c(ncol(design_mat),nrow(matA),nrow(matB)))
    
    
    for(groupn in 1:ncol(design_mat)){
      matAg = matA[,design_mat[,groupn]==1,drop=FALSE]
      matBg = matB[,design_mat[,groupn]==1,drop=FALSE]
      
      #normalize for expression differences between groups
      matA[,design_mat[,groupn]==1] = t(scale(t(matAg), scale = FALSE))
      matB[,design_mat[,groupn]==1] = t(scale(t(matBg), scale = FALSE))
      
      varsA_group[groupn,] = apply(matAg,1,function(x){var(x)}) #variance for genes in A
      varsB_group[groupn,] = apply(matBg,1,function(x){var(x)}) #variance for genes in B
    }
    
    varsA = apply(matA,1,function(x){var(x)})
    varsB = apply(matB,1,function(x){var(x)})
    
    pb <- progress::progress_bar$new(
      format = " [chowCor] [:bar] :current/:total (:percent) eta: :eta gene: :what",
      clear = FALSE, total = nrow(matA), width = 100)
    
    for(ix_row in 1:nrow(matA)){ #go row-by-row in matA and compare with matB, then reassemble results matrix
      
      pb$tick(tokens = list(what = sprintf(rownames(matA)[ix_row],fmt="%15s")))
      
      sse_groups = rep(0,nrow(matB))
      corrs_rg = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
      slopes_rg = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
      pvals_rg = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
      classes_rg = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
      for(groupn in 1:ncol(design_mat)){
        #calculate stats for row
        matAig = matA[ix_row,design_mat[,groupn]==1,drop=FALSE] #get row for matA
        matBg = matB[,design_mat[,groupn]==1,drop=FALSE]
        
        varA = varsA_group[groupn,ix_row,drop=FALSE]
        varB = varsB_group[groupn,,drop=FALSE]
        covsAB = cov(t(matAig),t(matBg))
        covs_group[groupn,ix_row,] = covsAB
        
        if(corrType=="bicor"){
          corrs = WGCNA::bicor(t(matAig),t(matBg),use="pairwise.complete.obs")
        } else {
          corrs = WGCNA::cor(t(matAig),t(matBg),,use="pairwise.complete.obs",method=corrType) #get Pearson rho, using Cpp fxn from WGCNA
        }
        
        slopes = corrs * (varB / varA[1])
        
        #calculate p-value for the individual groups
        deg_freedom = ncol(matAig)-2
        in_pt = (abs(corrs) * sqrt(deg_freedom) / sqrt(1 - corrs^2))
        pvals = 2 * (1 - pt(in_pt, deg_freedom))
        corrs_rg[groupn,] = corrs
        pvals_rg[groupn,] = pvals
        corrs[is.na(corrs)] = 0 #if NA, ignore
        pvals[is.na(pvals)] = 1 #if NA, ignore
        slopes_rg[groupn,] = slopes
        classes_rg[groupn,] = sapply(1:length(corrs),function(x){if((pvals[x]<0.05)&(corrs[x]>0)){"+"}else if((pvals[x]<0.05)&(corrs[x]<0)){"-"}else{"0"}})
        
        #calculate the mlr for each row
        #sse_groups = sse_groups + t(sapply(1:nrow(matBg),function(x){ (varB[x]-covsAB[x]*(corrs[x]/sqrt(varA/varB[x])))*(ncol(matBg)-1)}))
        sse_groups = sse_groups + (varB-covsAB*(corrs/sqrt(varA%*%(1/varB))))*(ncol(matBg)-1)
      } #end group
      matAi = matA[ix_row,,drop=FALSE]
      varA = varsA[ix_row]
      varB = varsB
      covsAB = cov(t(matAi),t(matB))
      if(corrType=="bicor"){
        corrs = WGCNA::bicor(t(matAi),t(matB),,use="pairwise.complete.obs")
      } else {
        corrs = WGCNA::cor(t(matAi),t(matB),,use="pairwise.complete.obs",method=corrType)
      }
      
      slopes = corrs * (varB / varA[1])
      
      #corr test
      deg_freedom = ncol(matAi)-2
      in_pt = (abs(corrs) * sqrt(deg_freedom) / sqrt(1 - corrs^2))
      corr_pvals = 2 * (1 - pt(in_pt, deg_freedom))
      
      #chow test
      #Chow test
      sse_all = array(0,nrow(matB))
      num_compare = 0
      num_groups = ncol(design_mat)
      for(group_A in 1:ncol(design_mat)){
        n_A = sum(design_mat[,group_A])
        var_xA = varsA_group[group_A,ix_row,drop=FALSE][1]
        var_yA = array(varsB_group[group_A,,drop=FALSE])
        cov_A = array(covs_group[group_A,ix_row,,drop=FALSE])
        if ((group_A+1)<=ncol(design_mat)){
          for(group_B in (group_A+1):ncol(design_mat)){
            n_B = sum(design_mat[,group_B])
            var_xB = varsA_group[group_B,ix_row,drop=FALSE][1]
            var_yB = array(varsB_group[group_B,,drop=FALSE])
            cov_B = array(covs_group[group_B,ix_row,,drop=FALSE])
            sse_ab = n_A*(var_yA-2*cov_B/var_xB*cov_A + cov_B^2/var_xB^2 * var_xA)
            sse_ba = n_B*(var_yB-2*cov_A/var_xA * cov_B + cov_A^2/var_xA^2 * var_xB)
            sse_all = sse_all + (sse_ab+sse_ba)
            num_compare = num_compare + 2
          }
        }
      }
      
      
      #print(sse_all)
      k = 2 #this is equivalent to the number of params for a linear relationship (y=mx+b)
      ngr <- ncol(design_mat) #number of groups
      nsamp <- nrow(design_mat) #number of samples
      sse_all = t(sse_all) / (num_compare/num_groups)
      f <- ((sse_all-sse_groups)/((k+1)))/(sse_groups/(nsamp-ngr*(k+1)))
      p <- pf(f/2, ngr/num_compare*k, nsamp-ngr*k, lower.tail=FALSE)
      p[is.na(p)] = 1
      
      corrs_mat[ix_row,] = apply(t(corrs_rg),1,paste,collapse="/")
      slopes_mat[ix_row,] = apply(t(slopes_rg),1,paste,collapse="/")
      pvals_mat[ix_row,] = apply(t(pvals_rg),1,paste,collapse="/")
      global_corrs_mat[ix_row,] = corrs
      global_pvals_mat[ix_row,] = corr_pvals
      global_slopes_mat[ix_row,] = slopes
      classes_mat[ix_row,] = apply(t(classes_rg),1,paste,collapse="/")
      ftest_mat[ix_row,] = f
      results_mat[ix_row,] = p
    } #end row
  } else { #if we are only looking at one input matrix (comparing all gene pairs)
    matB = matA
    
    results_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(results_mat) = rownames(matA) #gene1 names
    colnames(results_mat) = rownames(matB) #gene2 names
    
    ftest_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(ftest_mat) = rownames(matA) #gene1 names
    colnames(ftest_mat) = rownames(matB) #gene2 names
    
    corrs_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(corrs_mat) = rownames(matA) #gene1 names
    colnames(corrs_mat) = rownames(matB) #gene2 names
    
    pvals_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(pvals_mat) = rownames(matA) #gene1 names
    colnames(pvals_mat) = rownames(matB) #gene2 names
    
    slopes_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(slopes_mat) = rownames(matA) #gene1 names
    colnames(slopes_mat) = rownames(matB) #gene2 names
    
    global_corrs_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(global_corrs_mat) = rownames(matA) #gene1 names
    colnames(global_corrs_mat) = rownames(matB) #gene2 names
    
    global_slopes_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(global_slopes_mat) = rownames(matA) #gene1 names
    colnames(global_slopes_mat) = rownames(matB) #gene2 names
    
    global_pvals_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(global_pvals_mat) = rownames(matA) #gene1 names
    colnames(global_pvals_mat) = rownames(matB) #gene2 names
    
    classes_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(classes_mat) = rownames(matA) #gene1 names
    colnames(classes_mat) = rownames(matB) #gene2 names
    
    varsA_group = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matA))
    covs_group = array(NA,c(ncol(design_mat),nrow(matA),nrow(matA)))
    varsB_group = varsA_group
    
    for(groupn in 1:ncol(design_mat)){ #for each group in the design matrix
      # isolate those samples from the design matrix to calculate
      # the correlations, variances, slopes, etc. separately
      matAg = matA[,design_mat[,groupn]==1,drop=FALSE]
      matBg = matB[,design_mat[,groupn]==1,drop=FALSE]
      
      #normalize for expression differences between groups
      matA[,design_mat[,groupn]==1] = t(scale(t(matAg), scale = FALSE))
      matB[,design_mat[,groupn]==1] = matA[,design_mat[,groupn]==1]
      
      varsA_group[groupn,] = apply(matAg,1,function(x){var(x)}) #variance for genes in A
      varsB_group[groupn,] = varsA_group[groupn,] #variance for genes in B
    }
    
    varsA = apply(matA,1,function(x){var(x)})
    varsB = varsA
    
    pb <- progress::progress_bar$new(
      format = " [chowCor] [:bar] :current/:total (:percent) eta: :eta gene: :what",
      clear = FALSE, total = nrow(matA), width = 100)
    
    for(ix_row in 1:nrow(matA)){ #go row-by-row in matA and compare with matB, then reassemble results matrix
      
      pb$tick(tokens = list(what = sprintf(rownames(matA)[ix_row],fmt="%15s")))
      
      sse_groups = rep(0,nrow(matB)-ix_row+1)
      corrs_rg = matrix(NA,nrow=ncol(design_mat),ncol=(nrow(matB)-ix_row+1))
      slopes_rg = matrix(NA,nrow=ncol(design_mat),ncol=(nrow(matB)-ix_row+1))
      pvals_rg = matrix(NA,nrow=ncol(design_mat),ncol=(nrow(matB)-ix_row+1))
      classes_rg = matrix(NA,nrow=ncol(design_mat),ncol=(nrow(matB)-ix_row+1))
      
      for(groupn in 1:ncol(design_mat)){
        #calculate stats for row
        matAig = matA[ix_row,design_mat[,groupn]==1,drop=FALSE] #get row for matA
        matBg = matB[ix_row:nrow(matB),design_mat[,groupn]==1,drop=FALSE]
        
        varA = varsA_group[groupn,ix_row,drop=FALSE]
        varB = varsB_group[groupn,ix_row:nrow(matB),drop=FALSE]
        covsAB = cov(t(matAig),t(matBg))
        covs_group[groupn,ix_row,ix_row:nrow(matB)] = covsAB
        
        if(corrType=="bicor"){
          corrs = WGCNA::bicor(t(matAig),t(matBg),,use="pairwise.complete.obs")
        } else {
          corrs = WGCNA::cor(t(matAig),t(matBg),,use="pairwise.complete.obs",method=corrType) #get Pearson rho, using Cpp fxn from WGCNA
        }
        
        slopes = corrs * (varB / varA[1])
        
        #calculate p-value for the individual groups
        deg_freedom = ncol(matAig)-2
        in_pt = (abs(corrs) * sqrt(deg_freedom) / sqrt(1 - corrs^2))
        pvals = 2 * (1 - pt(in_pt, deg_freedom))
        corrs_rg[groupn,] = corrs
        pvals_rg[groupn,] = pvals
        slopes_rg[groupn,] = slopes
        corrs[is.na(corrs)] = 0 #if NA, ignore
        pvals[is.na(pvals)] = 1 #if NA, ignore
        classes_rg[groupn,] = sapply(1:length(corrs),
                                     function(x){if(
                                       (pvals[x]<0.05)&(corrs[x]>0)){"+"} else if ((pvals[x]<0.05)&(corrs[x]<0)) {"-"} else{"0"}})
        sse_groups = sse_groups + (varB-covsAB*(corrs/sqrt(varA%*%(1/varB))))*(ncol(matBg)-1)
      } #end group
      matAi = matA[ix_row,,drop=FALSE] #get row for matA
      matBi = matB[ix_row:nrow(matB),,drop=FALSE]
      
      varA = varsA[ix_row]
      varB = varsB[ix_row:nrow(matB)]
      covsAB = cov(t(matAi),t(matBi))
      if(corrType=="bicor"){
        corrs = WGCNA::bicor(t(matAi),t(matBi),,use="pairwise.complete.obs")
      } else {
        corrs = WGCNA::cor(t(matAi),t(matBi),,use="pairwise.complete.obs",method=corrType)
      }
      
      slopes = corrs * (varB / varA[1])
      
      #corr test
      deg_freedom = ncol(matAi)-2
      in_pt = (abs(corrs) * sqrt(deg_freedom) / sqrt(1 - corrs^2))
      corr_pvals = 2 * (1 - pt(in_pt, deg_freedom))
      
      #Chow test
      sse_all = array(0,c(length(ix_row:nrow(matB))))
      num_compare = 0
      num_groups = ncol(design_mat)
      for(group_A in 1:ncol(design_mat)){
        n_A = sum(design_mat[,group_A])
        var_xA = varsA_group[group_A,ix_row,drop=FALSE][1]
        var_yA = array(varsB_group[group_A,ix_row:nrow(matB),drop=FALSE])
        cov_A = array(covs_group[group_A,ix_row,ix_row:nrow(matB),drop=FALSE])
        if ((group_A+1)<=ncol(design_mat)){
          for(group_B in (group_A+1):ncol(design_mat)){
            n_B = sum(design_mat[,group_B])
            var_xB = varsA_group[group_B,ix_row,drop=FALSE][1]
            var_yB = array(varsB_group[group_B,ix_row:nrow(matB),drop=FALSE])
            cov_B = array(covs_group[group_B,ix_row,ix_row:nrow(matB),drop=FALSE])
            sse_ab = n_A*(var_yA-2*cov_B/var_xB*cov_A + cov_B^2/var_xB^2 * var_xA)
            sse_ba = n_B*(var_yB-2*cov_A/var_xA * cov_B + cov_A^2/var_xA^2 * var_xB)
            sse_all = sse_all + (sse_ab+sse_ba)
            num_compare = num_compare + 2
          }
        }
      }
      
      sse_all = t(sse_all) / (num_compare/num_groups)
      #print(sse_all)
      k = 2 #this is equivalent to the number of params for a linear relationship (y=mx+b)
      ngr <- ncol(design_mat) #number of groups
      nsamp <- nrow(design_mat) #number of samples
      f <- ((sse_all-sse_groups)/((k+1)))/(sse_groups/(nsamp-ngr*(k+1)))
      p <- pf(f/2, ngr/num_compare*k, nsamp-ngr*k, lower.tail=FALSE)
      p[is.na(p)] = 1
      
      corrs_mat[ix_row,ix_row:nrow(matB)] = apply(t(corrs_rg),1,paste,collapse="/")
      slopes_mat[ix_row,ix_row:nrow(matB)] = apply(t(slopes_rg),1,paste,collapse="/")
      pvals_mat[ix_row,ix_row:nrow(matB)] = apply(t(pvals_rg),1,paste,collapse="/")
      global_corrs_mat[ix_row,ix_row:nrow(matB)] = corrs
      global_pvals_mat[ix_row,ix_row:nrow(matB)] = corr_pvals
      global_slopes_mat[ix_row,ix_row:nrow(matB)] = slopes
      classes_mat[ix_row,ix_row:nrow(matB)] = apply(t(classes_rg),1,paste,collapse="/")
      ftest_mat[ix_row,ix_row:nrow(matB)] = f
      results_mat[ix_row,ix_row:nrow(matB)] = p
      #list(classes=apply(t(classes_rg),1,paste,collapse="/"),pvalues=p)
    } #end row
  } #end else
  diag(results_mat) = 1
  output = list(corrs=corrs_mat, slopes=slopes_mat, corrsP=pvals_mat, globalCor=global_corrs_mat,
                globalSlope=global_slopes_mat, globalCorP=global_pvals_mat, Ftest = ftest_mat,
                pvalues=results_mat,classes=classes_mat, secondMat=secondMat,
                groups_compared=paste(colnames(design_mat),collapse="/"))
  return(output)
}
