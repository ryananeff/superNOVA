flattenChow = function(chow_result){
  rownames = rownames(chow_result$corrs)
  colnames = colnames(chow_result$corrs)
  corrsflat = cbind(which(!is.na(chow_result$corrs),arr.ind = TRUE),na.omit(as.vector(chow_result$corrs)))
  corrsflat[,1] = sapply(corrsflat[,1],function(x){rownames[as.numeric(x)]})
  corrsflat[,2] = sapply(corrsflat[,2],function(x){colnames[as.numeric(x)]})
  columns = c("GeneA","GeneB","groupCor","groupCorPval","pValDiff","Classes")
  output = data.frame(matrix(NA,nrow=nrow(corrsflat),ncol=length(columns)))
  output[,1:3] = corrsflat
  output[,4] = na.omit(as.vector(chow_result$corrsP))
  output[,5] = as.numeric(na.omit(as.vector(chow_result$pvalues)))
  output[,6] = na.omit(as.vector(chow_result$classes))
  colnames(output) = columns
  output = output[output$GeneA!=output$GeneB,] #don't output self matching rows
  output = output[order(output$pValDiff,method = "radix"),] #order the rows with most significant first
  return(output)
}

chowCor = function(design_mat,matA,matB=NULL,compare=NULL,corrType="pearson"){
  if(!is.null(matB)){secondMat=TRUE} else {secondMat=FALSE}
  if(!is.null(compare)){
    if(length(compare)<2){stop("number of groups to compare must be 2 or larger")}
    if(length(intersect(compare,colnames(design_mat)))<length(compare)){stop("group(s) to compare not in design matrix")}
    design_mat = design_mat[,compare]
    design_mat = design_mat[rowSums(design_mat)>0,] #drop rows no longer relevant
    matA = matA[,rownames(design_mat)]
    matB = matB[,rownames(design_mat)]
  }
  if(secondMat){ #if we are comparing two matrices
    
    results_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(results_mat) = rownames(matA) #gene1 names
    colnames(results_mat) = rownames(matB) #gene2 names
    
    corrs_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(corrs_mat) = rownames(matA) #gene1 names
    colnames(corrs_mat) = rownames(matB) #gene2 names
    
    pvals_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(pvals_mat) = rownames(matA) #gene1 names
    colnames(pvals_mat) = rownames(matB) #gene2 names
    
    classes_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(classes_mat) = rownames(matA) #gene1 names
    colnames(classes_mat) = rownames(matB) #gene2 names
    
    varsA_group = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matA))
    varsB_group = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
    
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
    
    for(ix_row in 1:nrow(matA)){ #go row-by-row in matA and compare with matB, then reassemble results matrix
      sse_groups = rep(0,nrow(matB))
      corrs_rg = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
      pvals_rg = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
      classes_rg = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matB))
      for(groupn in 1:ncol(design_mat)){
        #calculate stats for row
        matAig = matA[ix_row,design_mat[,groupn]==1,drop=FALSE] #get row for matA
        matBg = matB[,design_mat[,groupn]==1,drop=FALSE]
        
        varA = varsA_group[groupn,ix_row,drop=FALSE]
        varB = varsB_group[groupn,,drop=FALSE]
        covsAB = cov(t(matAig),t(matBg))
        
        if(corrType=="bicor"){
        corrs = WGCNA::bicor(t(matAig),t(matBg),use="all.obs")
        } else {
        corrs = WGCNA::cor(t(matAig),t(matBg),use="all.obs",method=corrType) #get Pearson rho, using Cpp fxn from WGCNA
        }
        
        #calculate p-value for the individual groups
        deg_freedom = ncol(matAig)-2
        in_pt = (abs(corrs) * sqrt(deg_freedom) / sqrt(1 - corrs^2))
        pvals = 2 * (1 - pt(in_pt, deg_freedom))
        corrs_rg[groupn,] = corrs
        pvals_rg[groupn,] = pvals
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
      corrs = WGCNA::bicor(t(matAi),t(matB),use="all.obs")    
      } else {
      corrs = WGCNA::cor(t(matAi),t(matB),use="all.obs",method=corrType)
      }
      #sse_all = t(sapply(1:nrow(matB),function(x){ (varB[x]-covsAB[x]*(corrs[x]/sqrt(varA/varB[x])))*(ncol(matB)-1)}))
      sse_all = (varB-covsAB*(corrs/sqrt(varA%*%(1/varB))))*(ncol(matB)-1)
      
      df1 <- ncol(design_mat)
      df2 <- nrow(design_mat)-2*(df1)
      
      f <- (sse_all-sse_groups)*df2/(df1*sse_groups)
      p <- pf(f, df1, df2, lower.tail=FALSE)
      
      corrs_mat[ix_row,] = apply(t(corrs_rg),1,paste,collapse="/")
      pvals_mat[ix_row,] = apply(t(pvals_rg),1,paste,collapse="/")
      classes_mat[ix_row,] = apply(t(classes_rg),1,paste,collapse="/")
      results_mat[ix_row,] = p
    } #end row
  } else { #if we are only looking at one input matrix (comparing all gene pairs)
    matB = matA
    
    results_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(results_mat) = rownames(matA) #gene1 names
    colnames(results_mat) = rownames(matB) #gene2 names
    
    corrs_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(corrs_mat) = rownames(matA) #gene1 names
    colnames(corrs_mat) = rownames(matB) #gene2 names
    
    pvals_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(pvals_mat) = rownames(matA) #gene1 names
    colnames(pvals_mat) = rownames(matB) #gene2 names
    
    classes_mat = matrix(NA,nrow=nrow(matA),ncol=nrow(matB))
    rownames(classes_mat) = rownames(matA) #gene1 names
    colnames(classes_mat) = rownames(matB) #gene2 names
    
    varsA_group = matrix(NA,nrow=ncol(design_mat),ncol=nrow(matA))
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
    
    for(ix_row in 1:nrow(matA)){ #go row-by-row in matA and compare with matB, then reassemble results matrix
      print(paste0("working on row ",ix_row," ..."))
      sse_groups = rep(0,nrow(matB)-ix_row+1)
      corrs_rg = matrix(NA,nrow=ncol(design_mat),ncol=(nrow(matB)-ix_row+1))
      pvals_rg = matrix(NA,nrow=ncol(design_mat),ncol=(nrow(matB)-ix_row+1))
      classes_rg = matrix(NA,nrow=ncol(design_mat),ncol=(nrow(matB)-ix_row+1))
      for(groupn in 1:ncol(design_mat)){
        #calculate stats for row
        matAig = matA[ix_row,design_mat[,groupn]==1,drop=FALSE] #get row for matA
        matBg = matB[ix_row:nrow(matB),design_mat[,groupn]==1,drop=FALSE]
        
        varA = varsA_group[groupn,ix_row,drop=FALSE]
        varB = varsB_group[groupn,ix_row:nrow(matB),drop=FALSE]
        covsAB = cov(t(matAig),t(matBg))
        
        if(corrType=="bicor"){
        corrs = WGCNA::bicor(t(matAig),t(matBg),use="all.obs")
        } else {
        corrs = WGCNA::cor(t(matAig),t(matBg),use="all.obs",method=corrType) #get Pearson rho, using Cpp fxn from WGCNA
        }
        #calculate p-value for the individual groups
        deg_freedom = ncol(matAig)-2
        in_pt = (abs(corrs) * sqrt(deg_freedom) / sqrt(1 - corrs^2))
        pvals = 2 * (1 - pt(in_pt, deg_freedom))
        corrs_rg[groupn,] = corrs
        pvals_rg[groupn,] = pvals
        classes_rg[groupn,] = sapply(1:length(corrs),function(x){if((pvals[x]<0.05)&(corrs[x]>0)){"+"}else if((pvals[x]<0.05)&(corrs[x]<0)){"-"}else{"0"}})
        sse_groups = sse_groups + (varB-covsAB*(corrs/sqrt(varA%*%(1/varB))))*(ncol(matBg)-1)
      } #end group
      matAi = matA[ix_row,,drop=FALSE] #get row for matA
      matBi = matB[ix_row:nrow(matB),,drop=FALSE]
      
      varA = varsA[ix_row]
      varB = varsB[ix_row:nrow(matB)]
      covsAB = cov(t(matAi),t(matBi))
      if(corrType=="bicor"){
      corrs = WGCNA::bicor(t(matAi),t(matBi),use="all.obs")
      } else {
      corrs = WGCNA::cor(t(matAi),t(matBi),use="all.obs",method=corrType)
      }
      sse_all = (varB-covsAB*(corrs/sqrt(varA%*%(1/varB))))*(ncol(matB)-1)
      df1 <- ncol(design_mat)
      df2 <- nrow(design_mat)-2*(df1)
      f <- (sse_all-sse_groups)*df2/(df1*sse_groups)
      p <- pf(f, df1, df2, lower.tail=FALSE)
      corrs_mat[ix_row,ix_row:nrow(matB)] = apply(t(corrs_rg),1,paste,collapse="/")
      pvals_mat[ix_row,ix_row:nrow(matB)] = apply(t(pvals_rg),1,paste,collapse="/")
      classes_mat[ix_row,ix_row:nrow(matB)] = apply(t(classes_rg),1,paste,collapse="/")
      results_mat[ix_row,ix_row:nrow(matB)] = p
      #list(classes=apply(t(classes_rg),1,paste,collapse="/"),pvalues=p)
    } #end row
  } #end else
  diag(results_mat) = 1
  output = list(corrs=corrs_mat,corrsP=pvals_mat,pvalues=results_mat,classes=classes_mat,secondMat=secondMat)
  return(output)
}
