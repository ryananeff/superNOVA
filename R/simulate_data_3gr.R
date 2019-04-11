#' @title Simulate differentially linearly coexpressed genes between three groups.
#' @description Simulates sets of genes coexpressed with an eigengene which show differential coexpression among three groups. This tool simulates 9 possible classes of possible gene coexpression changes (positive, negative, null correlation in group A, group B, and groupC) but does not test all 27 possible classes.
#' @param samples_per_group The number of samples per group to simulate. Default=50
#' @param total_genes The number of genes among 9 differential coexpression classes to simulate. Will actually return a number of genes which is a multiple of 9 that is less than this number.
#' @param maxCor_true The absolute value of the maximum correlation rho value for highly correlated pairs.
#' @param minCor_true The absolute value of the minimum correlation rho value for highly correlated pairs.
#' @param maxCor_null The absolute value of the maximum correlation rho value for null correlated pairs. The minimum is 0 correlation.
#' @param verbose Whether or not to print additional information about the simulation when it runs.
#' @param corr_slope The slope of the correlations to simulate. All correlated pairs will have this slope.
#' @param corr_intercept The correlation intercept of the genes in matrix B.
#' @param group1_offset The correlation intercept of group1 (e.g. its mean expression value)
#' @param group2_offset The correlation intercept of group2 (e.g. its mean expression value)
#' @param group3_offset The correlation intercept of group3 (e.g. its mean expression value)
#' @return Returns a list of inputs for chowCor, including a matA, matB, design_mat, subgroup names (conditions), and the actual classes of differential coexpression.
#' @keywords superNOVA simulate three
#' @importFrom stats model.matrix
#' @export
simulate_data_3gr <- function(samples_per_group=50,total_genes=2000,
                              maxCor_true=0.95,minCor_true=0.3,maxCor_null=0.2,verbose=T,
                              corr_slope = 1, corr_intercept = 0, group1_offset=0, group2_offset=0,
                              group3_offset=0){

  samples = samples_per_group

  if(verbose){
    print("=============")
    print(paste("samples:",samples*3))
    print(paste("maxCor_true:",maxCor_true))
    print(paste("minCor_true:",minCor_true))
    print(paste("maxCor_null:",maxCor_null))
  }

  int <- function(n){as.integer(n)}

  genes_per_group = int(total_genes/9)

  #linear uniform right now, but we should try other relationship types (anscombe's quartet)
  eigengeneg1 = stats::rnorm(samples) + group1_offset #linear uniform right now, but we should try other relationship types (anscombe's quartet)
  eigengeneg2 = stats::rnorm(samples) + group2_offset
  eigengeneg3 = stats::rnorm(samples) + group3_offset

  #group1
  pos_corg1 = WGCNA::simulateModule(eigengeneg1,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL)*corr_slope+corr_intercept #module is positively correlated
  neg_corg1 = WGCNA::simulateModule(eigengeneg1,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1)*corr_slope+corr_intercept #module is negatively correlated
  null_corg1a = WGCNA::simulateModule(eigengeneg1,genes_per_group/2,nNearGenes=0,minCor=0.01,maxCor=maxCor_null,corPower=1,signed=T,geneMeans=NULL)*corr_slope+corr_intercept #no correlation
  null_corg1b = WGCNA::simulateModule(eigengeneg1,genes_per_group/2,nNearGenes=0,minCor=0.01,maxCor=maxCor_null,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1)*corr_slope+corr_intercept #no correlation
  null_corg1 = cbind(null_corg1a,null_corg1b)

  true_posg1 = attributes(pos_corg1)$trueKME
  true_negg1 = attributes(neg_corg1)$trueKME
  true_nullg1a = attributes(null_corg1a)$trueKME
  true_nullg1b = attributes(null_corg1b)$trueKME*-1
  true_nullg1 = c(true_nullg1a,true_nullg1b)

  #group2
  pos_corg2 = WGCNA::simulateModule(eigengeneg2,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL)*corr_slope+corr_intercept #module is positively correlated
  neg_corg2 = WGCNA::simulateModule(eigengeneg2,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1)*corr_slope+corr_intercept #module is negatively correlated
  null_corg2a = WGCNA::simulateModule(eigengeneg2,genes_per_group/2,nNearGenes=0,minCor=0.01,maxCor=maxCor_null,corPower=1,signed=T,geneMeans=NULL)*corr_slope+corr_intercept #no correlation
  null_corg2b = WGCNA::simulateModule(eigengeneg2,genes_per_group/2,nNearGenes=0,minCor=0.01,maxCor=maxCor_null,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1)*corr_slope+corr_intercept #no correlation
  null_corg2 = cbind(null_corg2a,null_corg2b)

  true_posg2 = attributes(pos_corg2)$trueKME
  true_negg2 = attributes(neg_corg2)$trueKME
  true_nullg2a = attributes(null_corg2a)$trueKME
  true_nullg2b = attributes(null_corg2b)$trueKME*-1
  true_nullg2 = c(true_nullg2a,true_nullg2b)

  #group3
  pos_corg3 = WGCNA::simulateModule(eigengeneg3,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL)*corr_slope+corr_intercept #module is positively correlated
  neg_corg3 = WGCNA::simulateModule(eigengeneg3,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1)*corr_slope+corr_intercept #module is negatively correlated
  null_corg3a = WGCNA::simulateModule(eigengeneg3,genes_per_group/2,nNearGenes=0,minCor=0.01,maxCor=maxCor_null,corPower=1,signed=T,geneMeans=NULL)*corr_slope+corr_intercept #no correlation
  null_corg3b = WGCNA::simulateModule(eigengeneg3,genes_per_group/2,nNearGenes=0,minCor=0.01,maxCor=maxCor_null,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1)*corr_slope+corr_intercept #no correlation
  null_corg3 = cbind(null_corg3a,null_corg3b)

  true_posg3 = attributes(pos_corg3)$trueKME
  true_negg3 = attributes(neg_corg3)$trueKME
  true_nullg3a = attributes(null_corg3a)$trueKME
  true_nullg3b = attributes(null_corg3b)$trueKME*-1
  true_nullg3 = c(true_nullg3a,true_nullg3b)

  order_control = cbind(pos_corg1,pos_corg1,pos_corg1,
                        null_corg1,null_corg1,null_corg1,
                        neg_corg1,neg_corg1,neg_corg1) #samples on rows, genes on columns, so bind columns

  order_diseaseA = cbind(pos_corg2,null_corg2,neg_corg2,
                         pos_corg2,null_corg2,neg_corg2,
                         pos_corg2,null_corg2,neg_corg2) #samples on rows, genes on columns, so bind columns

  order_diseaseB = cbind(pos_corg3,null_corg3,null_corg3,
                         pos_corg3,neg_corg3,neg_corg3,
                         pos_corg3,null_corg3,neg_corg3) #samples on rows, genes on columns, so bind columns

  matB = data.frame(rbind(order_control, order_diseaseA, order_diseaseB))
  eigengeneg1 = as.matrix(eigengeneg1)
  eigengeneg2 = as.matrix(eigengeneg2)
  eigengeneg3 = as.matrix(eigengeneg3)
  matA = data.frame(rbind(eigengeneg1,eigengeneg2,eigengeneg3))

  genenames = apply(as.matrix(1:ncol(matB)),1,function(x) paste0("gene",x))
  samplenames = apply(as.matrix(1:nrow(matA)),1,function(x) paste0("sample",x))
  colnames(matA) = c("eigengene")
  colnames(matB) = genenames
  rownames(matA) = samplenames
  rownames(matB) = samplenames

  matA = t(matA)
  matB = t(matB)

  ## calculate true pairs classes and control/exp true correlations
  colnames_truth = c("Gene1","Gene2","corA","corB","corC","true_class","true_signif")
  true_pairs = data.frame(matrix(nrow = nrow(matB),ncol=length(colnames_truth)))
  colnames(true_pairs) = colnames_truth #you can't handle the truth
  count = 1 #this is the row no. in the true pairs array
  # positive to positive to positive
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = true_posg1[iter]
    corB = true_posg2[iter]
    corC = true_posg3[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"+/+/+",F)
    count = count + 1
  }
  #positive to null to null
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = true_posg1[iter]
    corB = true_nullg2[iter]
    corC = true_nullg3[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"+/0/0",T)
    count = count + 1
  }
  #positive to negative to null
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = true_posg1[iter]
    corB = -1*true_negg2[iter]
    corC = true_nullg3[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"+/-/0",T)
    count = count + 1
  }
  #null to posiitve to positive
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = true_nullg1[iter]
    corB = true_posg2[iter]
    corC = true_posg3[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"0/+/+",T)
    count = count + 1
  }
  #null to null to negative
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = true_nullg1[iter]
    corB = true_nullg2[iter]
    corC = -1*true_negg3[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"0/0/-",T)
    count = count + 1
  }
  #null to negative to negative
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = true_nullg1[iter]
    corB = -1*true_negg2[iter]
    corC = -1*true_negg3[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"0/-/-",T)
    count = count + 1
  }
  #negative to posiitve to positive
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = -1*true_negg1[iter]
    corB = true_posg2[iter]
    corC = true_posg2[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"-/+/+",T)
    count = count + 1
  }
  #negative to null to null
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = -1*true_negg1[iter]
    corB = true_nullg2[iter]
    corC = true_nullg3[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"-/0/0",T)
    count = count + 1
  }
  #negative to negative to negative
  for(iter in 1:length(true_posg1)){
    gene_name = paste0("gene",count) #this is the gene name from genenames
    corA = -1*true_negg1[iter]
    corB = -1*true_negg2[iter]
    corC = -1*true_negg3[iter]
    true_pairs[count,] = c("eigengene",gene_name,corA,corB,corC,"-/-/-",F)
    count = count + 1
  }

  true_pairs[["corA"]] = as.numeric(true_pairs[["corA"]])
  true_pairs[["corB"]] = as.numeric(true_pairs[["corB"]])
  true_pairs[["corC"]] = as.numeric(true_pairs[["corC"]])
  rownames(true_pairs) = paste0(true_pairs$Gene1,"_",true_pairs$Gene2)

  #======
  #TODO: concat the matrices for control and disease, make matrices with information about the correlation before and after, etc.
  #======
  design_mat = data.frame(rbind(cbind(rep(1,samples),rep(0,samples),rep(0,samples)),
                                cbind(rep(0,samples),rep(1,samples),rep(0,samples)),
                                cbind(rep(0,samples),rep(0,samples),rep(1,samples))
  ))
  conditions = c("group_control","group_expA","group_expB")
  colnames(design_mat) = conditions
  design_mat = model.matrix(~0+group_control+group_expA+group_expB,design_mat)
  colnames(design_mat) = conditions

  matA_in = matA
  matB_in = matB

  return(list(matA=matA_in,
              matB=matB_in,
              design_mat=design_mat,
              conditions=conditions,
              true_pairs=true_pairs))
}
