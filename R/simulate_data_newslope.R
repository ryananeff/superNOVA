#' @title Simulate differentially linearly coexpressed genes between two groups, where only the slope is changed.
#' @description Simulates sets of genes coexpressed with an eigengene which show differential coexpression among two groups. This tool simulates all 9 possible classes of possible gene coexpression changes (positive, negative, null correlation in group A and group B).
#' @param samples_per_group The number of samples per group to simulate. Default=50
#' @param total_genes The number of genes among all 9 differential coexpression classes to simulate. Will actually return a number of genes which is a multiple of 9 that is less than this number.
#' @param maxCor_true The absolute value of the maximum correlation rho value for highly correlated pairs.
#' @param minCor_true The absolute value of the minimum correlation rho value for highly correlated pairs.
#' @param maxCor_null The absolute value of the maximum correlation rho value for null correlated pairs. The minimum is 0 correlation.
#' @param old_slope The slope of the first group (control)
#' @param new_slope The slope of the second group (experimental)
#' @param verbose Whether or not to print additional information about the simulation when it runs.
#' @return Returns a list of inputs for chowCor, including a matA, matB, design_mat, subgroup names (conditions), and the actual classes of differential coexpression.
#' @keywords superNOVA simulate slope
#' @importFrom stats model.matrix
#' @export
simulate_data_newslope <- function(samples_per_group,total_genes,
                                   maxCor_true,minCor_true,maxCor_null,
                                   old_slope=1, new_slope=1, verbose=T){

  samples = samples_per_group

  if(verbose){
    print("=============")
    print(paste("samples:",samples*2))
    print(paste("maxCor_true:",maxCor_true))
    print(paste("minCor_true:",minCor_true))
    print(paste("maxCor_null:",maxCor_null))
  }

  int <- function(n){as.integer(n)}

  genes_per_group = int(total_genes/9)

  #linear uniform right now, but we should try other relationship types (anscombe's quartet)
  eigengeneg1 = stats::rnorm(samples) #linear uniform right now, but we should try other relationship types (anscombe's quartet)
  eigengeneg2 = stats::rnorm(samples)

  #the second gene in all of the gene pairs
  pos_corg1 = WGCNA::simulateModule(eigengeneg1,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL) #module is positively correlated
  neg_corg1 = WGCNA::simulateModule(eigengeneg1,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1) #module is negatively correlated
  null_corg1a = WGCNA::simulateModule(eigengeneg1,genes_per_group/2,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL) #no correlation
  null_corg1b = WGCNA::simulateModule(eigengeneg1,genes_per_group/2,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1) #no correlation
  null_corg1 = cbind(null_corg1a,null_corg1b)

  true_posg1 = attributes(pos_corg1)$trueKME
  true_negg1 = attributes(neg_corg1)$trueKME
  true_nullg1a = attributes(null_corg1a)$trueKME
  true_nullg1b = attributes(null_corg1b)$trueKME*-1
  true_nullg1 = c(true_nullg1a,true_nullg1b)

  pos_corg2 = WGCNA::simulateModule(eigengeneg2,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL) #module is positively correlated
  neg_corg2 = WGCNA::simulateModule(eigengeneg2,genes_per_group,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1) #module is negatively correlated
  null_corg2a = WGCNA::simulateModule(eigengeneg2,genes_per_group/2,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL) #no correlation
  null_corg2b = WGCNA::simulateModule(eigengeneg2,genes_per_group/2,nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=F,geneMeans=NULL,propNegativeCor = 1) #no correlation
  null_corg2 = cbind(null_corg2a,null_corg2b)

  true_posg2 = attributes(pos_corg2)$trueKME
  true_negg2 = attributes(neg_corg2)$trueKME
  true_nullg2a = attributes(null_corg2a)$trueKME
  true_nullg2b = attributes(null_corg2b)$trueKME*-1
  true_nullg2 = c(true_nullg2a,true_nullg2b)

  order_control = cbind(pos_corg1,pos_corg1,pos_corg1,
                        null_corg1,null_corg1,null_corg1,
                        neg_corg1,neg_corg1,neg_corg1) #samples on rows, genes on columns, so bind columns
  
  new_slope = new_slope/old_slope

  order_disease = cbind(pos_corg2*new_slope, pos_corg2*new_slope, pos_corg2*new_slope,
                        null_corg2,null_corg2,null_corg2,
                        neg_corg2*new_slope, neg_corg2*new_slope, neg_corg2*new_slope) #samples on rows, genes on columns, so bind columns

  #order_disease = order_disease * new_slope #this changes the slope!!
  #order_disease[]

  matB = data.frame(rbind(order_control, order_disease))
  eigengeneg1 = as.matrix(eigengeneg1)
  eigengeneg2 = as.matrix(eigengeneg2)
  matA = data.frame(rbind(eigengeneg1,eigengeneg2))

  genenames = apply(as.matrix(1:ncol(matB)),1,function(x) paste0("gene",x))
  samplenames = apply(as.matrix(1:nrow(matA)),1,function(x) paste0("sample",x))
  colnames(matA) = c("eigengene")
  colnames(matB) = genenames
  rownames(matA) = samplenames
  rownames(matB) = samplenames

  matA = t(matA)
  matB = t(matB)

  ## calculate true pairs classes and control/exp true correlations
  colnames_truth = c("Gene1","Gene2","corA","corB","true_class","true_signif")
  true_pairs = data.frame(matrix(nrow = nrow(matB),ncol=length(colnames_truth)))
  colnames(true_pairs) = colnames_truth #you can't handle the truth
  count = 1 #this is the row no. in the true pairs array
  # positive to positive
  for(i in 1:3){
    for(iter in 1:length(true_posg1)){
      gene_name = paste0("gene",count) #this is the gene name from genenames
      corA = true_posg1[iter]
      corB = true_posg2[iter]
      true_pairs[count,] = c("eigengene",gene_name,corA,corB,"+/+",T)
      count = count + 1
    }
  }

  #null to null
  for(i in 1:3){
    for(iter in 1:length(true_posg1)){
      gene_name = paste0("gene",count) #this is the gene name from genenames
      corA = true_nullg2[iter]
      corB = true_nullg2[iter]
      true_pairs[count,] = c("eigengene",gene_name,corA,corB,"0/0",F)
      count = count + 1
    }
  }

  #negative to negative
  for(i in 1:3){
    for(iter in 1:length(true_negg1)){
      gene_name = paste0("gene",count) #this is the gene name from genenames
      corA = -1*true_negg1[iter]
      corB = -1*true_negg2[iter]
      true_pairs[count,] = c("eigengene",gene_name,corA,corB,"-/-",T)
      count = count + 1
    }
  }

  true_pairs[["corA"]] = as.numeric(true_pairs[["corA"]])
  true_pairs[["corB"]] = as.numeric(true_pairs[["corB"]])
  rownames(true_pairs) = paste0(true_pairs$Gene1,"_",true_pairs$Gene2)

  design_mat = data.frame(rbind(cbind(rep(1,samples),rep(0,samples)),cbind(rep(0,samples),rep(1,samples))))
  conditions = c("group_control","group_exp")
  colnames(design_mat) = conditions
  design_mat = model.matrix(~0+group_control+group_exp,design_mat)
  colnames(design_mat) = conditions

  matA_in = matA
  matB_in = matB

  return(list(matA=matA_in,
              matB=matB_in,
              design_mat=design_mat,
              conditions=conditions,
              true_pairs=true_pairs))
}