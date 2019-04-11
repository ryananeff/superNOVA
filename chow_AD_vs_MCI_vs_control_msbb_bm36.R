## module load R/3.4.3
rm(list = ls(all.names = TRUE))
set.seed(12345)

library(WGCNA)
disableWGCNAThreads()

options(stringsAsFactors = FALSE)
source("/sc/orga/projects/zhangb03a/neffr01/DGCA_MCI_control_U54/chowCor.R")

message("CHOW COR! hooray!")
message("loading data")
message(Sys.time())

#read data files and load it in

#### 1. expression file

data.file = "/sc/orga/projects/zhangb03a/shared/msbb-seq/Data/mixed_model_correction_updated/msbb.BM_36.PMI_AOD_race_sex_RIN_exonicRate_rRnaRate_batch_adj.tsv"
data.Df <- read.delim(file = data.file,sep = "\t",header = T)
datExpr <- as.matrix(data.Df);

rownames(datExpr) <- as.character(data.Df[,"Symbol"]);
datExpr = datExpr[,-1]
datExpr = datExpr[complete.cases(datExpr),]
rownames_exp = rownames(datExpr)
datExpr = apply(datExpr, 2, as.numeric)
rownames(datExpr) = rownames_exp

design_file = "msbb.meta.BM_36.tsv"
design_mat = read.delim(file=design_file, sep="\t", header=T)
rownames(design_mat) = design_mat[,"Sampleid"]
datExpr = datExpr[,rownames(design_mat),drop=FALSE]
design_mat = design_mat[colnames(datExpr),,drop=FALSE]

#create design matrix
design_mat = model.matrix(~0+control+MCI_1+AD_2,data=data.frame(design_mat))
	
print(paste0("starting for chowCor..."))
message(Sys.time())

chow_result = chowCor(design_mat,datExpr)
print("flattening...")
ddcor_res = flattenChow(chow_result)
print(paste0("completed for chowCor!! writing out..."))
message(Sys.time())

write.table(ddcor_res,file="chowCor_AD_vs_MCI_vs_control_msbb_bm36.tsv", sep = "\t",quote = F,row.names = F)

print(paste0("finished writing out."))
message(Sys.time())


