#' @title Worker thread for ddcorAllParallel (internal function).
#' @description Runs the discovery of differential correlation (ddcor) section for comparing pairwise correlations across conditions in the Differential Gene Correlation Analysis (DGCA) package on a worker thread for two blocks of data (A and B).
#' @param job A batchtools job instance.
#' @param data A named list containing the program kwargs.
#' @param instance Required by batchtools, but not used currently.
#' @return Typically, the returned object is a data frame of the table of differential correlations between conditions. In the case that dCorAvg is calculated, the returned object is instead a list containing that table as well as the object summarizing the difference in average correlation for the specified portion of the data set.
#' @export
chowParallelWorker <- function(job,data,instance){

	# Required fields in data:
	# matA
	# matB
	# design_mat
	# compare (present or NULL)
	# seed
	# n.cores
	# batchWarningLevel
	# libloc
	
	source("chowCor.R")
	source("getQValue.R")
	
	options(warn=data$batchWarningLevel)

	cl<-parallel::makeCluster(data$n.cores)
	doParallel::registerDoParallel(cl)

	set.seed(data$seed) #random seed for reproducibility

	if(data$verbose){
		cat("Starting run now...\n")
		cat(paste0(Sys.time(),"\n"))
	}
	if(rownames(instance$matA)==rownames(instance$matB)){
		chow_result = chowCor(design_mat=data$design,matA=instance$matA, matB=NULL, compare = data$compare,
		                      corrType=data$corrType)
	}else{
		chow_result = chowCor(design_mat=data$design,matA=instance$matA, matB=instance$matB, compare = data$compare,
		                      corrType=data$corrType)
	}
	ddcor_res = flattenChow(chow_result)
	#recalculate the q-values two ways
	if(data$verbose){
		message("calculating qvalues now.")
		cat(paste0(Sys.time(),"\n"))
	}
	ddcor_res[,"qValDiff"]=as.matrix(getQValue(ddcor_res$pValDiff)$qvalues)
	if(data$verbose){
		cat("Completed run\n")
		cat(paste0(Sys.time(),"\n"))
	}
	return(ddcor_res) #algorithm
}