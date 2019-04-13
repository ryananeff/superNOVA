#' @title Worker thread for chowParallel (internal function).
#' @description Runs chowCor between two or more subgroups of samples in the data versus the global model on a worker thread for two blocks of data (A and B).
#' @param job A batchtools job instance.
#' @param data A named list containing the program kwargs.
#' @param instance Required by batchtools, but not used currently.
#' @return Results from running chowCor() on data block A and B. See chowCor() documentation for more information.
#' @keywords superNOVA
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

	options(warn=data$batchWarningLevel)

	cl<-parallel::makeCluster(data$n.cores)
	doParallel::registerDoParallel(cl)

	set.seed(data$seed) #random seed for reproducibility


	cat("Starting run now...\n")
	cat(paste0(Sys.time(),"\n"))
	if(rownames(instance$matA)==rownames(instance$matB)){
		chow_result = superNOVA::chowCor(design_mat=data$design,matA=t(instance$matA), matB=NULL, compare = data$compare,
		                      corrType=data$corrType)
	}else{
		chow_result = superNOVA::chowCor(design_mat=data$design,matA=t(instance$matA), matB=t(instance$matB), compare = data$compare,
		                      corrType=data$corrType)
	}
	ddcor_res = superNOVA::flattenChow(chow_result)
	if(data$verbose){
		cat("Completed run\n")
		cat(paste0(Sys.time(),"\n"))
	}
	return(ddcor_res) #algorithm
}
