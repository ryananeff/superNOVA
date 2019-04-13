#' @title Calls the superNOVA pipeline, splitting input matrix into multiple batch jobs on an HPC.
#' @description Evaluates differential coexpression between two or more subgroups of samples in the data versus the global model, using multiple nodes in parallel in a batch environment.
#' @param inputMat The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your compare.
#' @param design A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.
#' @param compare Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.
#' @param outputFile Location to save the output. Required.
#' @param sigOutput Should we save the significant results in a separate file? Default = FALSE.
#' @param corrType The correlation type of the analysis, limited to "pearson","spearman",or "bicor". Default = "pearson".
#' @param sigThresh This numeric value specifies the p-value threshold at which a differential correlation p-value is deemed significant for differential correlation class calculation. Default = 1, as investigators may use different cutoff thresholds; however, this can be lowered to establish significant classes as desired.
#' @param verbose Option indicating whether the program should give more frequent updates about its operations. Default = FALSE.
#' @param perBatch Number of times to split the features of the input data into separate batches. A higher number creates a larger number of jobs, but may be less uniform. Default = 10.
#' @param coresPerJob Number of cores to use on each batch job run. Default = 2.
#' @param timePerJob Walltime to request for each batch job (e.g. in a HPC cluster), in minutes. Default = 60
#' @param memPerJob Memory to request for each batch job (e.g. in a HPC cluster), in MB. Default = 2000
#' @param batchConfig Location of the batchtools configuration file (e.g. to configure this tool to work with your HPC cluster). Defaults to one used at inst/config/batchConfig_Zhang.R.
#' @param batchDir Location to store temporary files, logs, and results of the batch run. This is the registry for the batchtools R package. Default = batchRegistry/
#' @param batchWarningLevel Warning level on remote nodes during chowCor calculation (equivalent to setting options(warn=batchWarningLevel). Default = 0.
#' @param batchSeed Random seed to use on all batch jobs. Default = 12345.
#' @param maxRetries Number of times to re-submit jobs that failed. This is helpful for jobs that failed due to transient errors on an HPC. Default = 3
#' @param testJob Test one job before running it? Default = FALSE
#' @param chunkSize Execute multiple splits sequentially on each node. Default = 1 (false)
#' @return Returns whether all jobs successfully executed or not. Output is in the output file.
#' @keywords superNOVA
#' @importFrom utils head
#' @import data.table
#' @export

chowParallel <- function(inputMat, design, outputFile, compare=NULL,
    sigOutput = FALSE, sigThresh = 0.05, verbose = FALSE, corrType="pearson",
	perBatch = 10, coresPerJob = 2, timePerJob = 60, memPerJob = 2000,
	batchConfig = system.file("config/batchConfig_Zhang.R",package="superNOVA"), batchDir = "batchRegistry",
	batchWarningLevel = 0, batchSeed = 12345, maxRetries = 3, testJob=FALSE,chunkSize=1){

	## REMOVED PARAMETERS

	plotFdr = FALSE # must be false because batch jobs would make multiple. not implemented yet.
	splitSet = NULL # not implemented.
	inputMatB = NULL # we use this to submit batch jobs, so not implemented yet
	heatmapPlot = FALSE # must be false because batch jobs would make multiple. not implemented yet.
	customize_heatmap = FALSE # must be false because batch jobs would make multiple. not implemented yet.
	color_palette = NULL # must be null because batch jobs would make multiple. not implemented yet.
	heatmapClassic = FALSE # must be false because batch jobs would make multiple. not implemented yet.

	#######
	if(verbose){
		message("Creating batch job registry (will overwrite previous registry at same location)")
		message(Sys.time())
	}
	if(dir.exists(batchDir)){
		system(paste0("rm -r ",batchDir)) #remove the old batchRegistry to make the new one
	}

	## Create simple registry:
	reg <- batchtools::makeExperimentRegistry(file.dir=batchDir, conf.file=batchConfig)
	res = list(measure.memory = TRUE,walltime=timePerJob,memory=memPerJob,cores=coresPerJob,chunks.as.array.jobs = TRUE)

	matrix_part <- function(job,data,startA, endA, startB, endB){
	    blockA = data$input[startA:endA,,drop=FALSE] #get first chunk of genes to compare
	    blockB = data$input[startB:endB,,drop=FALSE] #get second chunk of genes to compare
		nPairs = (nrow(blockA)*nrow(blockA)-nrow(blockA))/2+(nrow(blockB)*nrow(blockB)-nrow(blockB))/2
		list(matA=blockA, matB=blockB, nPairs=nPairs) #problem
	}

	batchtools::addAlgorithm(name="chow",fun=superNOVA::chowParallelWorker,reg=reg)
	batchtools::addProblem(name="input_data",fun=matrix_part,reg=reg,
	                       data = list(input=inputMat, compare=compare,design=design,
	                         n.cores=coresPerJob,corrType=corrType,
	                         verbose=verbose, batchWarningLevel=batchWarningLevel, seed=batchSeed))

	tot_len = nrow(inputMat) #total length of input genes
	split_size = round(tot_len/perBatch) #number of genes per split
	count = 0 #starting at 0
	total_runs = (perBatch**2-perBatch)/2 #total blocks we will run from the input data based on split
	algorithms = list(chow=NA)
	problems = list()

	message("======= INPUT MATRIX =======")
	message(paste0("split size: ",split_size))
	message(paste0("total length: ",tot_len))
	message(paste0("size of matrix: ",((split_size**2-split_size)/2)))

	if (((split_size**2-split_size)/2)>=2**31-1){
	  stop("You are comparing too many genes or doing too many permutations for this tool. Increase `split`.")
	}
	startAs = c()
	endAs = c()
	startBs = c()
	endBs = c()
	for(a in 0:(perBatch-1)){ #genes to compare from
	  for(b in a:(perBatch-1)){ #genes to compare to
	    count = count + 1
	    startA = (a*split_size+1) #start index for genesA...
	    endA = ((a+1)*split_size)
	    startB = (b*split_size+1) #start index for genesB...
	    endB = ((b+1)*split_size)
	    if ((endA>tot_len)&(startA<tot_len)){ #if it's not an even split
	      endA = tot_len
	    }
	    if ((endB>tot_len)&(startB<tot_len)){ #if it's not an even split
	      endB = tot_len
	    }
	    if((startA>=tot_len)|(startB>=tot_len)){ #if we've split too many times, skip
	      next
	    }
	    message(paste0("preparing block ",count,", a: ",a,", b: ",b))
	    startAs = c(startAs,startA)
	    startBs = c(startBs,startB)
	    endAs = c(endAs,endA)
	    endBs = c(endBs,endB)
	  }
	}
	pdes = list(input_data=data.table::data.table(startA=startAs,endA=endAs,startB=startBs,endB=endBs))
	batchtools::addExperiments(pdes) #add the problem x algorithm here

	id1 = head(batchtools::findExperiments(algo.name = "chow"), 1)

	if(testJob){
		message("testing one job that it will run")
		message(Sys.time())
		result <- batchtools::testJob(id = id1)
	}

	message("Submitting jobs to cluster...")
	message(Sys.time())
  ids = batchtools::findExperiments()
	batchtools::submitJobs(ids=ids, resources=res)

	message("Waiting for jobs to complete...")
	message(Sys.time())
	job_retries = sapply(batchtools::findJobs()$job.id, function(x) {0})
	while(length(batchtools::findNotDone()$job.id)>0){
		batchtools::waitForJobs(timeout=60)
		err = batchtools::findErrors()$job.id
		exp = batchtools::findExpired()$job.id
		if (length(err)>0){
			for(i in err){
				job_retries[i] = job_retries[i] + 1
				message(paste0("Found error in job ",i,", restarting, retry attempt ",job_retries[i]))
				print(batchtools::getErrorMessages(i))
				batchtools::submitJobs(i,resources=res)
			}
		}
		if (length(exp)>0){
			for(i in exp){
				job_retries[i] = job_retries[i] + 1
				message(paste0("Found expired job ",i,", restarting with ",1.25**job_retries[i],"x more resources, retry attempt ",job_retries[i]))
				res_job = res
				res_job$memory = round(res$memory*(1.25**job_retries[i]))
				res_job$walltime = round(res$walltime*(1.25**job_retries[i]))
				res_job$cores = round(res$cores*(1.25**job_retries[i]))
				print(batchtools::getLog(i))
				batchtools::submitJobs(i,resources=res_job)
			}
		}
		if(max(job_retries)>=maxRetries){
			message("Maximum number of retries exceeded, stopping jobs...")
			batchtools::killJobs()
			reg <<- reg
			stop("Automatic retry failed, registry available for debugging at `reg`.")
		}
	}

	message("Writing out results in chunks (unsorted)...")
	message(Sys.time())

	for (i in batchtools::findDone()$job.id){ #for now
		message(paste0("Writing chunk ",i," of ",count))
		message(Sys.time())
		result = batchtools::loadResult(i)
	   if (i==1){
		   utils::write.table(result,file = outputFile,sep = "\t",col.names = T,row.names=F,quote=F) #write sequentially to file, create file
		   if(sigOutput){
			   utils::write.table(result[result[,"pValDiff"]<=sigThresh,],file = paste0(outputFile,".signif.txt"),
			               sep = "\t",col.names = T,row.names=F,quote=F)
		   }
	   } else {
		   utils::write.table(result,file = outputFile,append = T, sep = "\t",col.names = F,row.names=F,quote=F) #write sequentially to file, append to file already made
		   if(sigOutput){
			   utils::write.table(result[result[,"pValDiff"]<=sigThresh,],file = paste0(outputFile,".signif.txt"),
			               append = T, sep = "\t",col.names = F,row.names=F,quote=F)
		   }
	   }
	}
	message("Completed!")
	message(Sys.time())
}
