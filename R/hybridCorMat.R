#' @title Create a hybrid correlation matrix for each subgroup
#' @description Create a hybrid correlation matrix from a superNOVA output, suitable for downstream MEGENA network creation. Differentially coexpressed edges use a subgroup-specific coexpression model, while others use a global model.
#' @param input_file An input filename, pointing to a superNOVA output file.
#' @param output_file_prefix The prefix to give to all of the subgroup-specific output files.
#' @param diffCor_cutoff The p-value at which to consider a subgroup-specific model versus a global model. Default = 0.05.
#' @param corPval_cutoff The p-value at which to consider an edge to have significant correlation between genes versus no correlation. Default = 0.05.
#' @param bufflen The number of buffer rows for reading the input file, increasing this will speed up the computation at the expense of memory. Default = 10000.
#' @keywords superNOVA hybrid
#' @return None (output to files)
#' @export
hybridCorMat <- function(input_file, output_file_prefix, diffCor_cutoff = 0.05, corPval_cutoff = 0.05, bufflen = 10000){

  fp = file(input_file,"rt")

  continue = TRUE
  header = TRUE
  need_subtypes = TRUE
  subtype_files  = c()
  sheader = NULL
  while(continue){
    buffer = readLines(con=fp,n=bufflen,ok=TRUE)
    if(length(buffer)!=bufflen) { continue=FALSE }
    for (line in buffer){
      splitline = unlist(strsplit(trimws(line), '\t'))

      if(header){ #to ignore the header line
        header=FALSE
        sheader = setNames(1:length(splitline),splitline)
        if(!("Gene1" in splitline)){
          stop("Please ensure your input file has a header in the first row from superNOVA.")
        }
        next
      }

      subtypes = unlist(strsplit(splitline[sheader[["group_order"]]],'/'))

      if(need_subtypes){
        subtype_files = list()
        for(ix in 1:length(subtypes)){
          subtype = subtypes[ix]
          filen = paste0(output_file_prefix,subtype,".ijw.tsv")
          con = file(filen,"wt", blocking=F)
          subtype_files[[ix]] = con
        }
        need_subtypes = FALSE
      }

      row = splitline[sheader[["Gene1"]]]
      col = splitline[sheader[["Gene2"]]]

      pValDiff = as.numeric(splitline[sheader[["pValDiff"]]]) #supernova diffcor pvalue

      if(pValDiff < diffCor_cutoff){
        groupCors = as.numeric(unlist(strsplit(splitline[sheader[["groupCor"]]],'/')))
        groupCorsPvals = as.numeric(unlist(strsplit(splitline[sheader[["groupCorPval"]]],'/'))) #subgroup cor pvalue for pair
        for(ix in 1:length(subtype_files)){
          if(groupCorsPvals[ix] < corPval_cutoff){
            writeLines(paste(row,col,groupCors[ix],sep='\t'),con=subtype_files[[ix]]) #write subtype-specific rho to file
          }
        }
      } else {
        globalCor = as.numeric(splitline[sheader[["globalCor"]]])
        globalCorP = as.numeric(splitline[sheader[["globalCorP"]]]) #global cor pvalue for pair
        if(globalCorP < corPval_cutoff){
          for(ix in 1:length(subtype_files)){
            writeLines(paste(row,col,globalCor,sep='\t'),con=subtype_files[[ix]]) #write global rho to all files
          } #end for
        } #end if
      } #end ifelse pval checking
    } #end for line in buffer
  } #end while

  #cleanup
  for(ix in 1:length(subtype_files)){close(subtype_files[[ix]])}
  close(fp)
}
