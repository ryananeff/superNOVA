#' @title SuperNOVA subgroup-specific hybrid MEGENA networks
#' @description Takes in an subgroup-specific hybrid correlation matrix from hybridCorMat and builds a MEGENA network from it.
#' @param data.file Path to the hybrid correlation matrix.
#' @param subtype.name Name of the subgroup for the output.
#' @param module.pval Module significance p-value. Recommended is 0.05 (default).
#' @param hub.pval Connectivity significance p-value based random tetrahedral networks. Recommended is 0.05 (default).
#' @param cor.perm Number of permutations for calculating FDRs for all correlation pairs. Default is 10.
#' @param hub.perm Number of permutations for calculating connectivity significance p-value. Default is 50.
#' @param n.cores Number of cores to use with MEGENA. Default is 1 (single-threaded).
#' @return None (output to files)
#' @keywords superNOVA MEGENA hybrid
#' @import MEGENA doParallel parallel igraph foreach
#' @export
snMEGENA <- function(data.file,subtype.name,
                     module.pval = 0.05, hub.pval = 0.05, cor.perm=10,
                     hub.perm=50, n.cores=1){

  annot.table = NULL
  id.col = 1
  symbol.col= 2

  wkdir = getwd()

  if(n.cores>1){ doPar=T } else { doPar=F }

  # create output directory
  out.dir <- paste(sep = "",wkdir,"/MEGENA_subtype_hybrid_",subtype.name);

  # load input data
  ijw = read.csv(data.file,header = TRUE,sep="\t") #output from hybridCorMat

  dir.create(out.dir)
  setwd(out.dir)

  if (doPar & getDoParWorkers() == 1) #create parallel worker backend
  {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:",foreach::getDoParWorkers(),"\n",sep = "")); flush.console()
  }

  el <- MEGENA::calculate.PFN(ijw,doPar = doPar,num.cores = n.cores,keep.track = TRUE)
  g <- igraph::graph.data.frame(el,directed = FALSE)

  MEGENA.output <- MEGENA::do.MEGENA(g,
                   mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                   min.size = 10,max.size = vcount(g)/2,
                   doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
                   save.output = TRUE)

  if (getDoParWorkers() > 1) #remove parallel workers we don't need anymore
  {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }

  summary.output <- MEGENA::MEGENA.ModuleSummary(MEGENA.output,
                                         mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                         min.size = 10,max.size = vcount(g)/2,
                                         annot.table = annot.table,
                                         id.col = id.col,
                                         symbol.col = symbol.col,
                                         output.sig = TRUE)

  igraph::write_graph(g,file=paste(out.dir,"/",subtype.name,".graph.graphml",sep=""),format = "graphml")
  igraph::write_graph(g,file=paste(out.dir,"/",subtype.name,".graph.tsv",sep=""),format = "edgelist")
  write.table(el,file=paste(out.dir,"/",subtype.name,".megena_out.tsv",sep=""),sep="\t")
  MEGENA::output.geneSet.file(summary.output$modules,paste(out.dir,"/",subtype.name,".multiscale_significant.modules.txt",sep=""))
  summary.table <- summary.output$module.table
  write.table(summary.table,file = paste(out.dir,"/",subtype.name,".module_summary.txt",sep=""),sep = "\t",row.names = F,col.names = T,quote = F)
  setwd(wkdir)
}
