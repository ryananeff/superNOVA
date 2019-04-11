## For configuration instructions and examples, see the batchtools R package reference on CRAN. 

cluster.functions = makeClusterFunctionsLSF(system.file("config/minerva_LSF.tmpl",package="DGCA")
                                            ,scheduler.latency = 1)
libpath = system.file("/../",package="DGCA")
default.resources = list(queue="premium", walltime=60, memory=2000, cores=1, project="acc_zhangb03a",libpath = libpath)
raise.warnings = FALSE
staged.queries = TRUE
