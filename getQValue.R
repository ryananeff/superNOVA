#' @title Calculates qvalues for any kind of data, failing gracefully. Wrapper for StoreyLab/qvalue function.
#' @description Runs the qvalue calculation step
#' @param pvalues A list of pvalues
#' @return A named list or qobj containing the qvalues at qobj$qvalues, depending on whether the calculation was successful.
#' @export
getQValue <- function(pvalues){
  try_lambda = 20
  qobj = tryCatch(
      {
        if ( max(0.01,round(min(pvalues)+0.01,2)) >= min(round(max(pvalues)-0.01,2),0.95) ) {
             qobj = list()
             qobj$qvalues = rep(NA, length.out=length(pvalues))
             return(qobj)
        }else{
        rangevals = (max(pvalues)-min(pvalues))
        qvalue::qvalue(p = pvalues, lambda = seq(max(0.01,round(min(pvalues)+0.01,2)), 
                                                 min(round(max(pvalues)-0.01,2),0.95), 
                                                 min(0.01,round(rangevals/try_lambda,4))), lfdr.out=FALSE)
        }
      }, error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        cat("\n")
        message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
        #if the qvalue computation returned without error, then its format should be a list; if not, there was an error.
        qobj = tryCatch(
          {
            qvalue::qvalue(p = pvalues, lambda = seq(0.1, 0.9, 0.01),lfdr.out=FALSE)
          }, error=function(cond) {
            message("Here's the original error message:")
            message(cond)
            cat("\n")
            message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
            qobj = tryCatch(
              {
                qvalue::qvalue(p = pvalues, lambda = seq(0.2, 0.8, 0.01),lfdr.out=FALSE)
              }, error=function(cond) {
                message("Here's the original error message:")
                message(cond)
                cat("\n")
                message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence...")
                qobj = tryCatch(
                  {
                    qvalue::qvalue(p = pvalues, lambda = seq(0.3, 0.7, 0.01),lfdr.out=FALSE)
                  }, error=function(cond) {
                    message("Here's the original error message:")
                    message(cond)
                    cat("\n")
                    message("estimated pi0 <= 0 sometimes happens with relatively small numbers of gene pairs. Using a more conservative lambda sequence... if this doesn't work, will report the empirical p-values and the adjusted q-values as NA values to indicate that q-value adjustment did not work.")
                    qobj = list()
                    qobj$qvalues = rep(NA, length.out=length(pvalues))
                    return(qobj)
                })
            })
        })
      })
  return(qobj)
}