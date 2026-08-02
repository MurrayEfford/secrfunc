###############################################################################
## package 'secrfunc'
## onLoad.R
## 2026-08-03
###############################################################################

.onLoad <- function (libname, pkgname) {
    ## also sets environment variable RCPP_PARALLEL_NUM_THREADS
    defaultncores <- RcppParallel::defaultNumThreads()
    if (defaultncores == 1) {
        RcppParallel::setThreadOptions(1)
    }
    else {
        RcppParallel::setThreadOptions(2)
    }
}

## .onLoad is preferred if actions are required for single functions 
## that may be called without attaching package