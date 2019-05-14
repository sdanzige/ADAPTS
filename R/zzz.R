#' @importFrom foreach %do% %dopar% getDoParWorkers foreach
#' @importFrom stats na.omit t.test var
#' @importFrom utils tail
.onLoad <- function(libname, pkgname){
  if (.Platform$OS.type == 'unix') {
    doParallel::registerDoParallel(cores = parallel::detectCores())
    options(mc.cores=parallel::detectCores()) 
    options(cores=parallel::detectCores())
  }
}
