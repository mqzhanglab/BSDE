.onLoad <- function(libname, pkgname) {
  tryCatch({
    reticulate::use_condaenv("r-reticulate", required = TRUE)
    reticulate::py_run_file(system.file("op_functions.py", package = "BSDE"))
    packageStartupMessage("Python environment loaded.")
    # enable parallelism
    n.cores <- max(parallel::detectCores() - 2, 1)
    doParallel::registerDoParallel(cores = n.cores)
    packageStartupMessage(sprintf("Running in parallel on %d cores.", n.cores))
  }, error = function(.e) {
    message(.e)
    message("\nPython environment cannot be loaded.\n")
    .install <- utils::askYesNo("Do you want to set up a conda environment for Python dependency?\nBefore you proceed, please make sure conda (or minoconda) is installed.\nSee https://docs.conda.io/en/latest/miniconda.html")
    if (.install) {
      tryCatch({
        reticulate::conda_create("r-reticulate", python_version = 3.7)
        reticulate::conda_install("r-reticulate", c("numpy", "cython"))
        reticulate::conda_install("r-reticulate", "POT")
        message("Environment `r-reticulate' is set up. Try `library(BSDE)' again.")
      }, error = function(.e) {
        stop("Environment cannot be loaded.\nSee installation instructions: https://github.com/mqzhanglab/BSDE")
      })
    }
  })
}
