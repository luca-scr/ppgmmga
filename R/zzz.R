.onAttach <- function(lib, pkg)
{
  unlockBinding(".ppgmmga.default", asNamespace("ppgmmga"))
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage(rep(" ", 13),"Welcome to ppgmmga"," " ,"version ",version)
  packageStartupMessage("Projection pursuit based on Gaussian mixtures and genetic algorithm for data visualisation and anomaly detection","\n")
  invisible()
}

.onUnload <- function (lib) {
  library.dynam.unload("ppgmmga", lib)
}


