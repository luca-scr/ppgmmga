# Default options used in 'ppgmmga' package

.ppgmmga <- list(
  # mclust 
  modelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV"),
  G = 1:9,
  initMclust = "SVD",
  # GA
  popSize = 100,
  pcrossover = 0.8,
  pmutation = 0.1,
  maxiter = 1000,
  run = 100,
  selection = GA::gareal_lsSelection,
  crossover = GA::gareal_laCrossover,
  mutation = GA::gareal_raMutation,
  parallel = FALSE,
  numIslands = 4,
  migrationRate = 0.1,
  migrationInterval = 10,
  # GA - local search
  optim = TRUE, 
  optimPoptim = 0.05,
  optimPressel = 0.5,
  optimMethod = "L-BFGS-B",
  optimMaxit = 100,
  # ppgmmga
  orth = "QR"
)

ppgmmga.options <- function(...)
{
  current <- get(".ppgmmga", envir = asNamespace("ppgmmga"))
  if(nargs() == 0) return(current)
  args <- list(...)
  if(length(args) == 1 && is.null(names(args))) 
  { arg <- args[[1]]
    switch(mode(arg),
           list = args <- arg,
           character = return(.ppgmmga[[arg]]),
           stop("invalid argument: ", dQuote(arg)))
  }
  if(length(args) == 0) return(current)
  n <- names(args)
  if(is.null(n)) stop("options must be given by name")
  current[n] <- args
  # browser()
  if(sys.nframe() == 1) 
    assign(".ppgmmga", current, envir = asNamespace("ppgmmga"))
  invisible(current)
}

