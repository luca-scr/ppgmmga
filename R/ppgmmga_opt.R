

.ppgmmga.default = list("Orth" = "QR",
                        "popSize" = 100,
                        "pcr" = 0.8,
                        "pm" = 0.1,
                        "maxiter" = 1000,
                        "run" = 100,
                        "optim" = TRUE,
                        "selection" = gareal_lsSelection,
                        "crossover" = gareal_laCrossover,
                        "mutation" = gareal_raMutation,
                        "parallel" = FALSE,
                        "numIslands" = 4,
                        "migrationRate" = 0.1,
                        "migrationInterval" = 10,
                        "optimMethod" = "L-BFGS-B",
                        "poptim" = 0.05,
                        "pressel" = 0.5,
                        "fnscale" = -1,
                        "maxit" = 100,
                        "initMclust" = "SVD",
                        "ModelNames" = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "EVV", "VVV","VEE"),
                        "G" = 1:9

)

ppgmmga.options <- function(...)
{

  temp <- .ppgmmga.default

  if (nargs()==0)
    return(temp)

  op <- list(...)
  if (length(op) == 1 && is.null(names(op))) {
    arg <- op[[1]]
    switch(mode(arg), list = op <- arg, character = return(.ppgmmga.default[[arg]]),
           stop("invalid argument: ", dQuote(arg)))
  }
  if(any(sapply(op, is.list))){
    opt.n <- names(op[[1]])
    temp[opt.n] <- op[[1]]
  } else {
    opt.n <- names(op)
    temp[opt.n] <- op
  }

  if (sys.parent() == 0) env <- asNamespace("ppgmmga") else env <- parent.frame()
  assign(".ppgmmga.default", temp, envir = env)
  invisible(temp)
}

ppgmmga.option.restore <- function()
{
  restorealg <- list("Orth" = "QR",
                     "popSize" = 100,
                     "pcr" = 0.8,
                     "pm" = 0.1,
                     "maxiter" = 1000,
                     "run" = 100,
                     "optim" = TRUE,
                     "selection" = gareal_lsSelection,
                     "crossover" = gareal_laCrossover,
                     "mutation" = gareal_raMutation,
                     "parallel" = FALSE,
                     "numIslands" = 4,
                     "migrationRate" = 0.1,
                     "migrationInterval" = 10,
                     "optimMethod" = "L-BFGS-B",
                     "poptim" = 0.05,
                     "pressel" = 0.5,
                     "fnscale" = -1,
                     "maxit" = 100,
                     "initMclust" = "SVD",
                     "ModelNames" = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "EVV", "VVV","VEE"),
                     "G" = 1:9

  )


assign(".ppgmmga.default", restorealg, envir = asNamespace("ppgmmga"))
invisible(restorealg)

}



