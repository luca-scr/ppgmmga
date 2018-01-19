ppgmmga <- function(data,
                    d,
                    approx = c("UT","VAR","SOTE"),
                    center = TRUE,
                    scale = FALSE,
	                  gmm = NULL,
                    gatype = c("ga", "gaisl"),
                    opt = list(),
                    monitor = TRUE,
                    ...)
{

  call <- match.call()

  if(center == FALSE & scale == FALSE)
    { warning("Input data are neither centred nor scaled.") }
  
  ### User can define its own fitness function (approximation)
  approxFun <- approx

  if(!is.function(approx))
  { 
    approx <- match.arg(approx, choices = eval(formals(ppgmmga)$approx))
    approxFun <- switch(approx,
                        "UT"   = NegentropyUT,
                        "SOTE" = NegentropySOTE,
                        "VAR"  = NegentropyVAR)
  }

  gatype <- match.arg(gatype, choices = eval(formals(ppgmmga)$gatype))

  data <-  data.matrix(data)
  p <- dim(data)[2]

  if(d > p)
    { stop(paste("d must be smaller or equal to", p)) }

  options <- ppgmmga.options(opt)

  ### center & scale the data
  Z <- scale(x = data, center = center, scale = scale)
  dimnames(Z) <- dimnames(data)

  ## GMM density estimation
  if(is.null(gmm))
  { gmm <- densityMclust(data = Z,
                         modelNames = options$ModelNames,
                         G = options$G,
                         initialization = list(hcPairs = hc(Z, use = options$initMclust)),verbose = FALSE,...)
  } else
  { if(!inherits(gmm, "densityMclust"))
    stop("If provided, argument 'gmm' must be an object of class 'densityMclust'.")
    if(!all(gmm$data == Z)) # (!identical(gmm$data, Z))
      stop("Input data appear to be different from that used for density estimation!")
  }

 ## GA negentropy maximization
  GA <- switch(gatype,
               "ga" = ga(type = "real-valued",
                         fitness = approxFun,
                         gmm = gmm,
                         d = d,
                         p = p,
                         popSize = options$popSize,
                         maxiter = options$maxiter,
                         run = options$run,
                         min = rep(rep(0, p-1), d),
                         max = rep(c(2*pi, rep(pi,p-2)), d),
                         optim = options$optim,
                         optimArgs = list(method = options$optimMethod,
                                          poptim = options$poptim,
                                          pressel = options$pressel,
                                          control = list(fnscale = options$fnscale, maxit = options$maxit)),
                         decomposition = options$Orth,
                         selection = options$selection,
                         crossover = options$crossover,
                         mutation = options$mutation,
                         elitism = base::max(1, round(options$popSize*0.05)),
                         pcrossover = options$pcr,
                         pmutation = options$pm,
                         parallel = options$parallel,
                         monitor = if(isTRUE(monitor)) {gaMonitor2},
                         ...),
               "gaisl" = gaisl(type = "real-valued",
                               fitness = approxFun,
                               gmm = gmm,
                               d = d,
                               p = p,
                               min = rep(rep(0, p-1), d),
                               max = rep(c(2*pi, rep(pi,p-2)), d),
                               popSize = options$popSize,
                               maxiter = options$maxiter,
                               run = options$run,
                               optim = options$optim,
                               optimArgs = list(method = options$optimMethod,
                                                poptim = options$poptim,
                                                pressel = options$pressel,
                                                control = list(fnscale = options$fnscale, maxit = options$maxit)),
                               decomposition = options$Orth,
                               selection = options$selection,
                               crossover = options$crossover,
                               mutation = options$mutation,
                               pcrossover = options$pcr,
                               pmutation = options$pm,
                               elitism = base::max(1, round(options$popSize*0.05)),
                               numIslands = options$numIsland,
                               migrationRate = options$migrationRate,
                               migrationInterval = options$migrationInterval,
                               parallel = options$parallel,
                               monitor = if(isTRUE(monitor)) {gaislMonitor2},
                               ...)
  )

  B <- encodeBasis(GA@solution[1,], p = p, d = d,
                       decomposition = options$Orth)

  colnames(B) <- paste0("PP", seq(d))
  rownames(B) <- colnames(data)

  out <- list(call = call,
              data = Z,
              d = d,
              approx = approx,
              GMM = gmm,
              GA = GA,
              Negentropy = GA@fitnessValue,
              basis = B,
              Z = Z %*% B)
  class(out) <- "ppgmmga"

  return(out)
}


print.ppgmmga <- function(x, ...)
{
  if(!is.null(cl <- x$call))
    { cat("Call:\n")
      dput(cl, control = NULL)
    }
  cat("\n'ppgmmga' object containing:","\n")
  print(names(x)[-1])
  invisible()
}

summary.ppgmmga <- function(object,
                            ppgmmgaCheck = FALSE,
                            ...)
{
  out <- list(approx = object$approx,
              Negentropy = object$Negentropy,
              basis = object$basis,
              p = object$GMM$d,
              d = object$d,
              n = nrow(object$Z),
              encodedBasis = object$GA@solution,
              modelName = object$GMM$modelName,
              G = object$GMM$G,
              data = object$data)

  if(ppgmmgaCheck)
  {

    MC <- ppgmmgaCheck(object, ...)
    out <- append(out, list(check = MC))

  }

  class(out) <- "summary.ppgmmga"
  return(out)
}

print.summary.ppgmmga <- function(object, digits = getOption("digits"), ...)
{


  cat("-------------------------------------------------------","\n")
  cat("Projection Pursuit GMM estimated via Genetic Algorithms","\n")
  cat("-------------------------------------------------------","\n\n")

  cat(paste("Data dimensions               =",
            object$n, "x", object$p, "\n"))
  x <- names(attributes(object$data))[names(attributes(object$data)) == "scaled:center" | names(attributes(object$data)) == "scaled:scale"]
  if(length(x) == 0)
  {
  cat(paste("Data transformation           = No preliminary transformation choosen in the argument","\n"))

  }else if (length(x)==2){

  cat(paste("Data transformation           = Center and Scale","\n"))

  }else if (x == "scaled:center"){

  cat(paste("Data transformation           = Center","\n"))

  }else{

  cat(paste("Data transformation           = Scale","\n"))

  }

  cat(paste("Projection subspace dimension =",
            object$d, "\n"))
  cat(paste0("GMM density estimate          = (",
             object$modelName,",", object$G,")", "\n"))
  cat(paste("Negentropy approximation      =",
            object$approx, "\n"))
  cat(paste("GA optimal negentropy         =",
            signif(object$Negentropy, digits), "\n"))

  cat(paste("GA encoded basis solution:","\n"))
  print(object$encodedBasis, digits = digits)
  cat("\n")
  cat(paste("Estimated projection basis:","\n"))
  print(object$basis, digits = digits)
  cat("\n")
  if(any(names(object)=="check"))
  {
    cat(paste("Monte Carlo Negentropy approximation with estimated basis:","\n"))
    print(t(object$check), digits = digits)
  }

  invisible()
}
