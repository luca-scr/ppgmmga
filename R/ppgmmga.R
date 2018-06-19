ppgmmga <- function(data,
                    d,
                    approx = c("UT","VAR","SOTE"),
                    center = TRUE,
                    scale = FALSE,
                    gmm = NULL,
                    gatype = c("ga", "gaisl", "de"),
                    opt = list(),
                    seed = NULL,
                    ...)
{

  call <- match.call()

  center <- as.logical(center)
  scale  <- as.logical(scale)
  if(!center & !scale)
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
  { 
    gmm <- densityMclust(data = Z,
                         modelNames = options$modelNames,
                         G = options$G,
                         initialization = list(hcPairs = hc(Z, use = options$initMclust)),
                         verbose = FALSE, ...)
  } else
  { 
    if(!inherits(gmm, "densityMclust"))
      stop("If provided, argument 'gmm' must be an object of class 'densityMclust'.")
    if(any((gmm$data - Z) > sqrt(.Machine$double.eps)))
      stop("Input data appears to be different from data used for density estimation!")
  }

  ## GA negentropy maximization
  GA <- switch(gatype,
               "ga" = ga(type = "real-valued",
                         fitness = approxFun,
                         gmm = gmm,
                         p = p,
                         d = d,
                         decomposition = options$orth,
                         lower = rep(rep(0, p-1), d),
                         upper = rep(c(2*pi, rep(pi,p-2)), d),
                         popSize = options$popSize,
                         maxiter = options$maxiter,
                         run = options$run,
                         optim = options$optim,
                         optimArgs = list(method = options$optimMethod,
                                          poptim = options$poptim,
                                          pressel = options$pressel,
                                          control = list(fnscale = options$fnscale, maxit = options$maxit)),
                         selection = options$selection,
                         crossover = options$crossover,
                         mutation = options$mutation,
                         elitism = base::max(1, round(options$popSize*0.05)),
                         pcrossover = options$pcr,
                         pmutation = options$pm,
                         parallel = options$parallel,
                         seed = seed, 
                         ...),
               "gaisl" = gaisl(type = "real-valued",
                               fitness = approxFun,
                               gmm = gmm,
                               p = p,
                               d = d,
                               decomposition = options$orth,
                               lower = rep(rep(0, p-1), d),
                               upper = rep(c(2*pi, rep(pi,p-2)), d),
                               popSize = options$popSize,
                               maxiter = options$maxiter,
                               run = options$run,
                               optim = options$optim,
                               optimArgs = list(method = options$optimMethod,
                                                poptim = options$poptim,
                                                pressel = options$pressel,
                                                control = list(fnscale = options$fnscale, maxit = options$maxit)),
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
                               seed = seed,
                               ...),
               "de" = de(type = "real-valued",
                         fitness = function(...) 
                           approxFun(...,
                                     gmm = gmm,
                                     p = p,
                                     d = d,
                                     decomposition = options$orth),
                         popSize = options$popSize,
                         maxiter = options$maxiter,
                         run = options$run,
                         optim = options$optim,
                         optimArgs = list(method = options$optimMethod,
                                          poptim = options$poptim,
                                          pressel = options$pressel,
                                          control = list(fnscale = options$fnscale, maxit = options$maxit)),
                         lower = rep(rep(0, p-1), d),
                         upper = rep(c(2*pi, rep(pi,p-2)), d),
                         pcrossover = options$pcr,
                         pmutation = 0,
                         stepsize = NA,
                         parallel = options$parallel,
                         seed = seed, 
                         ...)
  )

  B <- encodeBasis(GA@solution[1,], p = p, d = d,
                   decomposition = options$orth)
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
  { 
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  cat("\n'ppgmmga' object containing:","\n")
  print(names(x)[-1])
  invisible()
}

summary.ppgmmga <- function(object, check = FALSE, ...)
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

  if(check)
    { out$check <- NegentropyMCcheck(object, ...) }

  class(out) <- "summary.ppgmmga"
  return(out)
}

print.summary.ppgmmga <- function(x, digits = getOption("digits"), ...)
{
  cat("-------------------------------------------------------","\n")
  cat("Projection Pursuit GMM estimated via Genetic Algorithms","\n")
  cat("-------------------------------------------------------","\n\n")

  cat(paste("Data dimensions               =", x$n, "x", x$p, "\n"))
  cat(paste("Data transformation           = "))
  trasf <- names(attributes(x$data))
  trasf <- ifelse(any(trasf == "scaled:center"), 1, 0) +
           ifelse(any(trasf == "scaled:scale"),  2, 0)
  cat(switch(as.character(trasf), 
             "0" = "none", "1" = "center", 
             "2" = "scale", "3" = "center & scale"), "\n")

  cat(paste("Projection subspace dimension =", 
            x$d, "\n"))
  cat(paste0("GMM density estimate          = (", 
             x$modelName,",", x$G,")", "\n"))
  cat(paste("Negentropy approximation      =",
            x$approx, "\n"))
  cat(paste("GA optimal negentropy         =",
            signif(x$Negentropy, digits), "\n"))

  cat(paste("GA encoded basis solution:","\n"))
  print(x$encodedBasis, digits = digits)
  cat("\n")
  cat(paste("Estimated projection basis:","\n"))
  print(x$basis, digits = digits)
  
  if(!is.null(x$check))
  {
    cat(paste("\nMonte Carlo Negentropy approximation check:","\n"))
    tab <- x$check
    tab <- t(data.frame(tab[c(2:4,6)], check.names = FALSE))
    colnames(tab) <- x$check$approx
    print(tab, digits = digits)
  }

  invisible()
}

NegentropyMCcheck <- function(object, nsamples = 1e5, conf.level = NULL)
{
  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")
  
  MC <- NegentropyMC(par = object$GA@solution,
                     GMM = object$GMM,
                     p = object$GMM$d,
                     d = object$d,
                     nsamples = nsamples,
                     level = 1-conf.level)

  out <- list("approx"            = object$approx,
              "Approx Negentropy" = object$Negentropy,
              "MC Negentropy"     = MC$Negentropy, 
              "MC se"             = MC$se,
              "MC interval"       = MC$confint,
              "Relative accuracy" = object$Negentropy/MC$Negentropy)
  return(out)
}
