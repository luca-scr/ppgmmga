ppgmmga <- function(data,
                    d,
                    approx = c("UT","VAR","SOTE"),
                    center = TRUE,
                    scale = TRUE,
                    gmm = NULL,
                    gatype = c("ga", "gaisl"),
                    options = ppgmmga.options(),
                    seed = NULL,
                    verbose = interactive(),
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
  { 
    d <- p
    warning("dimension of projection subspace greater than the number of variables;\nautomatically set equal to ", d)
  }
  
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
                         verbose = verbose, 
                         ...)
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
                                          poptim = options$optimPoptim,
                                          pressel = options$optimPressel,
                                          control = list(fnscale = -1,
                                                         maxit = options$optimMaxit)),
                         selection = options$selection,
                         crossover = options$crossover,
                         mutation = options$mutation,
                         # elitism = base::max(1, round(options$popSize*0.05)),
                         elitism = 1,
                         pcrossover = options$pcrossover,
                         pmutation = options$pmutation,
                         parallel = options$parallel,
                         seed = seed, 
                         monitor = verbose,
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
                                                poptim = options$optimPoptim,
                                                pressel = options$optimPressel,
                                                control = list(fnscale = -1,
                                                               maxit = options$optimMaxit)),
                               selection = options$selection,
                               crossover = options$crossover,
                               mutation = options$mutation,
                               pcrossover = options$pcrossover,
                               pmutation = options$pmutation,
                               # elitism = base::max(1, round(options$popSize*0.05)),
                               elitism = 1,
                               numIslands = options$numIsland,
                               migrationRate = options$migrationRate,
                               migrationInterval = options$migrationInterval,
                               parallel = options$parallel,
                               seed = seed,
                               monitor = verbose,
                               ...)
               # "de" = de(type = "real-valued",
               #           fitness = function(...) 
               #             approxFun(...,
               #                       gmm = gmm,
               #                       p = p,
               #                       d = d,
               #                       decomposition = options$orth),
               #           popSize = options$popSize,
               #           maxiter = options$maxiter,
               #           run = options$run,
               #           optim = options$optim,
               #           optimArgs = list(method = options$optimMethod,
               #                            poptim = options$optimPoptim,
               #                            pressel = options$optimPressel,
               #                            control = list(fnscale = -1, 
               #                                           maxit = options$optimMaxit)),
               #           lower = rep(rep(0, p-1), d),
               #           upper = rep(c(2*pi, rep(pi,p-2)), d),
               #           pcrossover = options$pcrossover,
               #           pmutation = 0,
               #           stepsize = NA,
               #           parallel = options$parallel,
               #           seed = seed, 
               #           monitor = verbose,
               #           ...)
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
  dotargs <- list(...)
  if(is.null(dotargs$head)) dotargs$head <- 20
  if(is.null(dotargs$tail)) dotargs$tail <- 2
  if(is.null(dotargs$chead)) dotargs$chead <- 10
  if(is.null(dotargs$ctail)) dotargs$ctail <- 2
  
  # cat("-------------------------------------------------------","\n")
  # cat("Projection Pursuit GMM estimated via Genetic Algorithms","\n")
  # cat("-------------------------------------------------------","\n\n")
  cat(cli::rule(left = crayon::bold("ppgmmga"), 
                width = min(getOption("width"),40)), "\n\n")

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
  do.call(".printShortMatrix", 
          c(list(x$encodedBasis, digits = digits), 
            dotargs[c("head", "tail", "chead", "ctail")]))
  cat("\n")
  cat(paste("Estimated projection basis:","\n"))
  do.call(".printShortMatrix", 
          c(list(x$basis, digits = digits), 
            dotargs[c("head", "tail", "chead", "ctail")]))

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

#----------------------------------------------------------------------------#
# print a short version of a matrix by allowing to select the number of 
# head/tail rows and columns to display

.printShortMatrix <- function(x, head = 2, tail = 1, chead = 5, ctail = 1, ...)
{ 
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  if(is.na(head <- as.numeric(head))) head <- 2
  if(is.na(tail <- as.numeric(tail))) tail <- 1
  if(is.na(chead <- as.numeric(chead))) chead <- 5
  if(is.na(ctail <- as.numeric(ctail))) ctail <- 1
  
  if(nr > (head + tail + 1))
    { rnames <- rownames(x)
      if(is.null(rnames)) 
        rnames <- paste("[", 1:nr, ",]", sep ="")
      x <- rbind(x[1:head,,drop=FALSE], 
                 rep(NA, nc), 
                 x[(nr-tail+1):nr,,drop=FALSE])
      rownames(x) <- c(rnames[1:head], " ... ", rnames[(nr-tail+1):nr])
  }
  if(nc > (chead + ctail + 1))
    { cnames <- colnames(x)
      if(is.null(cnames)) 
        cnames <- paste("[,", 1:nc, "]", sep ="")
      x <- cbind(x[,1:chead,drop=FALSE], 
                 rep(NA, nrow(x)), 
                 x[,(nc-ctail+1):nc,drop=FALSE])
      colnames(x) <- c(cnames[1:chead], " ... ", cnames[(nc-ctail+1):nc])
  }
          
  print(x, na.print = "", ...)
}

