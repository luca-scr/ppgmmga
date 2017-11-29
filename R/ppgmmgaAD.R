
ppgmmgaAD <- function(data,
                       G = NULL,
                       modelNames = NULL,
                       volumeMethod = c("box", "pc","ConvHull"),
                       d = 2,
                       approx = c("UT", "VAR", "SOTE"),
                       center = TRUE,
                       scale = FALSE,
                       gmm = NULL,
                       gatype = c("ga","gaisl"),
                       nSub = NULL,
                       cutoff = 0.975,
                       noiseInit = NULL,
                       alpha = 0.5,
                       projData = NULL,
                       basis = NULL,
                       opt = list(),
                       ...)

{

  call <- match.call()

  volumeMethod <- match.arg(volumeMethod, choices = eval(formals(ppgmmgaAD)$volumeMethod))
  method <- match.arg(method, choices = eval(formals(ppgmmgaAD)$method))
  gatype <- match.arg(gatype, choices = eval(formals(ppgmmgaAD)$gatype))
  approx <- match.arg(approx, choices = eval(formals(ppgmmgaAD)$approx))

  dimension <- dim(data)


  out <- ppgmmga.out(data = data,
                     scale = scale,
                     center = center,
                     approx = approx,
                     volumeMethod = volumeMethod,
                     G = G,
                     modelNames = modelNames,
                     cutoff = cutoff,
                     d = d,
                     noiseInit = noiseInit,
                     alpha = alpha,
                     gmm = NULL,
                     gatype = gatype,
                     opt = opt,
                     projData = NULL,
                     basis = NULL,
                     ...)




  out$call <- call
  class(out) <- "ppgmmgaAD"
  return(out)

}





ppgmmga.out <- function(data,
                        approx,
                         d,
                         center,
                         scale,
                         volumeMethod,
                         G,
                         modelNames,
                         noiseInit,
                         cutoff,
                         alpha,
                        gmm,
                        gatype,
                        opt,
                        projData = NULL,
                        basis = NULL,
                         ...)

{



  data <- data.matrix(data)
  n <- nrow(data)


  if(is.null(projData)){

    pp <- ppgmmga(data = data, d = d,approx = approx ,scale = scale,opt = opt, gmm = gmm, gatype = gatype, ...)
    Z <- pp$Z

  }else{

    Z <- projData

  }

  if(!is.null(volumeMethod)){

    vol <- volume(Z,
                method = volumeMethod)
    V <- vol$volume
    attributes(V) <- list(method = vol$method, coord = vol$coord)


  }else{

    V <- NULL
  }

  modNoise <- Mclust.Noise(data = pp$Z,
                           vol = V,
                           G = G,
                           modelNames = modelNames,
                           noiseInit = noiseInit,
                           cutoff = cutoff,
                           alpha = alpha,...)

  pp$GMM <- modNoise


  return(append(pp,list(Class = modNoise$classification)))

}

  Noise <- function(data,
                    n,
                    cutoff = 0.975,
                    alpha)
  {

    robust <- CovMcd(x = data, alpha = alpha)
    cut <- sqrt(qchisq(cutoff, df = ncol(data)))
    if(is.null(cutoff))
    {
      cutoff <- 0.975
      cut <- sqrt(qchisq(cutoff, df = ncol(data)))
    }

    noise <- sqrt(robust@mah) > cut


    if(sum(noise) == n)
    {

      noise <- sample(c(TRUE,FALSE),size = n,replace = TRUE)

    }

    return(noise)


  }




Mclust.Noise <- function(data,
                         G = NULL,
                         modelNames = NULL,
                         vol = NULL,
                         noiseInit = NULL,
                         cutoff = 0.975,
                         alpha = 0.5,
                         ...)
{


  data <- data.matrix(data)
  n <- nrow(data)


  if(is.null(noiseInit)){

    noiseInit <- Noise(data = data,
                       n = n,
                       cutoff = cutoff,
                       alpha = alpha)

  }


  if(is.null(vol)){

    modNoise <-   tryCatch(Mclust(data, G = G, modelNames = modelNames,initialization = list(noise = noiseInit),...), warning = function(w) {invisible(w)},silent  =  TRUE)

  }else{

    modNoise <- tryCatch(Mclust(data,
                                initialization = list(noise = noiseInit),
                                Vinv = 1/vol, G = G,
                                modelNames = modelNames, ...), warning = function(w) {invisible(w)},silent  =  TRUE)

  }

  if(inherits(modNoise,"warning"))
  {

    out <- list(NULL,classification = NA)
    warning("Model does not converge, try to increse the cutoff.")

  }else{out <- modNoise}

  return(out)

}






print.ppgmmgaAD <- function(x, ...)
{
  if(!is.null(cl <- x$call))
  { cat("Call:\n")
    dput(cl, control = NULL)
  }
  cat("\n'ppgmmgaAD' object containing:","\n")
  print(names(x)[-1])
  invisible()
}



summary.ppgmmgaAD <- function(object)
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
              data = object$data,
              classification = object$Class,
              BIC = object$GMM$bic,
              loglik = object$GMM$loglik,
              volume = object$GMM$hypvol)


  class(out) <- "summary.ppgmmgaAD"
  return(out)
}

print.summary.ppgmmgaAD <- function(object, digits = getOption("digits"), ...)
{


  cat("-------------------------------------------------------","\n")
  cat("Projection Pursuit GMM estimated via Genetic Algorithms","\n")
  cat("-------------------------------------------------------","\n\n")

  cat(paste("Data dimensions               =",
            object$n, "x", nrow(object$basis), "\n"))
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

  cat(paste("Negentropy approximation      =",
            object$approx, "\n"))
  cat(paste("GA optimal negentropy         =",
            signif(object$Negentropy, digits), "\n"))

  cat(paste0("GMM noise density estimate    = (",
             object$modelName,",", object$G,")", "\n"))

  cat(paste0("Volume Method 	              = ",
             attributes(object$volume)$method,"\n"))

  cat(paste0("Volume        	              = ",
             signif(object$volume,digits),"\n"))

  cat(paste0("BIC 	                     = ",
             signif(object$BIC, digits) ,"\n"))

  cat(paste0("log-likelihood 	              = ",
             signif(object$loglik, digits),"\n"))

  cat("\n")

  cat(paste("Classification:","\n"))
  print(table(object$classification))

  cat("\n")

  cat(paste("Estimated projection basis:","\n"))
  print(object$basis, digits = digits)
  cat("\n")



  invisible()
}


