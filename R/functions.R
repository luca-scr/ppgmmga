###################################################
# MATRIX BASIS GIVEN AN ENCODED FORM OF THE BASIS #
###################################################

encodeBasis <- function(par, p, d, orth = TRUE, decomposition = c("QR","SVD"))
{
  decomposition <- match.arg(decomposition, 
                             choices = eval(formals(encodeBasis)$decomposition))
  B <- encodebasis(par = par, d = d, p = p)
  if(orth) 
    B <- orth(B, method = decomposition)
  return(B)
}


###########################################
#         MONTE CARLO ENTROPY             #
###########################################

EntropyMC <- function(G,
                      pro,
                      mean,
                      sigma,
                      d,
                      nsamples)
{
  if(d == 1)
  {
    par <- list(pro = as.vector(pro),
                mean = as.matrix(mean),
                variance = list(modelName = "V",
                                d = d,
                                G = G,
                                sigmasq = sigma))
    simdata <- sim("V", parameters = par, n = nsamples)[,-1,drop=FALSE]
  } else
  {
    par <- list(pro = pro,
                mean = mean,
                variance = list(modelName = "VVV",
                                d = d,
                                G = G,
                                sigma = sigma))
    simdata <- sim("VVV", parameters = par, n = nsamples)[,-1]
  }

  h <- EntropyMCapprox(data = simdata,
                       G = par$variance$G,
                       pro = par$pro,
                       mean = par$mean,
                       sigma = par$variance$sigma)
  return(h)
}

###########################################
#        MONTE CARLO NEGENTROPY           #
###########################################

NegentropyMC <- function(par,
                         GMM,
                         p,
                         d,
                         nsamples = 1e5,
                         level = 0.05,
                         decomposition = "QR")
{
  B <- encodeBasis(par, p = p, d = d, decomposition = decomposition)
  # Transform estimated parameters to projection subspace
  transfGMM <-  LinTransf(mean = GMM$parameters$mean,
                          sigma = GMM$parameters$variance$sigma,
                          B = B,
                          Z = GMM$data,
                          G = GMM$G,
                          d = d)
  # Negentropy
  EMC <- EntropyMC(G = GMM$G,
                   pro = GMM$parameters$pro,
                   mean = transfGMM$mean,
                   sigma = transfGMM$sigma,
                   d = d,
                   nsamples = nsamples)

  z <- qnorm(1-(level/2))
  Negentropy <- - EMC$Entropy + EntropyGauss(S = transfGMM$sz, d = d)

  output <- list(Negentropy = Negentropy,
                 se = EMC$se,
                 confint = c(Negentropy-z*EMC$se, Negentropy+z*EMC$se))
  return(output)
}

NegentropyMCcheck <- function(object, nsamples = 1e5, conf.level = NULL)
{
  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")
  
  MC <- NegentropyMC(par = object$GA@solution[1,],
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
