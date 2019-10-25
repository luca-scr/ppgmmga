##
## Matrix basis given an encoded form of the basis
##

encodeBasis <- function(par, p, d, orth = TRUE, decomp = c("QR","SVD"))
{
# Reference: Baragona, Battaglia, Poli (2011, p. 78)
  
  decomp <- match.arg(decomp, 
                      choices = eval(formals(encodeBasis)$decomp))
  B <- encodebasis(par = par, d = d, p = p)
  if(orth) 
    B <- orth(B, method = decomp)
  return(B)
}

##
## Exact GMM-based entropy
##

# For the Rcpp version see EntropyGMM()
# this is the companion R version (which is less efficient)
EntropyGMM_R <- function(data, pro, mean, sigma, ...)
{
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  pro <- as.vector(pro)
  G <- length(pro)
  mean <- array(mean, dim = c(d, G))
  par <- list(pro = pro, mean = mean, variance = NULL)
  if(d == 1) 
  {
    modelName <- "V"
    par$variance <- mclustVariance(modelName, d = d, G = G)
    par$variance$sigmasq <- as.vector(sigma)
  } else
  {
    modelName <- "VVV"
    par$variance <- mclustVariance(modelName, d = d, G = G)
    sigma <- array(sigma, dim = c(d,d,G))
    par$variance$sigma <- sigma
    cholsigma <- array(as.double(NA), dim(sigma))
    for(k in 1:G)
      cholsigma[,,k] <- chol(sigma[,,k])
    par$variance$cholsigma <- cholsigma
  }
  # component densities
  cden <- cdens(modelName = modelName, parameters = par, 
                data = data, logarithm = TRUE)
  cden <- cbind(cden) # drop redundant attributes
  cden <- sweep(cden, 2, FUN = "+", STATS = log(pro))
  # posterior probs
  z <- exp(sweep(cden, MARGIN = 1, FUN = "-", 
                 STATS = apply(cden, 1, logsumexp)))
  # log-densities
  maxlog <- apply(cden, 1, max)
  cden <- sweep(cden, 1, FUN = "-", STATS = maxlog)
  logdens <- log(apply(exp(cden), 1, sum)) + maxlog
  # entropy
  h <- -sum(sweep(z, MARGIN = 1,  FUN = "*", STATS = logdens))/n
    
  return(h)
}

##
## Exact GMM-based negentropy
##

NegentropyGMM <- function(par, GMM, p, d, decomp)
{

  # Encode orthonormal basis
  B <- encodeBasis(par = par, d = d, p = p, orth = TRUE, decomp = decomp)
  # Transform estimated parameters to projection subspace
  transfGMM <-  LinTransf(mean = GMM$parameters$mean,
                          sigma = GMM$parameters$variance$sigma,
                          B = B,
                          Z = GMM$data,
                          G = GMM$G,
                          d = d)
  # Negentropy
  Negentropy <-  -EntropyGMM(data = GMM$data %*% B,
                             pro = GMM$parameters$pro,
                             mean = transfGMM$mean,
                             sigma = transfGMM$sigma) + 
                  EntropyGauss(S = transfGMM$sz, d = d)

  return(Negentropy)
}


##
## Monte Carlo entropy
##

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

##
## Monte Carlo negentropy
##

NegentropyMC <- function(par,
                         GMM,
                         p,
                         d,
                         nsamples = 1e5,
                         level = 0.05,
                         decomp = "QR")
{
  B <- encodeBasis(par, p = p, d = d, decomp = decomp)
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

##
## Unscented Trasformation (UT) approximation negentropy
##

NegentropyUT <- function(par, GMM, p, d, decomp)
{

  # Encode orthonormal basis
  B <- encodeBasis(par = par, d = d, p = p, 
                   orth = TRUE, decomp = decomp)
  # Transform estimated parameters to projection subspace
  transfGMM <-  LinTransf(mean = GMM$parameters$mean,
                          sigma = GMM$parameters$variance$sigma,
                          B = B,
                          Z = GMM$data,
                          G = GMM$G,
                          d = d)
  # Negentropy
  Negentropy <-  -EntropyUT(G = GMM$G,
                            pro = GMM$parameters$pro,
                            mean = transfGMM$mean ,
                            sigma = transfGMM$sigma,
                            d = d) + 
                  EntropyGauss(S = transfGMM$sz, d = d)

  return(Negentropy)
}

##
## VARiational (VAR) approximation negentropy
##

NegentropyVAR <- function(par, GMM, p, d, decomp)
{

  # Encode orthonormal basis
  B <- encodeBasis(par = par, d = d, p = p, 
                   orth = TRUE, decomp = decomp)
  # Transform estimated parameters to projection subspace
  transfGMM <-  LinTransf(mean = GMM$parameters$mean,
                          sigma = GMM$parameters$variance$sigma,
                          B = B,
                          Z = GMM$data,
                          G = GMM$G,
                          d = d)
  # Negentropy
  Negentropy <- -EntropyVAR(G = GMM$G,
                            pro = GMM$parameters$pro,
                            mean = transfGMM$mean,
                            sigma= transfGMM$sigma,
                            d = d) + 
                 EntropyGauss(S = transfGMM$sz, d = d)

  return(Negentropy)
}

##
## Second Order Taylor Expansion (SOTE) approximation negentropy
##

NegentropySOTE <- function(par, GMM, d, p, decomp)
{
  # Encode orthonormal basis
  B <- encodeBasis(par = par, d = d, p = p, 
                   orth = TRUE, decomp = decomp)
  # Transform estimated parameters to projection subspace
  transfGMM <-  LinTransf(mean = GMM$parameters$mean,
                          sigma = GMM$parameters$variance$sigma,
                          B = B,
                          Z = GMM$data,
                          G = GMM$G,
                          d = d)
  # Negentropy
  Negentropy <- -EntropySOTE(data = transfGMM$mean,
                         G = GMM$G,
                         pro = GMM$parameter$pro,
                         mean = transfGMM$mean,
                         sigma = transfGMM$sigma) +
                EntropyGauss(S = transfGMM$sz, d = d)

  return(Negentropy)
}
