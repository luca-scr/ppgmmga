######################################################################################
#     FUNCTION RETURNS THE MATRIX BASIS GIVEN A ENCODED FORM OF THE BASIS            #
######################################################################################

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
# ENTROPY FOR THE GAUSSIAN MIXTURE MODELS #
###########################################

EntropyGMM <- function(G, 
                       d,
                       pro,
                       mean,
                       sigma,
                       method = c("UT", "VAR", "SOTE","MC"),
                       nsamples = 1e5)
{

  method <- match.arg(method, choices = eval(formals(EntropyGMM)$method))

  pro <- as.vector(pro)
  if(d == 1)
  {
    mean <- matrix(mean, ncol = G)
    sigma <- array(sigma, dim = c(d,d,G))
  }
  
  Entropy <- switch(method,
                    "UT" =  EntropyUT(G = G,
                                      pro = pro,
                                      mean = mean,
                                      sigma = sigma,
                                      d = d),
                    "VAR" =  EntropyVAR(G = G,
                                        pro = pro,
                                        mean = mean,
                                        sigma= sigma,
                                        d = d),
                    "SOTE" =  EntropySOTE(data = mean,
                                          G = G,
                                          pro = pro,
                                          mean = mean,
                                          sigma = sigma),
                    "MC" =   EntropyMC(G = G,
                                       pro = pro,
                                       mean = mean,
                                       sigma = sigma,
                                       d = d,
                                       nsamples = nsamples)
  )
  
  if(!is.list(Entropy))
    Entropy <- list(Entropy = Entropy)
  Entropy$method <- method
  return(Entropy)
}

##############################################
# NEGENTROPY FOR THE GAUSSIAN MIXTURE MODELS #
##############################################

NegentropyGMM <- function(G,
                          d,
                          pro,
                          mean,
                          sigma,
                          sigmaGauss,
                          method = c("UT", "VAR", "SOTE","MC"),
                          nsamples = 1e5)
{
  method <- match.arg(method, 
                      choices = eval(formals(NegentropyGMM)$method))
  if(d != nrow(sigmaGauss) | d != ncol(sigmaGauss))
    { stop("The dimension of sigmaGauss does not match that of GMM parameters!") }
  
  if(d == 1)
  {
    mean <- matrix(mean, ncol = G)
    sigma <- array(sigma, dim = c(d,d,G))
  }
    
  Entropy <- switch(method,
                    "UT" =  EntropyUT(G = G,
                                      pro = pro,
                                      mean = mean,
                                      sigma = sigma,
                                      d = d),
                    "VAR" =  EntropyVAR(G = G,
                                        pro = pro,
                                        mean = mean,
                                        sigma= sigma,
                                        d = d),
                    "SOTE" =  EntropySOTE(data = mean,
                                          G = G,
                                          pro = pro,
                                          mean = mean,
                                          sigma = sigma),
                    "MC" =   EntropyMC(G = G,
                                       pro = pro,
                                       mean = mean,
                                       sigma = sigma,
                                       d = d,
                                       nsamples = nsamples)
  )

  if(method != "MC")
  {
    Negentropy <- -Entropy[[1]] + EntropyGauss(S = sigmaGauss, d = d)
    se <- NA
  } else
  {
    Negentropy <- -Entropy[[1]] + EntropyGauss(S = sigmaGauss, d = d)
    se <- Entropy$se
  }

  Negentropy <- list(Negentropy = Negentropy,
                     se = se,
                     method = method)
  return(Negentropy)
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







