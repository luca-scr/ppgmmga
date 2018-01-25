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
                       pro,
                       mean,
                       sigma,
                       method = c("UT", "VAR", "SOTE","MC"),
                       nSamples = 1e5)
{

  method <- match.arg(method, choices = eval(formals(EntropyGMM)$method))
  sigma <- as.matrix(sigma)
  d <- dim(sigma)[1]

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
                                       nSample = nSamples)
  )
  
  attributes(Entropy) <- list(names = names(Entropy),
                              approximation = method)
  return(Entropy)
}

############################################
# NENTROPY FOR THE GAUSSIAN MIXTURE MODELS #
############################################

NegentropyGMM <- function(G,
                          pro,
                          mean,
                          sigma,
                          sigmaGauss,
                          method = c("UT", "VAR", "SOTE","MC"),
                          nSamples = 1e5)
{
  sigma <- as.matrix(sigma)
  d <- dim(sigma)[1]
  if(d != nrow(sigmaGauss)) 
    { stop("The dimension of GMM and the Gaussian distribution must be equal") }
  method <- match.arg(method, choices = eval(formals(NegentropyGMM)$method))

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
                                       nSample = nSamples)
  )

  if(method != "MC")
  {
    Negentropy <- - Entropy[[1]] + EntropyGauss(S = sigmaGauss,d = d)
  } else
  {
    Negentropy <- Entropy
    Negentropy[[1]] <- -Entropy[[1]] + EntropyGauss(S = sigmaGauss,d = d)
  }

  attributes(Negentropy) <- list(names = names(Negentropy),
                                 approximation = method)
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
                      nSample)
{

  if(d == 1)
  {
    par <- list(pro = as.vector(pro),
                mean = as.matrix(mean),
                variance = list(modelName = "V",
                                d = d,
                                G = G,
                                sigmasq = sigma))
    data <- sim("V", parameters = par, n = nSample)[,-1,drop=FALSE]

  } else
  {
    par <- list(pro = pro,
                mean = mean,
                variance = list(modelName = "VVV",
                                d = d,
                                G = G,
                                sigma = sigma))
    data <- sim("VVV", parameters = par, n = nSample)[,-1]
  }

  b <- Sim(data = data,
           G = par$variance$G,
           pro = par$pro,
           mean = par$mean,
           sigma = par$variance$sigma,
           S = nSample)
  return(b)
}

###########################################
#         MONTE CARLO NENTROPY           #
###########################################

NegentropyMC <- function(par,
                         GMM,
                         p,
                         d,
                         nSample = 1e5,
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

  # simulate data from GMM
  # x <- sim("VVV", parameters = par, n = S)[,-1]
  # f <- array(NA, c(nrow(x), G))
  # for(j in 1:G)
  # { f[,j] <- mclust:::mvdnorm(x, Mu[j,], Sigma[,,j], log = TRUE) }
  # f <- sweep(f, 2, FUN = "+", STATS = log(par$pro))
  # logf <- apply(f, 1, mclust:::logsumexp)

  # not needed to simulate from N(0,I) because closed form solution
  # par0 <- list(pro = 1, mean = matrix(0, d, 1),
  #              variance = list(modelName = "XII", d = d, G = 1, sigmasq = 1,
  #                              Sigma = diag(d), sigma = array(diag(d), c(d,d,1)),
  #                              scale = 1))
  # x0 <- mclust::sim("XII", n = 1e5, parameters = par0)
  # logf0 <- mclust:::dmvnorm(x, c(0,0), diag(2), log = TRUE)
  # -mean(logf0) # = (0.5*(1+log(2*pi)))*d

  # Negentropy
  EMC <- EntropyMC(G = GMM$G,
                    pro = GMM$parameters$pro,
                    mean = transfGMM$mean,
                    sigma = transfGMM$sigma,
                    d = d,
                    nSample = nSample)

  z <- qnorm(1-(level/2))
  Negentropy <- - EMC$Entropy + EntropyGauss(S = transfGMM$sz, d = d)

  output <- list(Negentropy = Negentropy,
                 se = EMC$se,
                 confint = c(Negentropy-z*EMC$se, Negentropy+z*EMC$se))
  return(output)
}


###########################################
#        VOLUME OF THE DATA              #
###########################################

volume <- function(data, 
                   method = c("box", "pc", "ConvHull"))
{

  data <- as.matrix(data)
  dim <- dim(data)
  n <- dim[1]
  d <- dim[2]
  method <- match.arg(method, choices = eval(formals(volume)$method))

  # if((n > 5000 | d > 5) & method == "ConvHull" )
  # {
  #
  #   method = "box"
  #   warning("The Convex Hull method is not available for dataset with dimension bigger than 5 or number of observation bigger that 2000. The box method is used.")
  #
  # }

  sumlogdifcol <- function(x) 
    sum(log(apply(x, 2, function(colm) diff(range(colm)))))

  conv <- NULL
  V <- switch(method,
              "box" = { exp(sumlogdifcol(data)) },
              "pc" =  { exp(sumlogdifcol(princomp(data)$scores)) },
              "ConvHull" = { conv <- convhulln(data,options = "FA")
                             conv$vol }
             )

  #out <- list(volume = V, method = method, coord = conv$hull)
  attributes(V) <- list(method = method, coord = conv$hull)
  return(V)
}

