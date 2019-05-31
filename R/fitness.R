
NegentropyUT <- function(par, gmm, p, d, decomposition)
{

  # Encode orthonormal basis
  B <- encodeBasis(par = par, d = d, p = p, 
                   orth = TRUE, decomposition = decomposition)
  # Transform estimated parameters to projection subspace
  transfGMM <-  LinTransf(mean = gmm$parameters$mean,
                          sigma = gmm$parameters$variance$sigma,
                          B = B,
                          Z = gmm$data,
                          G = gmm$G,
                          d = d)
  # Negentropy
  Negentropy <-  -EntropyUT(G = gmm$G,
                            pro = gmm$parameters$pro,
                            mean = transfGMM$mean ,
                            sigma = transfGMM$sigma,
                            d = d) + 
                  EntropyGauss(S = transfGMM$sz, d = d)

  return(Negentropy)
}


NegentropyVAR <- function(par, gmm, p, d, decomposition)
{

  # Encode orthonormal basis
  B <- encodeBasis(par = par, d = d, p = p, 
                   orth = TRUE, decomposition = decomposition)
  # Transform estimated parameters to projection subspace
  transfGMM <-  LinTransf(mean = gmm$parameters$mean,
                          sigma = gmm$parameters$variance$sigma,
                          B = B,
                          Z = gmm$data,
                          G = gmm$G,
                          d = d)
  # Negentropy
  Negentropy <- -EntropyVAR(G = gmm$G,
                            pro = gmm$parameters$pro,
                            mean = transfGMM$mean,
                            sigma= transfGMM$sigma,
                            d = d) + 
                 EntropyGauss(S = transfGMM$sz, d = d)

  return(Negentropy)
}


NegentropySOTE <- function(par, gmm, d, p, decomposition)
{
  # Encode orthonormal basis
  B <- encodeBasis(par = par, d = d, p = p, 
                   orth = TRUE, decomposition = decomposition)
  # Transform estimated parameters to projection subspace
  transfGMM <-  LinTransf(mean = gmm$parameters$mean,
                          sigma = gmm$parameters$variance$sigma,
                          B = B,
                          Z = gmm$data,
                          G = gmm$G,
                          d = d)
  # Negentropy
  Negentropy <- -EntropySOTE(data = transfGMM$mean,
                         G = gmm$G,
                         pro = gmm$parameter$pro,
                         mean = transfGMM$mean,
                         sigma = transfGMM$sigma) +
                EntropyGauss(S = transfGMM$sz, d = d)

  return(Negentropy)
}
