###########################################
#      MONTE CARLO NENTROPY FOR PCA       #
###########################################

NegentropyPCA <- function(object, d = object$d, nsamples = 1e5)
{
  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")
  
  PCA <- prcomp(object$data)
  B <- as.matrix(PCA$rotation)[,seq(d),drop=FALSE]
  transfGMM <-  LinTransf(mean = object$GMM$parameters$mean,
                          sigma = object$GMM$parameters$variance$sigma,
                          B = B,
                          Z = object$GMM$data,
                          G = object$GMM$G,
                          d = d)
  EMC <- EntropyMC(G = object$GMM$G, 
                   pro = object$GMM$parameters$pro,
                   mean = transfGMM$mean,
                   sigma = transfGMM$sigma,
                   d = d,
                   nsamples = nsamples)
  Negentropy <- -EMC$Entropy + EntropyGauss(S = transfGMM$sz, d = d)
  
  PCA$Negentropy <- Negentropy
  PCA$se <- EMC$se
  return(PCA)
}