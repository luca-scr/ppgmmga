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



###########################################
#      MONTE CARLO NENTROPY FOR ica       #
###########################################

NegentropyICA <- function(object, d = object$d, 
                          NegApprox = c("logcosh", "exp", "kur"), nsamples = 1e5)
{
  require(ica)
  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")
  
  fun <- match.arg(NegApprox, choices = eval(formals(NegentropyICA)$NegApprox))
  
  ICA <- icafast(X = object$data, nc = d, center = FALSE,fun = fun)
  B <- as.matrix(ICA$W)[,seq(d),drop=FALSE]
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
  
  ICA$Negentropy <- Negentropy
  ICA$se <- EMC$se
  return(ICA)
}


###########################################
#      MONTE CARLO NENTROPY FOR fastICA       #
###########################################

NegentropyfastICA <- function(object, d = object$d, 
                          NegApprox = c("logcosh", "exp"), nsamples = 1e5)
{
  require(fastICA)
  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")
  
  fun <- match.arg(NegApprox, choices = eval(formals(NegentropyfastICA)$NegApprox))
  
  fastICA <- fastICA(X = object$data, n.comp = d, fun = fun, method = "C")
  B <- as.matrix(fastICA$K)[,seq(d),drop=FALSE]
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
  
  fastICA$Negentropy <- Negentropy
  fastICA$se <- EMC$se
  return(fastICA)
}



###########################################
# MONTE CARLO NENTROPY FOR factor  analysis #
###########################################

NegentropyFA <- function(object, d = object$d, 
                              rotation = "varimax", nsamples = 1e5)
{

    if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")
  

  FA <- factanal(object$data, factors = d, rotation = rotation)
  B <- as.matrix(FA$loadings)[,seq(d),drop=FALSE]
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
  
  FA$Negentropy <- Negentropy
  FA$se <- EMC$se
  FA$Z <- object$data %*% B
  return(FA)
}