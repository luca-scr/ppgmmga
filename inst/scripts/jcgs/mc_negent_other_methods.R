###########################################
#     MONTE CARLO NEGENTROPY FOR PCA      #
###########################################

NegentropyPCA <- function(object, 
                          d = object$d, 
                          nsamples = 1e5)
{
  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")

  PCA <- prcomp(object$data)
  B <- as.matrix(PCA$rotation)[,seq(d),drop=FALSE]
  transfGMM <-  ppgmmga:::LinTransf(mean = object$GMM$parameters$mean,
                          sigma = object$GMM$parameters$variance$sigma,
                          B = B,
                          Z = object$GMM$data,
                          G = object$GMM$G,
                          d = d)
  EMC <- ppgmmga:::EntropyMC(G = object$GMM$G, 
                   pro = object$GMM$parameters$pro,
                   mean = transfGMM$mean,
                   sigma = transfGMM$sigma,
                   d = d,
                   nsamples = nsamples)
  Negentropy <- -EMC$Entropy + ppgmmga:::EntropyGauss(S = transfGMM$sz, d = d)
  
  PCA$basis <- B
  PCA$Z <- PCA$x[,seq(d),drop=FALSE]
  PCA$Negentropy <- Negentropy
  PCA$se <- EMC$se
  return(PCA)
}



###########################################
#      MONTE CARLO NEGENTROPY FOR ICA     #
###########################################

NegentropyICA <- function(object, 
                          d = object$d, 
                          NegApprox = c("kur", "logcosh", "exp"),
                          alg = c("par", "def"),
                          nsamples = 1e5)
{
  stopifnot(require(ica))
  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")
  
  fun <- match.arg(NegApprox, choices = eval(formals(NegentropyICA)$NegApprox))
  alg <- match.arg(alg, choices = eval(formals(NegentropyICA)$alg))
  
  ICA <- icafast(X = object$data, nc = d, fun = fun, alg = alg)
  B <- t(as.matrix(ICA$W))
  B <- apply(B, 2, normalize)[,seq(d),drop=FALSE]
  transfGMM <-  ppgmmga:::LinTransf(mean = object$GMM$parameters$mean,
                          sigma = object$GMM$parameters$variance$sigma,
                          B = B,
                          Z = object$GMM$data,
                          G = object$GMM$G,
                          d = d)
  EMC <- ppgmmga:::EntropyMC(G = object$GMM$G, 
                   pro = object$GMM$parameters$pro,
                   mean = transfGMM$mean,
                   sigma = transfGMM$sigma,
                   d = d,
                   nsamples = nsamples)
  Negentropy <- -EMC$Entropy + ppgmmga:::EntropyGauss(S = transfGMM$sz, d = d)
  
  ICA$basis <- B
  ICA$Z <- ICA$S; 
  colnames(B) <- colnames(ICA$Z) <- paste0("ICA", seq(d))
  ICA$Negentropy <- Negentropy
  ICA$se <- EMC$se
  return(ICA)
}


###########################################
#   MONTE CARLO NEGENTROPY FOR FASTICA    #
###########################################

NegentropyFASTICA <- function(object, 
                              d = object$d, 
                              NegApprox = c("logcosh", "exp"), 
                              alg.typ = c("parallel", "deflation"),
                              nsamples = 1e5)
{
  stopifnot(require(fastICA))
  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")
  
  fun     <- match.arg(NegApprox, 
                       choices = eval(formals(NegentropyFASTICA)$NegApprox))
  alg.typ <- match.arg(alg.typ, 
                       choices = eval(formals(NegentropyFASTICA)$alg.typ))
  
  fastICA <- fastICA(X = object$data, n.comp = d, method = "C",
                     fun = fun, alg.typ = alg.typ)
  B <- as.matrix(fastICA$K %*% fastICA$W)
  B <- apply(B, 2, normalize)[,seq(d),drop=FALSE]
  transfGMM <-  ppgmmga:::LinTransf(mean = object$GMM$parameters$mean,
                          sigma = object$GMM$parameters$variance$sigma,
                          B = B,
                          Z = object$GMM$data,
                          G = object$GMM$G,
                          d = d)
  EMC <- ppgmmga:::EntropyMC(G = object$GMM$G, 
                   pro = object$GMM$parameters$pro,
                   mean = transfGMM$mean,
                   sigma = transfGMM$sigma,
                   d = d,
                   nsamples = nsamples)
  Negentropy <- -EMC$Entropy + ppgmmga:::EntropyGauss(S = transfGMM$sz, d = d)
  
  fastICA$basis <- B
  fastICA$Z <- fastICA$S
  colnames(fastICA$basis) <- colnames(fastICA$Z) <- paste0("ICA", seq(d))
  fastICA$Negentropy <- Negentropy
  fastICA$se <- EMC$se
  return(fastICA)
}



##############################################
#  MONTE CARLO NENTROPY FOR FACTOR ANALYSIS  #
##############################################

NegentropyFA <- function(object, 
                         d = object$d, 
                         rotation = "varimax", 
                         nsamples = 1e5)
{

  if(!inherits(object, "ppgmmga"))
    stop("'object' must be of class 'ppgmmga'")

  FA <- factanal(object$data, factors = d, rotation = rotation)
  B <- as.matrix(FA$loadings)[,seq(d),drop=FALSE]
  transfGMM <-  ppgmmga:::LinTransf(mean = object$GMM$parameters$mean,
                          sigma = object$GMM$parameters$variance$sigma,
                          B = B,
                          Z = object$GMM$data,
                          G = object$GMM$G,
                          d = d)
  EMC <- ppgmmga:::EntropyMC(G = object$GMM$G, 
                   pro = object$GMM$parameters$pro,
                   mean = transfGMM$mean,
                   sigma = transfGMM$sigma,
                   d = d,
                   nsamples = nsamples)
  Negentropy <- -EMC$Entropy + ppgmmga:::EntropyGauss(S = transfGMM$sz, d = d)
  
  FA$Negentropy <- Negentropy
  FA$se <- EMC$se
  FA$Z <- object$data %*% B
  return(FA)
}