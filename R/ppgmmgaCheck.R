ppgmmgaCheck <- function(...,
                         AsList = NULL,
                         nSamples = 1e5,
                         level = 0.05,
                         names = NULL)
{
  
  if(is.null(AsList))
  {
    mods <- list(...)
  } else
  {
    mods <- AsList
  }
  n <- length(mods)
  
  approx <- rep(NA, n)
  Neg <- rep(NA, n)
  NegMC <- rep(NA,n)
  int <- matrix(NA, nrow = n, ncol = 2)
  se <-  rep(NA, n)
  for (i in 1:n)
  {  
    mod <- mods[[i]]
    approx[i] <- mods[[i]]$approx
    MC <- NegentropyMC(par = mod$GA@solution,
                       GMM =mod$GMM,
                       p = mod$GMM$d,
                       d = mod$d,
                       nSample = nSamples,
                       level = level)
    Neg[i] <- mod$Negentropy
    NegMC[i] <- MC$Negentropy
    int[i,] <- MC$confint
    se[i] <- MC$se    
  }
  
  if(is.null(names))
  {  
    names <- names(mods)
    if(is.null(names))  names <- approx
  }
  
  tab <- t(data.frame("Approx Negentropy" = Neg,
                     "MC Negentropy" = NegMC,
                     "Relative accuracy" = Neg/NegMC,
                     "MC se" = se,
                     # confint = int,
                     row.names = make.names(names,unique = TRUE),
                     check.names = FALSE))
  return(tab)
}