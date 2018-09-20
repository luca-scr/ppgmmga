\name{NegentropyGMM}
\alias{NegentropyGMM}
\title{
Approximated Negentropy for Gaussian mixture models (GMMs)
}
\description{
Computes the Negentropy for GMMs using different approximations.
}
\usage{
NegentropyGMM(G,
              pro, 
              mean, 
              sigma, 
              cov, 
              d,
              method = c("UT", "VAR", "SOTE", "MC"), 
              nsamples = 1e+05)
}
\arguments{
  \item{G}{
an integer value specifying the numbers of mixture components.
}
  \item{pro}{
a vector whose \emph{g}th component is the mixing proportion for the \emph{g}th component of the mixture model.
}
  \item{mean}{
a \eqn{p x G} matrix indicating the mean for each component of the mixture model. The \emph{g}th column represents the mean for the \emph{g}th
component.

If \code{d = 1}, also a vector of lenght \code{G} can be provided.}
  \item{sigma}{
an \eqn{p × p × G} array indicating the covariance matrix for the Gaussian mixture. The \emph{g}th matrix of this multidimensional array indicates the covariance matrix for the gth component.

If \code{d = 1}, also a vector of lenght \code{G} can be provided.}

  \item{cov}{
an \eqn{p × p} covariance matrix.}
  \item{d}{an integer value specifying the dimension of the data}
  \item{method}{
the type of approximation to compute the Negentropy for GMMs.
Possible values are:
  \tabular{ll}{
\code{"UT"} \tab for Unscented approximation. \cr
\code{"VAR"} \tab for Variational approximation. \cr
\code{"SOTE"} \tab for Second Order Taylor approximation. \cr
\code{"MC"} \tab Monte Carlo approximation. \cr
}
The default is "UT".}
  \item{nsamples}{
The sample size for the \code{"MC"} method. If a different method from \code{"MC"} is specified, this argument is ignored.}
}

\details{
The Negentropy index (Huber, 1985) describes the departure of a distribution from the Gaussian distribution. The Negentropy is defined as follows:
\deqn{J(z) = h(\phi(\mu_z, \Sigma_z)) - h(z),}
where \eqn{h(\phi(\mu_z, \Sigma_z))} is the entropy of the multivariate Normal distribution, \eqn{\phi(\mu_z, \Sigma_z)} is the multivariate Gaussian density of data \emph{z} with mean \eqn{\mu_z} and covariance matrix \eqn{\Sigma_z}, and \eqn{h(z)} is the entropy of the estimated density for the projected data. If GMMs are used to estimated the density of the data, the Negentropy does not have a closed formula, and then an approximation is needed. The first choice is to approximate the Negentropy using Monte Carlo method (Hershey and Olsen, 2007). Neverteless, such approximation is inefficient from a computational point of view. The Unscendent approximation (Goldberger and Aronowitz, 2005), the Variationa approximation (Hershey and Olsen, 2007), and the Second order Taylor approximation (Huber et al, 2008) provide a closed form approximation for the Negentropy. 
}
\value{
The function return a list with with the value of the approximated entropy and the method used. In addition, for the \code{"MC"} approximation the standard error is also reported. 
}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}

}

\references{
Goldberger, J. and Aronowitz, H. (2005). \emph{A distance measure between gmms based on
the unscented transform and its application to speaker recognition}. In INTERSPEECH,
pages 1985–1988. Citeseer.

Hershey, J. R. and Olsen, P. A. (2007). \emph{Approximating the Kullback Leibler divergence
between Gaussian mixture models}. In 2007 IEEE International Conference on Acoustics,
Speech and Signal Processing-ICASSP’07, volume 4, pages IV–317. IEEE.

Huber, P. J. (1985). \emph{Projection pursuit}. The Annals of Statistics, pages 435–475.

Huber, M. F., Bailey, T., Durrant-Whyte, H., and Hanebeck, U. D. (2008). \emph{On entropy
approximation for Gaussian mixture random vectors}. In Multisensor Fusion and Integra-
tion for Intelligent Systems, 2008. MFI 2008. IEEE International Conference on, pages
181–188. IEEE.}


\seealso{
\code{\link{ppgmmga}}, \code{\link{EntropyGMM}}
}

\examples{
\dontrun{
require(mclust)
data(iris)
X = iris[,-5]
mod <- densityMclust(data = X)

# Unscented Transformation
NentropyUT <- NegentropyGMM(G = mod$G,
                            pro = mod$parameters$pro, 
                            mean = mod$parameters$mean,
                            sigma = mod$parameters$variance$sigma, 
                            cov = cov(X), 
                            d = ncol(X),
                            method = "UT" )
# Variational approximation
NentropyVAR <- NegentropyGMM(G = mod$G,
                             pro = mod$parameters$pro, 
                             mean = mod$parameters$mean,
                             sigma = mod$parameters$variance$sigma, 
                             cov = cov(X), 
                             d = ncol(X),
                             method = "VAR" )

# Second order approximation
NentropySOTE <- NegentropyGMM(G = mod$G,
                              pro = mod$parameters$pro, 
                              mean = mod$parameters$mean,
                              sigma = mod$parameters$variance$sigma, 
                              cov = cov(X), 
                              d = ncol(X),
                              method = "SOTE" )

# Monte Carlo approximation
NentropyMC <- NegentropyGMM(G = mod$G,
                            pro = mod$parameters$pro, 
                            mean = mod$parameters$mean,
                            sigma = mod$parameters$variance$sigma, 
                            cov = cov(X),
                            d = ncol(X),
                            method = "MC" )
}
}

\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
