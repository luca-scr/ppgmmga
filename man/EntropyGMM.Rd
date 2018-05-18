\name{EntropyGMM}
\alias{EntropyGMM}
\title{
Approximated entropy for Gaussian mixture models
}
\description{
Entropy approximation for Gaussian mixture models (GMMs).

}
\usage{
EntropyGMM(G,
           d,
           pro,
           mean,
           sigma,
           method = c("UT", "VAR", "SOTE", "MC"),
           nsamples = 1e+05)
}
\arguments{

  \item{G}{
an integer value specifying the numbers of mixture components.}
  \item{d}{an integer value specifying the dimension of the data}
  \item{pro}{
a vector whose \emph{g}th component is the mixing proportion for the \emph{g}th component of the mixture model.
}
  \item{mean}{
a \eqn{p x G} matrix indicating the mean for each component of the
mixture model. The \emph{g}th column represent the mean for the gth
component. 

If \code{d = 1}, also a vector of lenght \code{G} can be provided.}
  \item{sigma}{
an \eqn{p × p × G} array indicating the covariance matrix for the Gaus-
sian mixture. The \emph{g}th matrix of this multidimensional array indicates the covariance matrix for the gth component.

If \code{d = 1}, also a vector of lenght \code{G} can be provided.}

  \item{method}{
the type of approximation to compute the Entropy for GMM
model.
Possible values are:
  \tabular{ll}{
\code{"UT"} \tab for Unscented approximation. \cr
\code{"VAR"} \tab for Variational approximation. \cr
\code{"SOTE"} \tab for Second Order Taylor approximation. \cr
\code{"MC"} \tab Monte Carlo approximation. \cr
}
The default is "UT".
}
  \item{nsamples}{
The sample size for the \code{"MC"} method. If a different method FROM \code{"MC"} is specified, this argument is ignored.}
}
\details{
Gaussian mixture models (GMMs) do not have a closed form for the entropy, and then an approximation needed. The first choice is to approximate the entropy using Monte Carlo methods (Hershey and Olsen, 2007). Neverteless such apprximation is the only one to converge to the real Entropy, the MC method is inefficient from a computational point of view. The Unscendent approximation (Goldberger and Aronowitz, 2005), the Variationa approximation (Hershey and Olsen, 2007), and the Second order Taylor approximation (Huber et al, 2008) approximate the entropy in a more efficient way.

If the number of components \code{G} is equal to 1, then the Gaussian entropy is computed.}
\value{
For \code{"UT"}, \code{"VAR"}, and \code{"SOTE"} methods the function returns the value of the approximated entropy. For the \code{"MC"} method the function returns a list with the value of entropy and standard error.
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

Huber, M. F., Bailey, T., Durrant-Whyte, H., and Hanebeck, U. D. (2008). \emph{On entropy
approximation for Gaussian mixture random vectors}. In Multisensor Fusion and Integra-
tion for Intelligent Systems, 2008. MFI 2008. IEEE International Conference on, pages
181–188. IEEE.

}

\examples{

require(mclust)
data(iris)
X = iris[,-5]
mod = Mclust(data = X)

# Unscented Transformation
entropyUT = EntropyGMM(G = mod$G,
                       d = ncol(X),
                       pro = mod$parameters$pro,
                       mean = mod$parameters$mean,
                       sigma = mod$parameters$variance$sigma,
                       method = "UT")

# Variational approximation
entropyVAR = EntropyGMM(G = mod$G,
                        d = ncol(X),
                        pro = mod$parameters$pro,
                        mean = mod$parameters$mean,
                        sigma = mod$parameters$variance$sigma,
                        method = "VAR")

# Second order approximation
entropySOTE = EntropyGMM(G = mod$G,
                         d = ncol(X),
                         pro = mod$parameters$pro,
                         mean = mod$parameters$mean,
                         sigma = mod$parameters$variance$sigma,
                         method = "SOTE")

# Monte Carlo approximation
entropyMC = EntropyGMM(G = mod$G,
                       d = ncol(X),
                       pro = mod$parameters$pro,
                       mean = mod$parameters$mean,
                       sigma = mod$parameters$variance$sigma,
                       method = "MC")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
