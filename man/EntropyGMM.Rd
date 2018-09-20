\name{EntropyGMM}
\alias{EntropyGMM}

\title{
Approximated entropy for Gaussian mixture models (GMMs)
}

\description{
Computes the entropy for GMMs using different approximations.
}

\usage{
EntropyGMM(G,
           pro,
           mean,
           sigma,
           d,
           method = c("UT", "VAR", "SOTE", "MC"),
           nsamples = 1e+05)
}
           
\arguments{
  \item{G}{
an integer value specifying the numbers of mixture components.}
  \item{pro}{
a vector with \emph{g}th component is the mixing proportion for the \emph{g}th component of the mixture model.
}
  \item{mean}{
a \eqn{p x G} matrix indicating the mean for each component of the
mixture model. The \emph{g}th column represents the mean for the \emph{g}th
component. 

If \code{d = 1}, also a vector of lenght \code{G} can be provided.}
  \item{sigma}{an \eqn{p × p × G} array indicating the covariance matrix for the Gaussian mixture. The \emph{g}th matrix of this multidimensional array indicates the covariance matrix for the gth component.

If \code{d = 1}, also a vector of lenght \code{G} can be provided.}
  \item{d}{an integer value specifying the dimension of the data.}

  \item{method}{
the type of approximation to compute the entropy for GMMs.
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
Gaussian mixture models (GMMs) do not have a closed form solution for the entropy, and then an approximation needed. The first choice is to approximate the entropy using Monte Carlo methods (Hershey and Olsen, 2007). Neverteless, such approximation is inefficient from a computational point of view. The Unscendent approximation (Goldberger and Aronowitz, 2005), the Variationa approximation (Hershey and Olsen, 2007), and the Second order Taylor approximation (Huber et al, 2008) provide a closed form approximation for the GMMs entropy.}

\value{
The function returns a list with the value of the approximated entropy and the method used. In addition, for the \code{"MC"} approximation the standard error is also reported. 
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

\seealso{
\code{\link{ppgmmga}}, \code{\link{NegentropyGMM}}
}

\examples{
\dontrun{
require(mclust)
data(iris)
X = iris[,-5]
mod = densityMclust(data = X)

# Unscented Transformation
entropyUT = EntropyGMM(G = mod$G,
                       pro = mod$parameters$pro,
                       mean = mod$parameters$mean,
                       sigma = mod$parameters$variance$sigma,
                       d = ncol(X),
                       method = "UT")

# Variational approximation
entropyVAR = EntropyGMM(G = mod$G,
                        pro = mod$parameters$pro,
                        mean = mod$parameters$mean,
                        sigma = mod$parameters$variance$sigma,
                        d = ncol(X),
                        method = "VAR")

# Second order approximation
entropySOTE = EntropyGMM(G = mod$G,
                         pro = mod$parameters$pro,
                         mean = mod$parameters$mean,
                         sigma = mod$parameters$variance$sigma,
                         d = ncol(X),
                         method = "SOTE")

# Monte Carlo approximation
entropyMC = EntropyGMM(G = mod$G,
                       pro = mod$parameters$pro,
                       mean = mod$parameters$mean,
                       sigma = mod$parameters$variance$sigma,
                       d = ncol(X),
                       method = "MC")
}
}

\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
