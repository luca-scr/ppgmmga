\name{EntropyGMM}
\alias{EntropyGMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximated Entropy for Gaussian mixture models
}
\description{
Approximation of the Entropy for Gaussian mixture models (GMMs).

}
\usage{
EntropyGMM(G,
           pro,
           mean,
           sigma,
           method = c("UT", "VAR", "SOTE", "MC"),
           nsamples = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{G}{
an integer value specifying the numbers of mixture components.}
  \item{pro}{
a vector whose \emph{g}th component is the mixing proportion for the \emph{g}th component of the mixture model.
}
  \item{mean}{
a \eqn{p x G} matrix indicating the mean for each component of the
mixture model. The \emph{g}th column represent the mean for the gth
component.}
  \item{sigma}{
an \eqn{p × p × G} array indicating the covariance matrix for the Gaus-
sian mixture. The \emph{g}th matrix of this multidimensional array indi-
cates the covariance matrix for the gth component.}

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
The sample size for the \code{"MC"} method.}
}
\details{
Gaussian mixture models (GMMs) do not have a closed form for the Entropy, and then an approximation needed. The first choice is to approximate the Entropy using Monte Carlo method (Hershey and Olsen, 2007). Neverteless such apprximation is the only one to converge to the real Entropy, the MC method is inefficient from a computational point of view. The Unscendent approximation (Goldberger and Aronowitz, 2005), the Variationa approximation (Hershey and Olsen, 2007), and the Second order Taylor approximation (Huber et al, 2008) approximate the Entropy in more eficient way.}
\value{
For "UT", "VAR", and "SOTE" method the function returns the value of the approximated Entropy. For the "MC" method the function returns a list with the value of the Entropy and the standard error.
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

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(iris)
X = iris[,-5]
mod = Mclust(data = X)

# Unscented Transformation
entropyUT = EntropyGMM(G = mod$G,
                     pro = mod$parameters$pro,
                     mean = mod$parameters$mean,
                     sigma = mod$parameters$variance$sigma,
                     method = "UT")

# Variational approximation
entropyVAR = EntropyGMM(G = mod$G,
                      pro = mod$parameters$pro,
                      mean = mod$parameters$mean,
                      sigma = mod$parameters$variance$sigma,
                      method = "VAR")

# Second order approximation
entropySOTE = EntropyGMM(G = mod$G,
                      pro = mod$parameters$pro,
                      mean = mod$parameters$mean,
                      sigma = mod$parameters$variance$sigma,
                      method = "SOTE")

# Second order approximation
entropyMC = EntropyGMM(G = mod$G,
                      pro = mod$parameters$pro,
                      mean = mod$parameters$mean,
                      sigma = mod$parameters$variance$sigma,
                      method = "MC")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
