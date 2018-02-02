\name{ppgmmgaAD}
\alias{ppgmmgaAD}

\title{
Projection pursuit based on Gaussian mixtures and evolutionary algorithms for anomaly detection 	through Gaussian mixture with uniform component
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ppgmmgaAD(data, 
          G = NULL, 
          modelNames = NULL, 
          volumeMethod = c("box", "pc", "ConvHull"), 
          d = 2, 
          approx = c("UT", "VAR", "SOTE"), 
          center = TRUE, 
          scale = FALSE, 
          gmm = NULL, 
          gatype = c("ga", "gaisl"),
          cutoff = 0.975, 
          noiseInit = NULL, 
          alpha = 0.5, 
          projData = NULL,
          opt = list(), 
          monitor = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a \emph{n x p} matrix containing the data with rows correspond to observations and columns correspond to variables.}
  \item{G}{
an integer vector or value specifying the numbers of mixture components. The default is \emph{G=1:9} and the BIC is computed to select the model.
}
  \item{modelNames}{
A string vector or value indicating the covariance model used to fit the mixture with the EM algorithm. See \code{\link[mclust]{mclustModelNames}} for the list of available models.}
  \item{volumeMethod}{
a string indicating the method to compute the volume of the data.
Possible values are:
   \tabular{ll}{
\code{"box"} \tab for the smallest hyperrectangle with sides parallel to the coordinate axes containing the data. \cr
\code{"pc"} \tab for the smallest hyperrectangle with sides parallel to the principal component axes. \cr
\code{"ConvHull"} \tab for the Convex Hull of the data.}
}
  \item{d}{
the dimension of the subspace onto which the data are projected.}
  \item{approx}{
the type of approximation to compute the Negentropy for GMM
model.
Possible values are:
  \tabular{ll}{
\code{"UT"} \tab for Unscented approximation. \cr
\code{"VAR"} \tab for Variational approximation. \cr
\code{"SOTE"} \tab for Second Order Taylor approximation. \cr
\code{"MC"} \tab Monte Carlo approximation. \cr}
}
  \item{center}{
a logical value set to \emph{TRUE} indicating whther or not the data are centred.
}
  \item{scale}{
a logical value set to \emph{FALSE} indicating whther or not the data are scaled.
}
  \item{gmm}{
%%     ~~Describe \code{gmm} here~~
}
  \item{gatype}{
the type of the genetic algoritm used to maximised the Negentropy.
Possible values are:
  \tabular{ll}{
\code{"ga"} \tab for genetic algorithm. \cr
\code{"gaisl"} \tab for Island genetic algorithm. \cr}
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{noiseInit}{
%%     ~~Describe \code{noiseInit} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{projData}{
%%     ~~Describe \code{projData} here~~
}

  \item{opt}{
%%     ~~Describe \code{opt} here~~
}
  \item{monitor}{
a logical value set to \emph{TRUE} to show the evolution of the evolution of the optimisation with the genetic algorith.
}
  \item{\dots}{
additional arguments.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an object of class \code{ppgmmgaAD-class}. See \code{\link{ppgmmgaAD-class}} for a description of the object.
}
\references{
Banfield, J. D. and Raftery, A. E. (1993). \emph{Model-based Gaussian and non-Gaussian clustering.} Biometrics, pages
803–821.
}
\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
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

ppAD = ppgmmgaAD(data = X)

# Fixing the number of components 
ppAD = ppgmmgaAD(data = X, G = 2)

# Fixing the model 
ppAD = ppgmmgaAD(data = X, modelNames = "EEE")

# Example changing the size of the initial population in the general options

ppAd = ppgmmgaAD(data = X, opt = list("popSize" = 10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
