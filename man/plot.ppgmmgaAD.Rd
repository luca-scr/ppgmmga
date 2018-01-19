\name{plot.ppgmmgaAD}
\alias{plot.ppgmmgaAD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot ppgmmgaAD results
}
\description{
Plot the data onto a selected subsapce and marking the observation according to the classification returned by appling the ppgmmgaAD algorithm.
}
\usage{
plot.ppgmmgaAD(object, 
               dim = 1:object$d, 
               drawRegion = FALSE, 
               drawAxis = FALSE, 
               class = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
output from \code{ppgmmgaAD}
}
  \item{dim}{
%%     ~~Describe \code{dim} here~~
}
  \item{drawRegion}{
a logical defaulting to \code{FALSE}. If set a \code{TRUE} the region of the data is drawn in the plot.}
  \item{drawAxis}{
a logical defaulting to \code{FALSE}. If set a \code{TRUE} the projected axis are include in the plot.}
  \item{class}{
A numeric or character vector indicating the classification of the observations.}
}

\details{
Plot the cloud of points in the chosen subspace. The observations are marked according to the classification resulting from fitting the Gaussian mixture with  uniform component noise.}
\value{
Returns a object of class \code{ggplot}.
}
\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}
\references{
%% ~put references to the literature/web site here ~
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
