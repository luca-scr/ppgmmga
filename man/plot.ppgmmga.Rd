\name{plot.ppgmmga}
\alias{plot.ppgmmga}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot ppgmmga results}
\description{
Plot the data onto a selected subsapce by appling the ppgmmga algorithm.}
\usage{
plot.ppgmmga(object, 
             class = NULL, 
             dim = 1:object$d, 
             drawAxis = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
output from \code{ppgmmga}.
}
  \item{class}{
A numeric or character vector indicating the classification of the observations.}
  \item{dim}{
%%     ~~Describe \code{dim} here~~
}
  \item{drawAxis}{
a logical defaulting to \code{TRUE}. If set a \code{TRUE} the projected axis are include in the plot.}
  \item{\dots}{
further arguments.}
}
\details{
Plot the cloud of points in the chosen subspace.}
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
data(iris)
X = iris[,-5]

pp = ppgmmga(data = X, d = 2)
plot(pp)

# Plot with classes

Class = iris[,5]
plot(pp, class = Class)

# 3d projection

pp = ppgmmga(data = X, d = 2)
plot(pp)

# Visualise only the second and third dimension

plot(pp, dim = c(2,3))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
