\name{plot.ppgmmga}
\docType{methods}
\alias{plot.ppgmmga}
\title{
Plot method for \code{ppgmmga}}
\description{
Plots the data onto the subspace by appling the ppgmmga algorithm.}
\usage{
\method{plot}{ppgmmga}(x, 
                       class = NULL, 
                       dim = 1:x$d, 
                       drawAxis = TRUE, 
                       bins = nclass.Sturges,...)
}
\arguments{
  \item{x}{output from \code{ppgmmga}.}
  \item{class}{a numeric or character vector indicating the classification of the observations.}
  \item{dim}{ numeric vector indicating the dimensions to plot when the subspace dimension is greater that 3 (\code{p > 3}). The defaulting value is set to \code{1:x$d}, i.e. a scatterplot matrix of projected data. }
  \item{drawAxis}{a logical defaulting to \code{TRUE}. If set a \code{TRUE} the projected axis are include in the plot.}
  \item{bins}{a R function to compute the number of classes for the histogram. \code{nclass.Sturges} is used by default. For more detail on \code{nclass.Sturges} and other methods look at \code{\link{nclass.Sturges}} help. User may provide a different function. 
For multivariate projections such argument is ignored.}
  \item{\dots}{further arguments.}
}
\details{
Plots the cloud of points onto a subspace after appling the Projection Pursuit algorithm based on Gaussian mixtures and Genetic algorithm implemented in \code{ppgmmga} function.}
\value{
Returns a object of class \code{\link{ggplot}}.
}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}

}

\seealso{
\code{\link{ppgmmga}}}


\examples{
\dontrun{
require(mclust)
data(iris)
X = iris[,-5]
Class = iris$Species

pp = ppgmmga(data = X, d = 2)
plot(pp)

# Plot with classes

Class = iris[,5]
plot(pp, class = Class)

# 3d projection

pp = ppgmmga(data = X, d = 3)
plot(pp)

# Second and third dimension

plot(pp, dim = c(2,3))
}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
