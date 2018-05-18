\name{plot.ppgmmga}
\alias{plot.ppgmmga}
\title{
Plot \code{ppgmmga} results}
\description{
Plot the data onto a selected subsapce by appling the ppgmmga algorithm.}
\usage{
\method{plot}{ppgmmga}(object, 
             class = NULL, 
             dim = 1:object$d, 
             drawAxis = TRUE, 
             bins = nclass.Sturges,...)
}
\arguments{
  \item{object}{
output from \code{ppgmmga}.
}
  \item{class}{
A numeric or character vector indicating the classification of the observations.}
  \item{dim}{A numeric vector indicating the dimension to plot when the subspace dimension \code{p > 3}. The defaulting value is set to \code{1:object$d}, i.e. a scatter plot for each pair of variables of the projected data. 
}
  \item{drawAxis}{
a logical defaulting to \code{TRUE}. If set a \code{TRUE} the projected axis are include in the plot.}
  \item{bins}{function to compute the number of classes for a histogram. \code{nclass.Sturges} is used by default. For more detail on \code{nclass.Sturges} and other methods look at \code{nclass} help. User may specified any type of one argument function to compute. For multivariate projections such argument is ignored.}
  \item{\dots}{
further arguments.}
}
\details{
Plot the cloud of points into chosen subspace after appling the projection pursuit algorithm based on Gaussian mixtures and Genetic algorithm, implemented in \code{ppgmmga} function.}
\value{
Returns a object of class \code{ggplot}.
}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
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

pp = ppgmmga(data = X, d = 3)
plot(pp)

# Visualise only the second and third dimension

plot(pp, dim = c(2,3))


}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
