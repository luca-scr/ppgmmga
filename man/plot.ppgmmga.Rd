\name{plot.ppgmmga}
\docType{methods}
\alias{plot.ppgmmga}

\title{Plots the data onto the projection subspace estimated by the ppgmmga algorithm}

\description{Plot method for objects of class \code{'ppgmmga'}.}

\usage{
\method{plot}{ppgmmga}(x, 
     class = NULL, 
     dim = seq(x$d), 
     drawAxis = TRUE, 
     bins = nclass.Sturges,
     \dots)
}

\arguments{
  \item{x}{An object of class \code{'ppgmmga'} obtained from a call to \code{\link{ppgmmga}} function.}
  \item{class}{A numeric or character vector indicating the classification of the observations/cases to be plotted.}
  \item{dim}{A numeric vector indicating the dimensions to use for plotting. 
  By default, all the dimensions of the projection subspace (i.e. \code{x$d}) are used. Subsets of all the available dimensions can also be provided (see example below.)\cr
  The resulting graph depends on the dimension: in the 1D case a histogram is provided, a scatterplot in the 2D case, a scatterplot matrix in higher dimensions.
  }
  \item{drawAxis}{A logical value specifying whether or not the axes should be included in the 2D scatterplot. By default is to \code{TRUE}.}
  \item{bins}{An R function to be used for computing the number of classes for the histogram. By default \code{nclass.Sturges} is used; see \code{\link{nclass.Sturges}} for more details. Users may provide a different function. This argument only applies to 1D graphs.}
  \item{\dots}{further arguments.}
}

\details{
Plots the cloud of points onto a subspace after appling the Projection Pursuit algorithm based on Gaussian mixtures and Genetic algorithm implemented in \code{ppgmmga} function.}

\value{
Returns a object of class \code{\link{ggplot}}.
}

\author{
Serafini A. \email{srf.alessio@gmail.com}\cr
Scrucca L. \email{luca.scrucca@unipg.it}
}

\references{
Scrucca, L. and Serafini, A. (2019) Projection pursuit based on Gaussian mixtures and evolutionary algorithms. *Journal of Computational and Graphical Statistics*, DOI: 10.1080/10618600.2019.1598871 
}

\seealso{\code{\link{ppgmmga}}}

\examples{
\dontrun{
data(iris)
X <- iris[,-5]
Class <- iris$Species

# 1D
pp1 <- ppgmmga(data = X, d = 1, approx = "UT")
summary(pp1, check = TRUE)
plot(pp1)
plot(pp1, Class)

# 2D
pp2 <- ppgmmga(data = X, d = 2, approx = "UT")
summary(pp2, check = TRUE)
plot(pp2)
plot(pp2, Class)

# 3D
pp3 <- ppgmmga(data = X, d = 3)
summary(pp3, check = TRUE)
plot(pp3)
plot(pp3, Class)
plot(pp3, Class, dim = c(1,3))
plot(pp3, Class, dim = c(2,3))
}
}
