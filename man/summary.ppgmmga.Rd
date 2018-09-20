\name{summary.ppgmmga}
\docType{methods}
\alias{summary.ppgmmga}
\alias{print.summary.ppgmmga}

\title{
Summary for Projection Pursuit based on Gaussian Mixtures and Evolutionary Algorithms}
\description{
Summary method for class \code{ppgmmga}.}
\usage{
\method{summary}{ppgmmga}(object, check = FALSE, ...)
                
\method{print}{summary.ppgmmga}(x, digits = getOption("digits"), ...)
}
\arguments{
  \item{object}{an object of class \code{ppgmmga}.}
  \item{check}{a logical value set to \code{FALSE}. If \code{TRUE} a Monte Carlo Negentropy approximation check is performed.}
\item{x}{an object of class \code{summary.ppgmmga}.}
\item{digits}{the number of significant digits.}
\item{\dots}{further arguments.}
}

\value{
The summary function returns an object of class \code{summary.ppgmmga} which can be printed by the corresponding print method. A list with the information from the \code{ppgmmga} algorithm is returned.

The argument \code{check = TRUE} compares the value of approximated Negentropy and the Monte Carlo Negentropy considering the same basis generated from the algorithm. Since the Monte Carlo approximation is the only to converge to the real Negentropy, it allows to check if approximated values are closed to the real Negentropy. The ratio between the approximated value from the algorithm and the value computed with Monte Carlo approximation is called Relative Accuracy. Such value should be close to \code{1} for a good approximation. 
}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}

}

\seealso{
\code{\link{ppgmmga}}
}
\examples{
\dontrun{
require(mclust)
data(iris)
X = iris[,-5]
Class <- iris$Species

pp <- ppgmmga(data = X, d = 2)  
summary(pp, check = TRUE)
}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
