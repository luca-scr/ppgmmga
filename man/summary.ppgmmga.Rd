\name{summary.ppgmmga}
\docType{methods}
\alias{summary.ppgmmga}
\alias{print.summary.ppgmmga}

\title{Summary for projection pursuit based on Gaussian mixtures and evolutionary algorithms for data visualisation}

\description{Summary method for objects of class \code{'ppgmmga'}.}

\usage{
\method{summary}{ppgmmga}(object, check = (object$approx != "none"), \dots)
                
\method{print}{summary.ppgmmga}(x, digits = getOption("digits"), \dots)
}

\arguments{
  \item{object}{An object of class \code{'ppgmmga'} as returned by \code{\link{ppgmmga}}.}
  \item{check}{A logical value specifying whether or not a Monte Carlo negentropy approximation check should be performed. By default is \code{FALSE} for exact negentropy calculation and \code{TRUE} for approximated negentropy.}
\item{x}{An object of class \code{summary.ppgmmga}.}
\item{digits}{The number of significant digits.}
\item{\dots}{Further arguments passed to or from other methods.}
}

\value{
The summary function returns an object of class \code{summary.ppgmmga} which can be printed by the corresponding print method. A list with the information from the \code{ppgmmga} algorithm is returned.

If the optional argument \code{check = TRUE} then the value of negentropy is compared to the Monte Carlo negentropy calculated for the same optimal projection basis selected by the algorithm. 
By default, it allows to check if the value returned by the employed approximation is closed to the Monte Carlo approximation of to the "true" negentropy. 
The ratio between the approximated value returned by the algorithm and the value computed with Monte Carlo is called Relative Accuracy. Such value should be close to 1 for a good approximation. 
}

\author{
Serafini A. \email{srf.alessio@gmail.com}\cr
Scrucca L. \email{luca.scrucca@unipg.it}
}

\seealso{\code{\link{ppgmmga}}}
