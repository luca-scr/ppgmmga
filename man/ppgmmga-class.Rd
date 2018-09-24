\name{ppgmmga-class}
\alias{ppgmmga-class}
\docType{class}

\title{Class \code{'ppgmmga'}}

\description{An S3 class object for ppgmmga algorithm}

\section{Objects from the class}{
Object can be created by calls to the \code{\link{ppgmmga}} function.
}

\section{Values}{
\describe{
\item{data}{The input data matrix.}
\item{d}{The dimension of the projection subspace.}
\item{approx}{The type of approximation used for computing negentropy.}
\item{GMM}{An object of class \code{'densityMclust'} containing the Gaussian mixture density estimation. See \code{\link[mclust]{densityMclust}} for details.}
\item{GA}{An object of class \code{'ga'} containing the Genetic Algorithm search. See \code{\link[GA]{ga}} for details.}
\item{Negentropy}{The value of maximised negentropy.}
\item{basis}{The matrix basis of the projection subspace.}
\item{Z}{The matrix of projected data.}
}
}

\author{
Serafini A. \email{srf.alessio@gmail.com}\cr
Scrucca L. \email{luca.scrucca@unipg.it}
}

\keyword{classes}

\seealso{\code{\link{ppgmmga}}, \code{\link{plot.ppgmmga}}, \code{\link{summary.ppgmmga}}}

