\name{ppgmmgaAD-class}
\alias{ppgmmgaAD-class}
\docType{class}
\title{
Class \code{"ppgmmgaAD"}}
\description{
An S3 class object for ppgmmga algorithm for anomaly detection trought Gaussian mixture with uniform component.}
\usage{
ppgmmgaAD-class(x)
}
\section{Objects from the class}{
Objects can be created by calls to the \code{\link{ppgmmgaAD}} function.
}

\section{Values}{

\describe{
\item{data}{The original data.}
\item{d}{The subspace dimension.}
\item{approx}{The type of approximation used.}
\item{GMM}{Object of class \code{"mclust"} with the infomration from the estimated GMM density.}
\item{GA}{Object of class \code{"ga"} with the information from the genetic algorithm search. }
\item{Negentropy}{The value of the maximised Negentropy.}
\item{basis}{The basis of the subspace.}
\item{Z}{The projected data.}
\item{Class}{The final classification of the observation.}
}

}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}

\keyword{classes}

\seealso{See \code{\link{ppgmmgAD}}.}

