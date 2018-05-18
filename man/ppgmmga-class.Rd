\name{ppgmmga-class}
\alias{ppgmmga-class}
\docType{class}
\title{
Class \code{"ppgmmga"}}
\description{
An S3 class object for ppgmmga algorithm.}

\section{Objects from the class}{
Objects can be created by calls to the \code{\link{ppgmmga}} function.
}

\section{Values}{

\describe{
\item{data}{The original data.}
\item{d}{The subspace dimension.}
\item{approx}{The type of approximation used.}
\item{GMM}{Object of class \code{"mclust"} with the infomration for  GMM density estimations.}
\item{GA}{Object of class \code{"ga"} with the information from the genetic algorithm search. }
\item{Negentropy}{The value of maximised Negentropy.}
\item{basis}{The basis of subspace.}
\item{Z}{The projected data.}
}

}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}

\keyword{classes}

\seealso{\code{\link{ppgmmga}}}

