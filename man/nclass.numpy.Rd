\name{nclass.numpy}
\alias{nclass.numpy}

\title{
Compute the Number of Classes for a Histogram}

\description{
Compute the number of classes for a histogram as the maximum of the "Sturges" and "FD" (Freedman Diaconis) estimators as in numpy library for Python. 
}

\usage{nclass.numpy(x, \dots)}

\arguments{
\item{x}{A vector of values.}
\item{\dots}{Further arguments passed to or from other methods.}
}

\author{
Scrucca L. \email{luca.scrucca@unipg.it}
}

\seealso{\code{\link{nclass.Sturges}}, \code{\link{nclass.FD}}}

\examples{
\dontrun{
library(ggplot2)
x <- rnorm(100)
ggplot() + geom_histogram(aes(x), col = "grey92", bins = nclass.numpy(x))
x <- rnorm(1000)
ggplot() + geom_histogram(aes(x), col = "grey92", bins = nclass.numpy(x))
n = c(50, seq(100,1000,by=100))
brks = rep(NA, length(n))
for(i in seq(n)) brks[i] = nclass.numpy(rnorm(n[i]))
ggplot() + geom_point(aes(x = n, y = brks))
}
}

