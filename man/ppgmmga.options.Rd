\name{ppgmmga.options}
\alias{ppgmmga.options}

\title{
Default values for \code{ppgmmga} package}

\description{
Set or retrieve default values for use with \code{ppgmmga} package.
}

\usage{ppgmmga.options(\ldots)}

\details{
The function assigns values to the \code{ga} function for genetic algorithms, \code{Mclust} function to estimate Guassim mixture, and to other values in the \code{ppgmmga} function. 

\describe{

\item{}{}


}
}


\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}

\seealso{
\code{\link{ppgmmga.options.restore}}
}
\examples{
ppgmmga.options()

# Print option

ppgmmga.options("popSize")

# Change option

ppgmmga.options("popSize" = 10)
}
\keyword{datasets}
