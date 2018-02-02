\name{ppgmmga.options}
\alias{ppgmmga.options}

\title{
Default values for \code{ppgmmga} package}

\description{
Set or retrieve default values for use with \code{ppgmmga} package.
}

\usage{ppgmmga.options(\ldots)}

\details{
The function assigns values to the \code{ga} function for genetic algorithms, \code{Mclust} function to estimate the Guassima mixture, and to other value in the \code{ppgmmga} function. 

\describe{

\item{}{}


}
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}

\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
ppgmmga.options()

# Single option

ppgmmga.options("popSize")

# Change option

ppgmmga.options("popSize" = 10)
}
\keyword{datasets}
