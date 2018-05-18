\name{ppgmmga.options}
\alias{ppgmmga.options}

\title{
Default values for \code{ppgmmga} package}

\description{
Set or retrieve default values for use with \code{ppgmmga} package.
}

\usage{ppgmmga.options(\ldots)}
\arguments{
\items{\dots}{one argument in the form \code{name = value} to change a single option or more argument in the form \code{list(name1 = value1, name2 = value2)}. If no arguments are provided, the function print all the options. For more details on the argument look at Detail section below.}
}
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
