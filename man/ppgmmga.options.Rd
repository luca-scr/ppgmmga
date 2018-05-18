\name{ppgmmga.options}
\alias{ppgmmga.options}

\title{
Default values for \code{ppgmmga} package}

\description{
Set or retrieve default values for use with \code{ppgmmga} package.
}

\usage{ppgmmga.options(\ldots)}
\arguments{
\item{\dots}{one argument in the form \code{name = value} to change a single option or more argument in the form \code{list(name1 = value1, name2 = value2)}. If no arguments are provided, the function print all the options. For more details on the argument look at Detail section below.}
}

\details{
The function assigns values to the \code{ga} function for genetic algorithms, \code{Mclust} function to estimate Guassim mixture, and to other values in the \code{ppgmmga} function. 

Available options are:

\describe{
  \item{\code{orth}}{a}
  \item{popSize}{a}
  \item{pcr}{a}
  \item{pm}{a}
  \item{\code{maxiter}}{an integer value indicating the maximum number of iteration before stopping the Genetic algorithm. The value is set by default to 100}
\item{run}{an integer value indicating the number of generation without improvment in the best value of fitness fuction befor stopping the algorithm. The value is set by default to 100}
\item{optim}{}
\item{selection}{}
\item{crossover}{}
\item{mutation}{}
\item{parallel}{a logical value, set to \code{FALSE} bu=y default, for parallel Genetic algorithms.}
\item{numIslands}{}
\item{migrationRate}{}
\item{migrationInterval}{}
\item{optimMethod}{}
\item{poptim}{}
\item{pressel}{}
\item{fnscale}{}
\item{maxit}{a logical value indicating whether or not lcoal search is performed. the value is set to \code{TRUE}}
\item{initMclust}{a string indicating the type of initialisation for the mixture model. Look \code{mclust.options} for the available initialisation.}
\item{modelNames}{a string value or vector containing the type of covariance parametrisation for the mixtures. Look at \code{mclustModelNames} for the available models.}
\item{G}{an integer value or vactor containing the number of mixtures component to estimate}

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
