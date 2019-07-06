\name{ppgmmga.options}
\alias{ppgmmga.options}

\title{
Default values for \pkg{ppgmmga} package}

\description{
Set or retrieve default values to be used by the \pkg{ppgmmga} package.
}

\usage{ppgmmga.options(\dots)}

\arguments{
\item{\dots}{A single character vector, or a named list with components. In the one argument case, the form \code{name = value} can be used to change a single option. In the multiple arguments case, the form \code{list(name1 = value1, name2 = value2)} can be used to change several arguments. If no arguments are provided, then the function returns all the current options. 
For the available options see the Details section below.}
}

\details{
This function can be used to set or retrieve the values to be used by the \pkg{ppgmmga} package.

The function globally sets the arguments for the current session of R. The default options are restored with a new R session. To temporarily change the options for a single call to \code{ppgmmga} function, look at \code{options} argument in \code{\link{ppgmmga}}. 

Available options are:

\describe{
\item{\code{modelNames}}{A string specifying the GMM to fit. See \code{\link[mclust]{mclustModelNames}} for the available models.}

\item{\code{G}}{An integer value or a vector of integer values specifying the number of mixture components. If more than a single value is provided, the best model is selected using the BIC criterion. By default \code{G = 1:9}.}

\item{\code{initMclust}}{A string specifying the type of initialisation to be used for the EM algorithm. See \code{\link[mclust]{mclust.options}} for more details.}
  
\item{\code{popSize}}{The GA population size. By default \code{popSize = 100}.}
  
\item{\code{pcrossover}}{The probability of crossover. By default \code{pcrossover = 0.8}.}
  
\item{\code{pmutation}}{The probability of mutation. By default \code{pmutation = 0.1}.}
  
\item{\code{maxiter}}{An integer value specifying the maximum number of iterations before stopping the GA. By default \code{maxiter = 1000}.}

\item{\code{run}}{An integer value indicating the number of generations without improvment in the best value of fitness fuction. \code{run = 100}.}

\item{\code{selection}}{An \code{R} function performing the selection genetic operator. See \code{\link{ga_Selection}} for details. By default \code{selection = gareal_lsSelection}.}

\item{\code{crossover}}{An \code{R} function performing the crossover genetic operator. See \code{\link{ga_Crossover}} for details. By default \code{crossover = gareal_laCrossover}.}

\item{\code{mutation}}{An \code{R} function performing the mutation genetic operator. See \code{\link{ga_Mutation}} for details. By default \code{mutation = gareal_raMutation}.}

\item{\code{parallel}}{A logical value specifying whether or not GA should be run in parallel. By default \code{parallel = FALSE}.}

\item{\code{numIslands}}{An integer value specifying the number of islands to be used in the Island Genetic Algorithm. By default \code{numIslands = 4}.}

\item{\code{migrationRate}}{A value specifying the fraction of migration between islands. By default \code{migrationRate = 0.1}.}

\item{\code{migrationInterval}}{An integer values specifying the number of generations to run before each migration. By default \code{migrationInterval = 10}.}

\item{\code{optim}}{A logical value specifying whether or not a local search should be performed. By default \code{optim = TRUE}.}

\item{\code{optimPoptim}}{A value specifying the probability a local search is performed at each GA generation. By default \code{optimPoptim = 0.05}.}

\item{\code{optimPressel}}{A value in \eqn{[0,1]} specifying the pressure selection. Values close to 1 tend to assign higher selection probabilities to solutions with higher fitness, whereas values close to 0 tend to assign equal selection probability to any solution. By default \code{optimPressel = 0.5}.}

\item{\code{optimMethod}}{A string specifying the general-purpose optimisation method to be used for local search. See \code{\link{optim}} for the available algorithms. By default \code{optimMethod = "L-BFGS-B"}.}

\item{\code{optimMaxit}}{An integer value specifying the number of iterations for the local search algorithm. By default \code{optimMaxit = 100}.}

\item{\code{orth}}{A string specifying the method employed to orthogonalise the matrix basis. Available methods are the QR decomposition \code{"QR"}, and the Singular Value Decomposition \code{"SVD"}. By default \code{orth = "QR"}.}
}

For more details about options related to Gaussian mixture modelling see \code{\link[mclust]{densityMclust}}, and for those related to genetic algorithms see \code{\link[GA]{ga}} and \code{\link[GA]{gaisl}}. 
}

\author{
Serafini A. \email{srf.alessio@gmail.com}\cr
Scrucca L. \email{luca.scrucca@unipg.it}
}

\references{
Scrucca, L., Fop, M., Murphy, T. B., & Raftery, A. E. (2016) mclust 5: Clustering, classification and density estimation using gaussian finite mixture models. \emph{The R journal}, 8(1), 205-233. https://journal.r-project.org/archive/2016-1/scrucca-fop-murphy-etal.pdf

Scrucca, L. (2013) GA: A Package for Genetic Algorithms in R. \emph{Journal of Statistical Software}, 53(4), 1-37. http://www.jstatsoft.org/v53/i04/

Scrucca, L. (2017) On some extensions to GA package: hybrid optimisation, parallelisation and islands evolution. \emph{The R Journal}, 9/1, 187-206. https://journal.r-project.org/archive/2017/RJ-2017-008

Scrucca, L. and Serafini, A. (2019) Projection pursuit based on Gaussian mixtures and evolutionary algorithms. *Journal of Computational and Graphical Statistics*, DOI: 10.1080/10618600.2019.1598871 
}

\seealso{\code{\link{ppgmmga}}}

\examples{
\dontrun{
ppgmmga.options()

# Print a single option
ppgmmga.options("popSize")

# Change (globally) an option
ppgmmga.options("popSize" = 10)
ppgmmga.options("popSize")
}
}
\keyword{options}
