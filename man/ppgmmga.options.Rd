\name{ppgmmga.options}
\alias{ppgmmga.options}

\title{
Default values for \code{ppgmmga} package}

\description{
Set default values for use with \code{ppgmmga} package.
}

\usage{ppgmmga.options(\ldots)}
\arguments{
\item{\dots}{one argument in the form \code{name = value} to change a single option or \code{list(name1 = value1, name2 = value2)} to change more arguments. If no arguments are provided, the function prints all the options. For the available options look at Details section below.}}

\details{
The function assigns values  \code{densityMclust()} function to estimate Gaussian mixture, and to \code{ga} function for Genetic Algorithms. For the arguments related to GMMs (\code{densityMclust()}) see \code{\link{Mclust}} and \code{\link{mclust.options}} for more details. For arguments related to Genetic Algorithms (\code{ga}) see \code{\link{ga}} and \code{\link{gaisl}} for more details. 

Available options are:

\describe{
\item{\code{modelNames}}{A string indicating the GMMs to fit. See \code{\link{mclustModelNames}} for the available models.}

\item{\code{G}}{A integer vector or an integer value indicating the number of GMMs component. If a integer vector is provided, the best model is selected by BIC criterion. By default \code{G = 1:9}.}

  \item{\code{initMclust}}{A sting indicating the type of initialisation for EM algorithm used to fit the GMMs. See \code{\link{mclust.options}} for more details.}
  
  \item{\code{popSize}}{The population size. By default \code{popSize = 100}}
  
  \item{pcrossover}{The probability of crossover. By default \code{pcrossover = 0.8}.}
  
  \item{pmutation}{The probability of mutation. By default \code{pmutation = 0.1}}
  
  \item{\code{maxiter}}{An integer value indicating the maximum number of iterations before stopping the Genetic Algorithm. By default \code{maxiter = 1000}.}
\item{run}{An integer value indicating the number of generations without improvment in the best value of fitness fuction. \code{run = 100}}

\item{selection}{An \code{R} function performing the selection step in \code{ga}. See \code{\link{ga_Selection}} for the available functions. By default \code{selection = gareal_lsSelection}.}

\item{crossover}{An \code{R} function performing the crossover step in \code{ga}. See \code{\link{ga_Crossover}} for the available functions. By default \code{crossover = gareal_laCrossover}.}

\item{mutation}{An \code{R} function performing the mutation step in \code{ga}. See \code{\link{ga_Mutation}} for the available functions. By default \code{mutation = gareal_raMutation}.}

\item{parallel}{A logical value for parallel GA. By default \code{parallel = FALSE}.}

\item{numIslands}{An integer value indicating the number of island to be used in the Island Genetic Algorithm. By default \code{numIslands = 4}.}
\item{migrationRate}{A value indicating the percentage of migration between islands. By default \code{migrationRate = 0.1}.}

\item{migrationInterval}{An integer values indicating the number of generations with migrations. By default \code{migrationInterval = 10}.}
\item{optim}{A logical value for the local search. By default \code{optim = TRUE} for local search.}

\item{optimPoptim}{A value indicating the probability of local search. \code{optimPoptim = 0.05}.}

\item{optimPressel}{A value indicating the pressure selection. \code{optimPressel = 0.5}.}

\item{optimMethod}{A string indicating the general-purpose optimisation method to be used for local search. See \code{\link{optim}} fo the available algorithms. By default \code{optimMethod = "L-BFGS-B"}}

\item{optimMaxit}{An integer value indicating the number of iterations for the local search algorithm used. By default \code{optimMaxit = 100}.}

\item{orth}{A string indicating the method to orthogonalised the matrix with the basis. The available methods are the QR decomposition \code{"QR"} and Singular Value Decomposition \code{"SVD"}. By default \code{orth = "QR"}.}
}
The function globally sets the arguments for the current session of R. The default options are restored with a new R session. To change locally the options for only a single call to \code{ppgmmga} function,look at \code{options} argument in \code{\link{ppgmmga}}. 
}



\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}

}

\references{
Scrucca, L., Fop, M., Murphy, T. B., & Raftery, A. E. (2016). \emph{mclust 5: Clustering, classification and density estimation using gaussian finite mixture models.} The R journal, 8(1), 289.

Luca Scrucca (2013). \emph{GA: A Package for Genetic Algorithms in R.} Journal of Statistical Software, 53(4), 1-37. http://www.jstatsoft.org/v53/i04/

  Luca Scrucca (2017). \emph{On some extensions to GA package: hybrid optimisation, parallelisation and islands evolution.} The R Journal,
  9/1, 187-206. https://journal.r-project.org/archive/2017/RJ-2017-008
}

\seealso{
\code{\link{ppgmmga}}
}
\examples{
\dontrun{
ppgmmga.options()

# Print a single option

ppgmmga.options("popSize")

# Change option

ppgmmga.options("popSize" = 10)
ppgmmga.options("popSize")
}
}
\keyword{options}
