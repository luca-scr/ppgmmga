\name{ppgmmga}
\alias{ppgmmga}
\alias{print.ppgmmga}
\title{
Projection Pursuit based on Gaussian Mixtures and Evolutionary Algorithms for data visualisation}

\description{
Cluster visualisation of multivariate data using an tailored Projection Pursuit algorithm. Gaussian Mixture Models and Genetic Algorithms are employed to project the data underling cluster structures otherwise impossible to visualise.}
\usage{
ppgmmga(data, 
        d, 
        approx = c("UT", "VAR", "SOTE"), 
        center = TRUE, 
        scale = TRUE, 
        gmm = NULL, 
        gatype = c("ga", "gaisl"), 
        options = ppgmmga.options(),
        seed = NULL, 
        verbose = interactive(), ...)
}
\arguments{
  \item{data}{a \eqn{n x p} matrix containing the data with rows corresponding to observations and columns corresponding to variables. }
  \item{d}{
the dimension of the subspace onto which the data are projected and visualised.}
  \item{approx}{
the type of approximation to compute the Negentropy for GMMs.
Possible values are:
  \tabular{ll}{
\code{"UT"} \tab for Unscented approximation. \cr
\code{"VAR"} \tab for Variational approximation. \cr
\code{"SOTE"} \tab for Second Order Taylor approximation. \cr}
}
  \item{center}{a logical value set to \code{TRUE} indicating wheter or not the data are centred.}
  
  \item{scale}{a logical value set to \code{TRUE} indicating wheter or not the data are scaled.}
  
  \item{gmm}{an object of class \code{\link{densityMclust}}. This object contain an external estimate of the Gaussian mixture refered to the input data. The user can specify its own mixture estimates providing a \code{densityMclust} object to pass to the function. See \code{\link{Mclust}} documentation for more details.}
  
  \item{gatype}{the type of the genetic algoritm used to maximised the Negentropy.
Possible values are:
  \tabular{ll}{
\code{"ga"} \tab for genetic algorithm (\code{\link{ga}}). \cr
\code{"gaisl"} \tab for Island genetic algorithm \code{\link{gaisl}}. \cr}}
  
  \item{options}{a list of options containing all the important arguments to pass to \code{\link{densityMclust}} function of \code{\link{Mclust}} package, and to \code{\link{ga}} function of \code{GA} package. See \code{\link{ppgmmga.options}} for the available options. The \code{options} argument does not change the global options in \code{ppgmmga.options}, but only the options for the single call to \code{ppgmmga}.}
  
  \item{seed}{an integer value with the random number generator state. It may be used to replicate the results of ppgmmga algorithm.}
  
  \item{verbose}{a logical value controlling if the evolution of GA search is shown. By default is \code{TRUE} reporting the number of iteration, average and best fitness value.}
  
  \item{\dots}{
additional arguments.
}
}
\details{
Projection pursuit (PP) is a features extraction method for analysing high-dimensional data with low-dimension projections by maximising a projection index to find out the best orthogonal projections. A general PP procedure can be summarised in few steps: the data may be transformed, the PP index is chosen and the subspace dimension is fixed. Then, the PP index is optimised.

For clusters visualisation the Negentropy index has been considerd. Since such index requires an estimation of the underling data density, Gaussian mixture models (GMMs) has been used to approximate such density.  GMMs do not have a closed formula for the Negentropy and different closed formula approximations have been implemented. Genetic Algorithms have been employed to maximised the approximated Negentropy respect to the system of basis in the desidered subspace.
}
\value{
Returns an object of class \code{ppgmmga}. See \code{\link{ppgmmga-class}} for a description of the object.
}


\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}

}

\references{
Goldberger, J. and Aronowitz, H. (2005). \emph{A distance measure between gmms based on
the unscented transform and its application to speaker recognition}. In INTERSPEECH,
pages 1985–1988. Citeseer.

Hershey, J. R. and Olsen, P. A. (2007). \emph{Approximating the Kullback Leibler divergence
between Gaussian mixture models}. In 2007 IEEE International Conference on Acoustics,
Speech and Signal Processing-ICASSP’07, volume 4, pages IV–317. IEEE.

Huber, P. J. (1985). \emph{Projection pursuit}. The Annals of Statistics, pages 435–475.

Huber, M. F., Bailey, T., Durrant-Whyte, H., and Hanebeck, U. D. (2008). \emph{On entropy
approximation for Gaussian mixture random vectors}. In Multisensor Fusion and Integra-
tion for Intelligent Systems, 2008. MFI 2008. IEEE International Conference on, pages
181–188. IEEE.
Jones, M. C. and Sibson, R. (1987). \emph{What is projection pursuit? (with discussion).} Journal of the Royal Statistical Society. Series A (General), pages 1–37.

Kruskal, J. B. (1969). \emph{Toward a practical method which helps uncover the structure of a set of multivariate observations by finding the linear transformation which optimizes a new index of condensation.} In Milton, R. C. and Nelder, J. A., editors, Statistical Computation, pages 427–440. Academic Press, New York.
}

\seealso{
\code{\link{summary.ppgmmga}}, \code{\link{plot.ppgmmga}}, \code{\link{ppgmmga-class}}
}
\examples{
\dontrun{
require(mclust)
data(iris)
X = iris[,-5]
Class <- iris$Species

# Unscented Transformation
pp1 <- ppgmmga(data = X, d = 2, approx = "UT")  # default
summary(pp1, check = TRUE)
plot(pp1, Class)


# Variational approximation
pp2 <- ppgmmga(data = X, d = 2, approx = "VAR")
summary(pp2, check = TRUE)
plot(pp2, Class)

# Second order Taylor approximation
pp3 <- ppgmmga(data = X, d = 2, approx = "SOTE")
summary(pp3, check = TRUE)
plot(pp3, Class)

# Changing the initial population size
pp4 <- ppgmmga(data = X, d = 2, opt = list("popSize" = 10))
}
}

\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
