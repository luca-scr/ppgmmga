\name{ppgmmga}
\alias{ppgmmga}
\alias{print.ppgmmga}

\title{
Projection pursuit based on Gaussian mixtures and evolutionary algorithms for data visualisation}
\description{

Cluster visualisation in multivariate data using an tailored projection pursuit algorithm. Gaussian mixture models and Genetic algorithms are employed to project the data and to underline cluster structures otherwise impossible to visualise.

TODO: quasi sempre non usare \verb{\emph} ma \verb{\code} \cr
Vedi \url{https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Marking-text}
e più in generale 
\url{https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-R-documentation-files}
}
\usage{
ppgmmga(data, 
        d, 
        approx = c("UT", "VAR", "SOTE"), 
        center = TRUE, 
        scale = FALSE, 
        gmm = NULL, 
        gatype = c("ga", "gaisl"), 
        opt = list(),
        seed = NULL, ...)
}

\arguments{
  \item{data}{
a \code{n x p} matrix containing the data with rows correspond to observations and columns correspond to variables. }
  \item{d}{
the dimension of the subspace onto which the data are projected and visualised.}
  \item{approx}{
the type of approximation to compute the Negentropy for GMM
model.
Possible values are:
  \tabular{ll}{
\code{"UT"} \tab for Unscented approximation. \cr
\code{"VAR"} \tab for Variational approximation. \cr
\code{"SOTE"} \tab for Second Order Taylor approximation. \cr
\code{"MC"} \tab Monte Carlo approximation. LS: non sembra essere previsto... \cr}
}
  \item{center}{
a logical value set to \code{TRUE} indicating whther or not the data are centred.
\verb{LS: potrebbe aver senso passare un vettore di medie da usare in scale ?}
}
  \item{scale}{
a logical value set to \code{FALSE} indicating whther or not the data are scaled.
\verb{LS: potrebbe aver senso passare un vettore di std dev da usare in scale ?}
}
  \item{gmm}{
An object of class \code{densityMclust}. This object contain an external estimate of the Gaussian mixture refered to the input data. The user can specify its own mixture estimates, providing a \code{densityMclust} object to pass to the function. See \code{Mclust} documentation for more details.}
  \item{gatype}{
the type of the genetic algoritm used to maximised the Negentropy.
Possible values are:
  \tabular{ll}{
\code{"ga"} \tab for genetic algorithm. \cr
\code{"gaisl"} \tab for Island genetic algorithm. \cr}
}
  \item{opt}{
A list of options for the density estimation and Genetic algorithms embedded in the function. This parameter collect all the important options to pass to \code{densityMclust} function of \code{Mclust} package for the density estimation, and to code{ga} function of \code{GA} package for genetic operator of the Genetic algorithm.}
  \item{seed}{an integer value with the random number generator state. It may be used to replicate the results of ppgmmga algorithm.}
  \item{\dots}{
additional arguments.
}
}
\details{
Projection pursuit (PP) is a features extraction method for analysing high-dimensional data with low-dimension projections, by maximising a projection index to find out the best orthogonal projections. A general PP procedure can be summarised in few steps: the data may be transformed, usually sphered, the PP index is chosen and the subspace dimension is fixed. Then, the PP index is optimised.

For clusters visualisation the Negentropy index has been considerd. Since such index requires an estimation of the underling data density, Gaussian mixture models (GMMs) has been used to approximate such density.  GMMs do not have a closed formula for the Negentropy and different closed formula approximations have been implemented. Genetic algorithms have been employed to maximised the approximated Negentropy respect to the system of basis in the desidered subspace.}
\value{
Returns an object of class \code{ppgmmga-class}. See \code{\link{ppgmmga-class}} for a description of the object.
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
181–188. IEEE.}
\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}


\seealso{
\code{\link{summary.ppgmmga},\code{\link{plot.ppgmmga}},\code{\link{ppgmmga-class}}}
}
\examples{
data(iris)
X <- iris[,-5]
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

# Changing the size of the initial population in the general options
pp4 <- ppgmmga(data = X, d = 2, opt = list("popSize" = 10))
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
