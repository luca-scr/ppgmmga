\name{ppgmmga}
\alias{ppgmmga}
\alias{print.ppgmmga}

\title{Projection pursuit based on Gaussian mixtures and evolutionary algorithms for data visualisation}

\description{
A Projection Pursuit (PP) method for dimension reduction seeking "interesting" data structures in low-dimensional projections. An approximated negentropy index is computed from the density estimated using Gaussian Mixture Models (GMMs). Then, the PP index is maximised by Genetic Algorithms (GAs) to find the optimal projection basis.
}

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
        verbose = interactive(), \dots)
}
\arguments{
  \item{data}{A \eqn{n x p} matrix containing the data with rows corresponding to observations and columns corresponding to variables. }
  \item{d}{An integer specifying the dimension of the subspace onto which the data are projected and visualised.}
  \item{approx}{A string specifying the type of approximation to use for computing the negentropy for GMMs. Possible values are:
  \tabular{ll}{
\code{"UT"} \tab for Unscented Trasformation approximation. \cr
\code{"VAR"} \tab for VARiational approximation. \cr
\code{"SOTE"} \tab for Second Order Taylor Expansion approximation.}
}
  \item{center}{A logical value indicating whether or not the data are centred. By default is set to \code{TRUE}.}
  
  \item{scale}{A logical value indicating whether or not the data are scaled. By default is set to \code{TRUE}.}
  
  \item{gmm}{An object of class \code{'densityMclust'} specifying a Gaussian mixture density estimate as returned by \code{\link{densityMclust}}.}
  
  \item{gatype}{A string specifying the type of genetic algoritm to be used to maximised the negentropy. Possible values are:
  \tabular{ll}{
\code{"ga"} \tab for simple genetic algorithm (\code{\link{ga}}). \cr
\code{"gaisl"} \tab for island genetic algorithm (\code{\link{gaisl}}). \cr}
}
  
  \item{options}{A list of options containing all the important arguments to pass to \code{\link{densityMclust}} function of the \pkg{mclust} package, and to \code{\link{ga}} function of the \pkg{GA} package. See \code{\link{ppgmmga.options}} for the available options. Note that by setting the \code{options} argument does not change the global options provided by \code{ppgmmga.options}, but only the options for a single call to \code{ppgmmga}.}
  
  \item{seed}{An integer value with the random number generator state. It may be used to replicate the results of ppgmmga algorithm.}
  
  \item{verbose}{A logical value controlling if the evolution of GA search is shown. By default is \code{TRUE} reporting the number of iteration, average and best fitness value.}
  
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{
Projection pursuit (PP) is a features extraction method for analysing high-dimensional data with low-dimension projections by maximising a projection index to find out the best orthogonal projections. A general PP procedure can be summarised in few steps: the data may be transformed, the PP index is chosen and the subspace dimension is fixed. Then, the PP index is optimised.

For clusters visualisation the negentropy index has been considerd. Since such index requires an estimation of the underling data density, Gaussian mixture models (GMMs) has been used to approximate such density.  GMMs do not have a closed formula for the Negentropy and different closed formula approximations have been implemented. Genetic Algorithms have been employed to maximised the approximated negentropy respect to the system of basis in the desidered subspace.
}
\value{
Returns an object of class \code{'ppgmmga'}. See \code{\link{ppgmmga-class}} for a description of the object.
}

\author{
Serafini A. \email{srf.alessio@gmail.com}\cr
Scrucca L. \email{luca.scrucca@unipg.it}
}

\references{
Scrucca, L. and Serafini, A. (2019) Projection pursuit based on Gaussian mixtures and evolutionary algorithms. *Journal of Computational and Graphical Statistics*, DOI: 10.1080/10618600.2019.1598871 
}

\seealso{
\code{\link{summary.ppgmmga}}, \code{\link{plot.ppgmmga}}, \code{\link{ppgmmga-class}}
}
\examples{
\dontrun{
data(iris)
X <- iris[,-5]
Class <- iris$Species

# 1-dimensional PPGMMGA

pp1D <- ppgmmga(data = X, d = 1)
summary(pp1D, check = TRUE)
plot(pp1D, Class)

# 2-dimensional PPGMMGA

## Unscented Transformation approximation (default)

pp2D_1 <- ppgmmga(data = X, d = 2, approx = "UT")
summary(pp2D_1, check = TRUE)
plot(pp2D_1, Class)

## VARiational approximation

pp2D_2 <- ppgmmga(data = X, d = 2, approx = "VAR")
summary(pp2D_2, check = TRUE)
plot(pp2D_2, Class)

## Second Order Taylor Expansion approximation

pp2D_3 <- ppgmmga(data = X, d = 2, approx = "SOTE")
summary(pp2D_3, check = TRUE)
plot(pp2D_3, Class)

# 3-dimensional PPGMMGA

pp3D <- ppgmmga(data = X, d = 3, approx = "UT")
summary(pp3D, check = TRUE)
plot(pp3D, Class)

# A rotating 3D plot can be obtained using:
# if(!require("msir")) install.packages("msir")
# msir::spinplot(pp3D$Z, markby = Class,
#                col.points = ggthemes::tableau_color_pal("Classic 10")(3))
}
}
