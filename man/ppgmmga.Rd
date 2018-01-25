\name{ppgmmga}
\alias{ppgmmga}

\title{
Projection pursuit based on Gaussian mixtures and evolutionary algorithms for data visualisation}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

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
        monitor = TRUE, ...)
}

\arguments{
  \item{data}{
a \emph{n x p} matrix containing the data with rows correspond to observations and columns correspond to variables. }
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
a logical value set to \emph{TRUE} indicating whther or not the data are centred.
\verb{LS: potrebbe aver senso passare un vettore di medie da usare in scale ?}
}
  \item{scale}{
a logical value set to \emph{FALSE} indicating whther or not the data are scaled.
\verb{LS: potrebbe aver senso passare un vettore di std dev da usare in scale ?}
}
  \item{gmm}{
%%     ~~Describe \code{gmm} here~~
}
  \item{gatype}{
the type of the genetic algoritm used to maximised the Negentropy.
Possible values are:
  \tabular{ll}{
\code{"ga"} \tab for genetic algorithm. \cr
\code{"gaisl"} \tab for Island genetic algorithm. \cr}
}
  \item{opt}{
%%     ~~Describe \code{opt} here~~
}
  \item{monitor}{
a logical value set to \emph{TRUE} to show the evolution of the evolution of the optimisation with the genetic algorith.
}
  \item{\dots}{
additional arguments.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(iris)
X = iris[,-5]

# Unscented Transformation
pp = ppgmmga(data = X, d = 2, approx = "UT")

# Variational approximation
pp = ppgmmga(data = X, d = 2, approx = "VAR")


# Second order approximation
pp = ppgmmga(data = X, d = 2, approx = "SOTE")

# Example changing the size of the initial population in the general options

pp = ppgmmga(data = X, d = 2, opt = list("popSize" = 10))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
