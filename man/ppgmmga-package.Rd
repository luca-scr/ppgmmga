\name{ppgmmga-package}
\alias{ppgmmga-package}
\alias{ppgmmga}
\docType{package}

\title{
Projection pursuit based on Gaussian mixtures and evolutionary
        algorithms for data visualisation and anomaly detection
}
\description{
An R packege implementing a projection pursuit algorithm based on an Finite Gaussian mixture for density estimation, and genetic algorithm to maximise the approximate Negentropy index for Gaussian mixture. The ppgmmga algorithm implemented in the package provides a method to visualise high-dimensional data. In adddition the package provide an anomaly detection algorithm based on the ppgmmga algorithm and Gaussian mixture with Uniform component to detect the anomalies in a dataset. 
}


\details{
For a quick intro to ppgmmga package see the \href{}{vignette}. 
%%Further details are provided in the papers referenced below.
}

\references{
%% ~put references to the literature/web site here ~
}
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

# Data visualisation with ppgmmga

pp = ppgmmga(data = iris[,-5], d = 2)
plot(pp)

# Anomaly detection

ppAD = ppgmmgaAD(data = iris[,-5])
plot(ppAD)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
