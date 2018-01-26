\name{summary.ppgmmga}
\alias{summary.ppgmmga}
\alias{print.summary.ppgmmga}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary for projection pursuit based on Gaussian mixtures and evolutionary algorithms}
\description{
Summary method for class \code{ppgmmga}.}
\usage{
\method{summary}{ppgmmga}(object, ppgmmgaCheck = FALSE, ...)
                
\method{print}{summary.ppgmmga}(x, digits = getOption("digits"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{ppgmmga}.}
  \item{ppgmmgaCheck}{
%%     ~~Describe \code{ppgmmgaCheck} here~~
}
\item{x}{

an object of class \code{summary.ppgmmga}.

}

\item{digits}{

number of significant digits.
}
  \item{\dots}{
further arguments.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The summary function returns an object of class summary.ppgmmga which can be printed by the corresponding print method. A list with the information from the ppgmmga algorithm is returned.
}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}

\references{
%% ~put references to the literature/web site here ~
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

pp = ppgmmga(data = X, d = 2)
summ = summary(pp)
print(summ)
str(summ)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
