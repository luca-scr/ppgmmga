\name{ppgmmga.option.restore}
\alias{ppgmmga.option.restore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Restore the default values for \code{ppgmmga} package}
\description{
Restore the default values for use with \code{ppgmmga} package.}
\usage{
ppgmmga.option.restore()
}
%- maybe also 'usage' for other objects documented here.
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
ppgmmga.options("popSize")
ppgmmga.options("popSize" = 10)
ppgmmga.option.restore()
ppgmmga.options("popSize")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
