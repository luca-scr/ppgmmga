\name{volume}
\alias{volume}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Volume of the data
}
\description{
Compute the volume of the data using the smallest hyperrectangle with sides parallel to the coordinate axes containing the data (\code{"box"}), the smallest hyperrectangle with sides parallel to the principal component axes (\code{"pc"}), and the Convex Hull of the data (\code{"ConvHull"})}
\usage{
volume(data,
       method = c("box", "pc", "ConvHull"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A numeric matrix of observations.
}
  \item{method}{A string indicating the method to build the region of the data. Possible values are:
  \tabular{ll}{
\code{"box"} \tab for the smallest hyperrectangle with sides parallel to the coordinate axes containing the data. \cr
\code{"pc"} \tab for the smallest hyperrectangle with sides parallel to the principal component axes. \cr
\code{"ConvHull"} \tab for the Convex Hull of the data.
}
The default is \code{"box"}.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the value of the volume.
}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}
\references{
%% ~put references to the literature/web site here ~
}

%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%5\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{

data("faithful")
data <- faithful

### Box method

volume(data, method = "box")

### Volume of the Convex Hull

volume(data, method = "ConvHull")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
