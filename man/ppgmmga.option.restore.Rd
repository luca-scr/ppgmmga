\name{ppgmmga.option.restore}
\alias{ppgmmga.option.restore}
\title{
Restore the default values for \code{ppgmmga} package}
\description{
Such function restores the default values used in \code{ppgmmga} package.}
\usage{
ppgmmga.option.restore()
}
\details{
Using the option function \code{ppgmmga.options} function to change options of \code{ppgmmga} package, the value are set until a new R session is start. Then, each time that \code{ppgmmga} function is called, the options are applied to the procedure. This function is implemented to restore the all the initial options without restart a new R session or manually changing them.}

\author{
Serafini A. \email{srf.alessio@gmail.com}

Scrucca L. \email{luca.scrucca@unipg.it}
}

\seealso{
\code{\link{ppgmmga.options}}
}
\examples{
ppgmmga.options("popSize")
ppgmmga.options("popSize" = 10)
ppgmmga.option.restore()
ppgmmga.options("popSize")

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
