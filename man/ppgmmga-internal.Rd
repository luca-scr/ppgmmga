\name{ppgmmga-internal}
\alias{EntropyGMM}
\alias{EntropyUT}
\alias{EntropyVAR}
\alias{EntropySOTE}
\alias{EntropyMC}
\alias{EntropyGauss}
\alias{logsumexp}

\title{Internal ppgmmga functions}

\description{
Internal functions not intended to be called directly by users.
}
\usage{
EntropyGMM(data, pro, mean, sigma)
EntropyUT(G, pro, mean, sigma, d)
EntropyVAR(G, pro, mean, sigma, d)
EntropySOTE(data, G, pro, mean, sigma)
EntropyMC(G, pro, mean, sigma, d, nsamples)
EntropyGauss(S, d)
logsumexp(x)
}


\keyword{internal}
