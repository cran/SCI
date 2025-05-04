\name{pe3}
\alias{pe3}
\alias{dpe3}
\alias{ppe3}
\alias{qpe3}

\title{
  Pearson Type III distribution
}
\description{
  Density, distribution and quantile function of the Pearson Type III
  distribution
}
\usage{
dpe3(x, shape, scale, location)
ppe3(q, shape, scale, location)
qpe3(p, shape, scale, location)
}

\arguments{
  \item{x, q}{
    vector of quantiles.
  }
  \item{p}{
    vector of probabilities.
  }
  \item{shape}{
    shape parameter
  }
  \item{scale}{
    scale parameter
  }
  \item{location}{
    location parameter
  }

}

\details{
  The functions of the \code{pe3} family are a reimplementation of
  the Pearson Type III Distribution in the \pkg{lmomco}
  package, making the code compatible with the standard nomenclature for
  distributions in \R. The original functions in \pkg{lmomco} are 
  \code{\link[lmomco]{pdfpe3}} (density function), \code{\link[lmomco]{quape3}}
  (quantile function) and \code{\link[lmomco]{cdfpe3}} (distribution function).  
  
}


\value{
  \code{dpe3} gives the density (pdf), \code{ppe3} gives
  the distribution function (cdf), and \code{qpe3} gives the
  quantile function (inverse cdf).
}

\references{
  Asquith, W.H., 2013: lmomco -- L-moments, trimmed L-moments,                                                                      
  L-comoments, censored L-moments, and many distributions. R package                                                               
  version 1.7.8 , Tech University, Lubbock, Texas.  
}

\author{
James Stagge & Lukas Gudmundsson
}


\examples{
dpe3(1, shape=1, scale=2, location=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
