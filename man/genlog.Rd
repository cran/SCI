\name{genlog}
\alias{genlog}
\alias{dgenlog}
\alias{pgenlog}
\alias{qgenlog}

\title{
  Generalized Logistic Distribution
}
\description{
  Density, distribution and quantile function of the generalized
  logistic distribution
}

\usage{
pgenlog(q, shape, scale, location)
dgenlog(x, shape, scale, location)
qgenlog(p, shape, scale, location)
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
  The functions of the \code{genlog} family are a reimplementation of
  the Generalized Logistic Distribution in the \pkg{lmomco}
  package, making the code compatible with the standard nomenclature for
  distribution in \R. The original functions in \pkg{lmomco} are 
  \code{\link{pdfglo}} (density function), \code{\link{quaglo}}
  (quantile function) and \code{\link{cdfglo}} (distribution function).  
}


\value{
  \code{dgenlog} gives the density (pdf), \code{pgenlog} gives
  the distribution function (cdf), and \code{qgenlog} gives the
  quantile function (inverse cdf).
}

\references{
  Asquith, W.H., 2013, lmomco -- L-moments, trimmed L-moments,                                                                      
  L-comoments, censored L-moments, and many distributions. R package                                                               
  version 1.7.8, \url{http://www.cran.r-project.org/package=lmomco}, Texas                                                               
  Tech University, Lubbock, Texas.  
}

\author{
  James Stagge & Lukas Gudmundsson
}


\examples{
dgenlog(1, shape=1, scale=2, location=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}

