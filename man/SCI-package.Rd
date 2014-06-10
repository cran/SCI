\name{SCI-package}
\alias{SCI-package}
\alias{SCI}
\docType{package}
\title{
  Standardized Climate Indices such as SPI, SRI or SPEI
}

\description{
  Functions for generating Standardized Climate Indices (SCI). SCI is a
  transformation of (smoothed) climate (or environmental) time series
  that  removes seasonality and forces the data to take values of the
  standard normal distribution. SCI was originally developed for
  precipitation. In this case it is known as the Standardized
  Precipitation Index (SPI).
}
\details{
\tabular{ll}{
  Package: \tab SCI\cr
  Type: \tab Package\cr
  Version: \tab 1.0-1\cr
  Date: \tab 2014-06-10\cr
  License: \tab GPL (>= 2)\cr
}

}
\author{
  Lukas Gudmundsson & James Stagge
  
  Maintainer: Lukas Gudmundsson <lukas.gudmundsson@env.ethz.ch>
}

\references{
  Stagee, J.H. ; Tallaksen, L.M.; Gudmundsson, L.; van Loon, A.;
  Stahl, K. Candidate Distributions for Climatological Drought Indices
  (SPI and SPEI), 2014, submitted. 
}

\keyword{ package }

\examples{
## create artificial data, resembling precipitation
set.seed(101)
n.years <- 60
date <- rep(1:n.years,each=12) + 1950 + rep((0:11)/12,times=n.years)
PRECIP <- (0.25*sin( 2 * pi * date) + 0.3)*rgamma(n.years*12, shape = 3, scale = 1)
PRECIP[PRECIP<0.1] <- 0

## apply SCI transformation
spi.para <- fitSCI(PRECIP,first.mon=1,time.scale=6,distr="gamma",p0=TRUE)
spi <- transformSCI(PRECIP,first.mon=1,obj=spi.para)
plot(date,spi,t="l")
}
