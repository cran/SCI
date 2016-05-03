\name{dist.start}
\alias{dist.start}
\alias{lmom.start}
\alias{mom.start}

\title{
  Rough estimates for parameters of selected distributions
}
\description{
  Produces rough parameter estimates for specific distributions
  (\code{distr}) that are useful as starting values for maximum
  likelihood estimation. 
}
\usage{
dist.start(x, distr,...)
lmom.start(x, distr = c("gamma","genlog","gev","gumbel",
          "lnorm","norm","pe3","weibull"),...)
mom.start(x, distr = c("gamma","gumbel","logis","lnorm","norm",
          "weibull"),...)
}

\arguments{
  \item{x}{
    \code{numeric} vector
  }
  \item{distr}{
    A character string \code{"name"} naming a distribution for which the
    corresponding density function (\code{dname}), the corresponding
    distribution function (\code{pname}) and the quantile function
    (\code{qname}) must be defined (see for example
    \code{\link{GammaDist}})
  }
  \item{\dots}{
    arguments passed to other functions, currently not used.
  }
}

\details{
  \code{lmom.start} uses L-moments for parameter estimation. In most
  cases it relies on functionality of the \pkg{lmomco}
  package. Currently available distributions are: \code{"\link{gamma}"},
  \code{"\link{genlog}"}, \code{"\link[evd:dgev]{gev}"},
  \code{"\link[evd:dgumbel]{gumbel}"}, \code{"\link[=Logistic]{logis}"}, 
  \code{"\link[=Lognormal]{lnorm}"}, \code{"\link{norm}"},
  \code{"\link{pe3}"}, \code{"\link[=Weibull]{weibull}"}.   

  \code{mom.start} uses moments (e.g. mean, standard deviation) for
  parameter estimation. Some estimates are precise, others only
  approximations that provide reasonable starting values. Currently
  available distributions are: \code{"\link{gamma}"},
  \code{"\link[evd:dgumbel]{gumbel}"}, \code{"\link[=Logistic]{logis}"},
  \code{"\link[=Lognormal]{lnorm}"}, \code{"\link{norm}"}, \code{"\link[=Weibull]{weibull}"}. 
  
  \code{dist.start} calls first \code{lmom.start} to estimate
  parameters. In case of failure \code{mom.start} is called, hopefully
  producing reasonable parameter estimates.
}

\value{
  named list, names correspond to distribution parameters. In case of
  failure, the same list with \code{NA} values is returned.
}


\author{
  Lukas Gudmundsson & James Stagge 
}


\examples{
lmom.start(rgamma(100,shape=0.5,rate=1),"gamma")

mom.start(rgamma(100,shape=0.5,rate=1),"gamma")

dist.start(rgamma(100,shape=0.5,rate=1),"gamma")
}
