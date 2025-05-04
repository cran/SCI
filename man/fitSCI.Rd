\name{fitSCI}
\alias{fitSCI}
\alias{fitSCI.default}
\alias{transformSCI}
\alias{transformSCI.default}


\title{
  Standardized Climate Index (SCI)
}
\description{
  \code{fitSCI} identifies parameters for the Standardized Climate Index
  (SCI) transformation. \code{transformSCI} applies the transformation
}
\usage{
fitSCI(x, ...)
\method{fitSCI}{default}(x, first.mon, time.scale, distr, p0, 
       p0.center.mass=FALSE, scaling=c("no","max","sd"),mledist.par =  list(),
       start.fun = dist.start, start.fun.fix = FALSE, warn = TRUE, ...)

transformSCI(x, ...)
\method{transformSCI}{default}(x, first.mon, obj, sci.limit = Inf, warn=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    \code{numeric} vector, representing a monthly univariate time series.
  }
  \item{first.mon}{
    value in 1:12 indicating the month of the first element of x
  }
  \item{time.scale}{
    The time scale (\code{integer}) of  the SCI calculation. The time
    scale is the window length of an backward looking running mean.
  }
  \item{distr}{
    A character string \code{"name"} naming a distribution for which the
    corresponding density function (\code{dname}), the corresponding
    distribution function (\code{pname}) and the quantile function
    (\code{qname}) must be defined (see for example
    \code{\link{GammaDist}})
  }
  \item{p0}{
    if TRUE, model Probability of zero (precipitation) months is
    modeled with a mixed distribution as \eqn{D(x) = p_0 +
    (1-p_0)G(x)}{D(x) = p0 + (1-p0)G(x)}, where  \eqn{G(x) > 0} is the 
    reference distribution (e.g. Gamma) and \eqn{p0} is the probability of a 
    zero (precipitation) month. 
  }
  \item{p0.center.mass}{
    If TRUE, the Probability of zero (precipitation) is estimated using
    a "center of mass" estimate based on the Weibull plotting position
    function (see details). Only applies if \code{p0=TRUE}. 
  }  
  \item{scaling}{
    Indicates whether to do some scaling of \code{x} prior to parameter
    identification. \code{"no"} (the default) indicates no scaling.
    \code{"max"} indicates scaling by the maximum of \code{x}, such that
    \code{x <- x/max(x,na.rm=TRUE)}. \code{"sd"} stands for scaling by
    the standard deviation. Scaling can stabilize parameter estimation.
  }  
  \item{mledist.par}{
    named \code{list} that can be used to pass parameters to
    \code{\link[fitdistrplus]{mledist}} in package \pkg{fitdistrplus}.
  }  
  \item{start.fun}{
    Function with arguments \code{x} and \code{distr} estimating initial
    parameters of the function \code{distr} for each month. The function
    should return a named list corresponding to the parameters of
    \code{distr}. (See also \code{\link{dist.start}})
  }
  \item{start.fun.fix}{
    \code{logical} argument, indicating if parameter estimates of
    \code{start.fun} should be used if maximum likelihood estimation
    breaks down. This stabilizes the implementation but can introduce
    biases in the resulting SCI.
  }  
  \item{obj}{
    an object of class \code{fitSCI}, output from \code{\link{fitSCI}}.
  }  
  \item{sci.limit}{
    Truncate absolute values of SCI that are lage than
    \code{sci.limit}. See details.
  }  
  \item{warn}{
    Issue warnings if problems in parameter estimation occur.
  }     
  \item{\dots}{
    further arguments passed to methods
  }
}
\details{
  \code{fitSCI} estimates the parameters for transforming a
  meteorological and environmental time series to a Standardized
  Climate Index (SCI). \code{transformSCI} applies the
  standardisation. Typical SCI are the Standardized Precipitation
  Index (SPI), the Standardized Runoff Index (SRI) or the Standardized
  Precipitation Evapotranspiration Index (SPEI). 
  
  To reduce biases in the presence of many zero (precipitation) events,
  the probability of these events (\eqn{p_0}{p0}) can be estimated using
  a "center of mass" estimate based on the Weibull plotting position
  function (\code{p0.center.mass=TRUE}). Following Stagge et al. (2014)
  the probability of zero events is then estimated as \eqn{p_0 =
  \frac{n_p}{n + 1}}{p0 = (n_p)/(n + 1)}, where \eqn{n_p}{np} refers to
  the number of zero events and \eqn{n} is the sample size. The
  resulting mixed distribution used fro SCI transformation is then   
  %%
  \deqn{D(x) = \left\{
    \begin{array}{l l}
    p_0 + (1-p_0) G(x) & \quad \mbox{if } x > 0 \\
    \frac{n_p + 1}{2(n+1)} & \quad \mbox{if } x = 0    
    \end{array} \right.
  }{g(x) = if(x > 0) p0 + (1 - p0)G(x)
    else if(x == 0) (np + 1)/(2(n + 1))}
  %%
  where \eqn{G(x)>0} is a model (e.g. gamma) distribution.
  
  Uncertainty in distribution parameters can cause unrealistically large
  (small) SCI values if values in \code{x} exceed the values used for
  parameter estimation (see \code{fitSCI}). Therefore
  \code{transformSCI} allows for a truncation of the SCI series such
  that \code{abs(sci)<=sci.limit}. The truncation can be disabled by
  setting \code{sci.limit=Inf}.
}
\value{
  \code{fitSCI} returns an object of class \code{"fitSCI"} with the
  following components:
  \item{dist.para}{
    A column \code{matrix} containing the parameters of distribution
    \code{distr} for each month. Row names correspond to the distribution
    parameters. If \code{p0=TUE} an additional row named \code{P0} is
    introduced, indicating the probability of zero (precipitation)
    events. 
  }
  \item{dist.para.flag}{
    an vector indicating possible issues occurring throughout parameter
    estimation. Possible values are: 0. no problems occurred;
    1. starting values could not be estimated; 2. \code{\link[fitdistrplus]{mledist}}
    crashed with unknown error; 3. \code{\link[fitdistrplus]{mledist}} did not
    converge; 4. all values in this month are \code{NA}; 5. all values
    in this month are constant, distribution not defined.
  }
  \item{time.scale}{
    The time scale (\code{integer}) of  the SCI calculation. 
  }
  \item{distr}{
     A character string \code{"name"} naming a distribution used
  }
  \item{p0}{
    \code{logical} indicating whether probability of zero
    (precipitation) events is estimated separately.
  }
  \item{p0.center.mass}{
        \code{logical} indicating whether probability of zero
	(precipitation) events is estimated using the "centre of mass"
        estimator (see Stagge et al. (2014) for details).

  }
  \item{scaling}{
    \code{numeric} value that has been used to scale \code{x} (see
    argument \code{scaling}). A value of 1 results from
    \code{scaling="no"}, other values are the maximum value or the
    standard deviation of \code{x}, depending on the choice of the
    parameter \code{scaling}.
  }   
  \item{call}{
    the function call
  }  

  \code{transformSCI} returns a \code{numeric} vector
  containing the SCI, having values of the standard normal
  distribution. 
}

\references{
  Stagee, J.H. ; Tallaksen, L.M.; Gudmundsson, L.; van Loon, A.;
  Stahl, K.: Candidate Distributions for Climatological Drought Indices  (SPI and SPEI),
  2015, International Journal of Climatology, 35, 4027-4040, doi:10.1002/joc.4267.

  Stagee, J.H. ; Tallaksen, L.M.; Gudmundsson, L.; van Loon, A.;
  Stahl, K.: Response to comment on "Candidate Distributions for Climatological
  Drought Indices (SPI and SPEI)", 2016, International Journal of
  Climatology, 36, 2132-2138, doi:10.1002/joc.4564.
  
  McKee, T.; Doesken, N. & Kleist, J.: The relationship of drought
  frequency and duration to time scales Preprints, 8th Conference on
  Applied Climatology, 1993, 179-184.
  
  Shukla, S. & Wood, A. W.: Use of a standardized runoff index for
  characterizing hydrologic drought Geophysical Research Letters, 2008,
  35, L02405.
  
  Vicente-Serrano, S. M.; Begueria, S. & Lopez-Moreno, J. I.: A Multiscalar
  Drought Index Sensitive to Global Warming: The Standardized
  Precipitation Evapotranspiration Index J. Climate, Journal of Climate,
  American Meteorological Society, 2009, 23, 1696-1718.
}

\author{
  Lukas Gudmundsson & James Stagge
}

\note{
  This function is intended to be used together with
  \code{\link{transformSCI}}. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{dist.start}}
}

\examples{
##
## generate artificial data
##
set.seed(101)
n.years <- 60
date <- rep(1:n.years,each=12) + 1950 + rep((0:11)/12,times=n.years)
## Precipitation
PRECIP <- (0.25*sin( 2 * pi * date) + 0.3)*rgamma(n.years*12, shape = 3, scale = 1)
PRECIP[PRECIP<0.1] <- 0
## Potential Evapotranspiration
PET <- 0.5*sin( 2 * pi * date)+1.2+rnorm(n.years*12,0,0.2)
## display test data
matplot(date,cbind(PRECIP,PET),t=c("h","l"),col=c("blue","red"),lty=1)
legend("topright",legend=c("PRECIPitation","temperature"),fill=c("blue","red"))

##
## example SPI
##
spi.para <- fitSCI(PRECIP,first.mon=1,distr="gamma",time.scale=6,p0=TRUE)
spi.para
spi <- transformSCI(PRECIP,first.mon=1,obj=spi.para)
plot(date,spi,t="l")

##
## effect of time.scale on SPI
##
spi.1.para <- fitSCI(PRECIP,first.mon=1,time.scale=1,distr="gamma",p0=TRUE)
spi.12.para <- fitSCI(PRECIP,first.mon=1,time.scale=12,distr="gamma",p0=TRUE)
spi.1 <- transformSCI(PRECIP,first.mon=1,obj=spi.1.para)
spi.12 <- transformSCI(PRECIP,first.mon=1,obj=spi.12.para)
matplot(date,cbind(spi.1,spi.12),t="l",lty=1,col=c("red","blue"),lwd=c(1,2))
legend("topright",legend=c("time.scale=1","time.scale=12"),fill=c("red","blue"))

##
## example SPEI
##
if(require(evd)){
    spei.para <- fitSCI(PRECIP-PET,first.mon=1,time.scale=6,distr="gev",p0=FALSE)
    spei <- transformSCI(PRECIP-PET,first.mon=1,obj=spei.para)
    plot(date,spei,t="l")
}

##
## effect of changing different distribution for SPEI computation
##
spei.genlog.para <- fitSCI(PRECIP-PET,first.mon=1,time.scale=6,distr="genlog",p0=FALSE)
spei.genlog <- transformSCI(PRECIP-PET,first.mon=1,obj=spei.genlog.para)
if(require(evd)){lines(date,spei.genlog, col="red")} else {plot(date,spei.genlog,t="l")}
## in this case: only limited effect.
## generally: optimal choice of distribution: user responsibility.

##
## use a 30 year reference period for SPI parameter estimation
##
sel.date <- date>=1970 & date<2000
spi.ref.para <- fitSCI(PRECIP[sel.date],first.mon=1,distr="gamma",time.scale=6,p0=TRUE)
## apply the the parameters of the reference period to all data
## also outside the reference period
spi.ref <- transformSCI(PRECIP,first.mon=1,obj=spi.ref.para)
plot(date,spi.ref,t="l",col="blue",ylim=c(-5,5),lwd=2)
lines(date[sel.date],spi.ref[sel.date],col="red",lwd=3)
legend("bottom",legend=c("reference period","extrapolation period"),fill=c("red","blue"),
       horiz=TRUE)

##
## use "start.fun.fix" in instances where maximum likelyhood estimation fails
##
## force failure of maximum likelyhood estimation by adding "strange" value
## a warning should be issued
xx <- PRECIP-PET; xx[300] <- 1000
spei.para <- fitSCI(xx,first.mon=2,time.scale=1,p0=FALSE,distr="gev")
spei.para$dist.para
## use start.fun, usually ment for estimating inital values for
## parameter optimisation if maximum likelihood estimation fails
spei.para <- fitSCI(xx,first.mon=2,time.scale=1,p0=FALSE,distr="gev",
                    start.fun.fix=TRUE)
spei.para$dist.para

##
## usage of sci.limit to truncate unrealistic SCI values
## 
PRECIP.mod <- PRECIP
PRECIP.mod[300] <- 100 ## introduce spuriously large value
spi.mod.para <- fitSCI(PRECIP.mod,first.mon=1,time.scale=3,p0=TRUE,distr="gamma")
plot(transformSCI(PRECIP.mod,first.mon=1,obj=spi.mod.para,sci.limit=Inf),
     t="l",col="blue",lwd=2)
lines(transformSCI(PRECIP.mod,first.mon=1,obj=spi.mod.para,sci.limit=4),col="red")

##
## how to modify settings of function "mledist" used for parameter identification
## 
## identify parameters with standard settings
spi.para <- fitSCI(PRECIP,first.mon=1,distr="gamma",time.scale=6,p0=TRUE)
## add lower and upper limits for parameter identification
lower.lim <- apply(spi.para$dist.para,1,min) - 0.5*apply(spi.para$dist.para,1,sd)
upper.lim <- apply(spi.para$dist.para,1,max) + 0.5*apply(spi.para$dist.para,1,sd)
spi.para.limit <- fitSCI(PRECIP,first.mon=1,distr="gamma",time.scale=6,p0=TRUE,
                         mledist.par=list(lower=lower.lim, upper=upper.lim))

##
## how to write an own start.fun
## (required if distributions not mentioned in "dist.start" are used)
##
## function with same arguments as "dist.start" 
my.start <- function(x,distr="gamma"){
### code based on "mmedist" in package "fitdistrplus"
    ppar <- try({
        n <- length(x)
        m <- mean(x)
        v <- (n - 1)/n * var(x)
        shape <- m^2/v
        rate <- m/v
        list(shape = shape, rate = rate)},TRUE)
    if (class(ppar) == "try-error") ## function has to be able to return NA parameters
        ppar <- list(shape = NA, rate = NA)
    return(ppar)
}
my.start(PRECIP)
spi.para <- fitSCI(PRECIP,first.mon=1,time.scale=6,p0=TRUE,
                   distr="gamma",start.fun=my.start)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}



