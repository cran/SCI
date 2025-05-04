lmom.start <- function(x,distr=c("gamma","genlog","gev","gumbel",
                                 "lnorm","norm","pe3","weibull"),...){
  ### estimates parameters of distributions using L-moments
  ###
  ### x : data vector
  ### distr : A character string "name" naming a distribution for which
  ### the corresponding density function ('dname'), the
  ### corresponding distribution function ('pname') and the
  ### quantile function ('qname') must be defined (see for example 'GammaDist'
  distr <- match.arg(distr)
  
  
  x.lmom <- try(lmom.ub(x),TRUE)
  if(distr == "gamma"){
    ppar <- tryCatch({
      result <- pargam(x.lmom, checklmom = FALSE)
      setNames(result$para, c("shape", "rate"))
    }, error = function(e) {
      c(shape = NA, rate = NA)
    })
    
  } else if(distr == "genlog"){
    ppar <- tryCatch({
      k.est <- -x.lmom$TAU3
      a.est <- x.lmom$L2 / (gamma(1 + k.est) * gamma(1 - k.est))
      e.est <- x.lmom$L1 + (x.lmom$L2 - a.est) / k.est
      c(shape = k.est, scale = a.est, location = e.est)
    }, error = function(e) {
      c(shape = NA, scale = NA, location = NA)
    })
    
    
  } else if(distr == "gev"){
    if(requireNamespace("evd", quietly = TRUE)){
      
      ppar <- tryCatch({
        result <- pargev(x.lmom, checklmom = FALSE)
        para <- setNames(result$para,c("loc","scale","shape"))
        para["shape"] <- -para["shape"]
        para
      }, error = function(e) {
        c(loc = NA, scale = NA, shape = NA)
      })
      
      
    } else {
      stop("package 'evd' needed fro 'distr=gev'")
    }
  } else if(distr == "gumbel"){
    if(requireNamespace("evd", quietly = TRUE)){
      
      ppar <- tryCatch({
        result <- pargum(x.lmom, checklmom = FALSE)
        setNames(result$para,c("loc","scale"))
      }, error = function(e) {
        c(loc = NA, scale = NA)
      })
      
      
      
    } else {
      stop("package 'evd' needed fro 'distr=gumbel'")
    }
  } else if(distr == "lnorm"){
    # 
    ppar <- tryCatch({
      result <- parln3(x.lmom, checklmom = FALSE)
      setNames(result$para[c("mulog", "sigmalog")], c("meanlog", "sdlog"))
    }, error = function(e) {
      c(meanlog = NA, sdlog = NA)
    })
    # 
    
  } else if(distr == "norm"){
    ppar <- tryCatch({
      result <- parnor(x.lmom, checklmom = FALSE)
      setNames(result$para, c("mean", "sd"))
    }, error = function(e) {
      c(mean = NA, sd = NA)
    })
    
    
    
  } else if(distr == "pe3"){
    ppar <- tryCatch({
      estimate <- parpe3(x.lmom, checklmom = FALSE)
      shape <- estimate$para[[3]]
      # Clamp shape parameter within [-1.75, 1.75]
      shape <- max(min(shape, 1.75), -1.75)
      c(
        shape = shape,
        scale = estimate$para[[2]],
        location = estimate$para[[1]]
      )
    }, error = function(e) {
      c(shape = NA, scale = NA, location = NA)
    })
    
    
  } else if(distr == "weibull"){
    ppar <- tryCatch({
      result <- parwei(x.lmom, checklmom = FALSE)
      setNames(result$para[c("delta", "beta")], c("shape", "scale"))
    }, error = function(e) {
      c(shape = NA, scale = NA)
    })
    
  }
  
  
  ppar <- as.list(ppar)
  return(ppar)
}

mom.start <- function(x,distr=c("gamma","gumbel","logis","lnorm","norm",
                                "weibull"),...){
  ### estimates parameters of distributions using moments
  ###
  ### x : data vector
  ### distr : A character string "name" naming a distribution for which
  ### the corresponding density function ('dname'), the
  ### corresponding distribution function ('pname') and the
  ### quantile function ('qname') must be defined (see for example 'GammaDist'
  distr <- match.arg(distr)
  if(distr %in% c("norm","lnorm","gamma","logis")){
    ppar <- tryCatch({
      result <- fitdistrplus::mmedist(x, distr)
      result$estimate
    }, error = function(e) {
      switch(distr,
             norm = c(mean = NA, sd = NA),
             lnorm = c(meanlog = NA, sdlog = NA),
             gamma = c(shape = NA, rate = NA),
             logis = c(location = NA, scale = NA),
             stop("Unsupported distribution type")
      )
    })
    
    
  } else if(distr == "gumbel"){
    if(requireNamespace("evd", quietly = TRUE)){
      ppar <- tryCatch({
        x <- x[x > 0]  # Filter out non-positive values
        m <- mean(x)  # Compute mean
        sd.data <- sd(x)  # Compute standard deviation
        scale <- (sd.data * (6) ^ 0.5) / pi  # Scale calculation
        location <- m - 0.5772 * scale  # Location calculation
        c(loc = location, scale = scale)
      }, error = function(e) {
        c(location = NA, scale = NA)
      })
    } else {
      stop("package 'evd' needed fro 'distr=gumbel'")
    }
  } else if(distr == "weibull"){
    ## adapted from 'mledist' in package 'fitdistrplus'
    ppar <- tryCatch({
      m <- mean(log(x))  # Mean of the log of x
      v <- var(log(x))   # Variance of the log of x
      shape <- 1.2 / sqrt(v)  # Shape calculation
      scale <- exp(m + 0.572 / shape)  # Scale calculation
      c(shape = shape, scale = scale)
    }, error = function(e) {
      c(shape = NA, scale = NA)
    })
    
  }
  ppar <- as.list(ppar)
  return(ppar)
}



dist.start <- function(x,distr,...){
  ### estimates starting values for distributions
  ###
  ### first 'lmom.start' is called. In case of failure of 'lmom.start'
  ### an attempt is undertaken to use 'mom.start'
  ###
  ### x : data vector
  ### distr : A character string "name" naming a distribution for which
  ### the corresponding density function ('dname'), the
  ### corresponding distribution function ('pname') and the
  ### quantile function ('qname') must be defined (see for example 'GammaDist'
  
  
  par1 <- try(suppressWarnings(lmom.start(x=x,distr=distr,...)),silent=TRUE)
  # if(class(par1)!="try-error"&all(is.finite(unlist(par1)))){
  if((!inherits(par1,"try-error"))&all(is.finite(unlist(par1)))){
    return(par1)
  } else {
    par2 <- try(suppressWarnings(mom.start(x=x,distr=distr,...)),silent=TRUE)
    if(inherits(par2,"try-error")){
      if(!inherits(par1,"try-error")){
        return(par1)
      } else {
        stop("parameters of distribution",distr,"could not be identifyed")
      }
    } else {
      return(par2)
    }
  }
  
}

