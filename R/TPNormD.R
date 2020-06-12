# Arash : 12/05/2020
#====================================================================
#         dtwopnorm
#====================================================================

dtwopnorm <- function(x, location = 0, scale = 1, skewpar = 0.5,
                 log.arg = FALSE) {
# Reference: Arash handnotes

  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
           na.rm = TRUE))
  stop("some parameters out of bound")

# Recycle the vectors to equal lengths
  LLL = max(length(x), length(location),
            length(scale), length(skewpar))
  if (length(x)        != LLL)        x = rep(x,        length.out = LLL)
  if (length(location) != LLL) location = rep(location, length.out = LLL)
  if (length(scale)    != LLL)    scale = rep(scale,    length.out = LLL)
  if (length(skewpar)  != LLL)  skewpar = rep(skewpar,  length.out = LLL)

  zedd <- (x - location) / scale
  log.s1 <-  -zedd^2 / (8 * skewpar^2)
  log.s2 <-  -zedd^2 / (8 * (1 - skewpar)^2)
  logdensity <- log.s1
  logdensity[zedd > 0] <- log.s2[zedd > 0]
  logdensity <- logdensity - log(scale) - 0.5 * log(2 * pi)

  if (log.arg) logdensity else exp(logdensity)
}



#=========================================================================
#         ptwopnorm
#=========================================================================


ptwopnorm <- function(q, location = 0, scale = 1, skewpar = 0.5) {
  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
          na.rm = TRUE))
    stop("some parameters out of bound")

# Reference: Arash handnotes

# Recycle the vectors to equal lengths
  LLL = max(length(q), length(location),
            length(scale), length(skewpar))
  if (length(q)        != LLL)        q = rep(q,        length.out = LLL)
  if (length(location) != LLL) location = rep(location, length.out = LLL)
  if (length(scale)    != LLL)    scale = rep(scale,    length.out = LLL)
  if (length(skewpar)  != LLL)  skewpar = rep(skewpar,  length.out = LLL)

  zedd <- (q - location) / scale
  s1 <- 2 * skewpar * pnorm(zedd, sd = 2 * skewpar)
  s2 <- skewpar + (1 - skewpar) *
        pgamma(zedd^2 / (8 * (1-skewpar)^2), 0.5)
  ans <- rep(0.0, length.out = length(zedd))
  ans[zedd <= 0] <- s1[zedd <= 0]
  ans[zedd >  0] <- s2[zedd >  0]
  ans
}


#=========================================================================
#         qtwopnorm
#=========================================================================


qtwopnorm <- function(p, location = 0, scale = 1, skewpar = 0.5){

# 20120328; this function has not been checked fully.

# posfun <- function(x) ifelse(x >  0, x, 0.0)
  posfun <- function(x) pmax(x, 0)
  lone   <- function(x) ifelse(x <= 1, x, 0.0)

  pp = p
  if (any(pp      <= 0 |
          pp      >= 1 |
          skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
             na.rm = TRUE))
  stop("some parameters out of bound")

# Recycle the vectors to equal lengths
  LLL = max(length(pp), length(location),
            length(scale), length(skewpar))
  if (length(pp)       != LLL)       pp = rep(pp,       length.out = LLL)
  if (length(location) != LLL) location = rep(location, length.out = LLL)
  if (length(scale)    != LLL)    scale = rep(scale,    length.out = LLL)
  if (length(skewpar)  != LLL)  skewpar = rep(skewpar,  length.out = LLL)


  qtpn <- rep(as.numeric(NA), length.out = LLL)
  qtpn <- qnorm(lone(pp / (2 * skewpar)), sd = 2 * skewpar)
  qtpn[pp > skewpar] <- ((1 - skewpar) *
                        sqrt(8 *
                             qgamma(posfun(pp - skewpar) / (1 - skewpar),
                                    0.5)))[pp > skewpar]
  location + qtpn * scale
}


#=========================================================================
#         rtwopnorm
#=========================================================================

rtwopnorm <- function(n, location = 0, scale = 1, skewpar = 0.5) {
  p = runif(n)
  rtpn <- qtwopnorm(p, location = location, scale = scale,
          skewpar = skewpar)
  rtpn
}





