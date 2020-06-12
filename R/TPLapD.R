

### Two-piece Laplace (twoplaplace) family ###

#=========================================================================
#         dtwoplaplace
#=========================================================================

dtwoplaplace <- function(x, location = 0, scale = 1, skewpar = 0.5,
                 log.arg = FALSE) {
 # Reference: Arash handnotes
   #if (any(
   #skewpar <= 0 |
   #      skewpar >= 1 |
    #      scale   <= 0 ,
    #       na.rm = TRUE)
    #       )
   # stop("some parameters out of bound")
 # Recycle the vectors to equal lengths
  LLL = max(length(x), length(location), length(scale),
            length(skewpar))
  if (length(x) != LLL) x = rep(x, length = LLL)
  if (length(location) != LLL) location = rep(location, length = LLL)
  if (length(scale) != LLL) scale = rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar = rep(skewpar, length = LLL)

  zedd <- (x - location) / scale
    s1 <-  exp(zedd / (2 * skewpar))
    #sech2(zedd / (4 * skewpar))
    s2 <- exp(-zedd/(2 * (1 - skewpar)))
    #sech2(zedd / (4 * (1 - skewpar)))
    dens <- s1
    dens[zedd > 0] <- s2[zedd > 0]
    dens <- 1 / (2 * scale) * dens

  if (log.arg) log(dens) else dens
}
#=========================================================================
#         ptwoplaplace
#=========================================================================
ptwoplaplace <- function(q, location = 0, scale = 1, skewpar = 0.5) {

     if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
          na.rm = TRUE))
    stop("some parameters out of bound")
   # Reference: Arash hand notes
    zedd <- (q - location) / scale
    p2 <- (1-skewpar)
  s1 <-  skewpar * (exp(zedd/ (2 * skewpar) - 0))
  s2 <- skewpar + p2 * (1 - exp(-zedd/ (2 * p2)))
 ans <- rep(0.0, length(zedd))
 ans[zedd <= 0] <- s1[zedd <= 0]
 ans[zedd > 0] <- s2[zedd > 0]
  ans
}
#=========================================================================
#         qtwoplaplace
#=========================================================================
qtwoplaplace <- function(p, location = 0, scale = 1, skewpar = 0.5){

  pp = p
  if (any(pp      <= 0 |
          pp      >= 1 |
          skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
             na.rm = TRUE))
    stop("some parameters out of bound")
    # Recycle the vectors to equal lengths
  LLL = max(length(pp), length(location), length(scale),
            length(skewpar))
  if (length(pp) != LLL) pp = rep(pp, length = LLL)
  if (length(location) != LLL) location = rep(location, length = LLL)
  if (length(scale) != LLL) scale = rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar = rep(skewpar, length = LLL)

   QL = pp / skewpar; QL = ifelse(QL <=0, 0.1, QL)
   QR = (1-skewpar) / (1 - pp); QR = ifelse(QR <= 0, 0.1, QR)

   qtp.L <- (2 * skewpar * log(QL)) * (pp <= skewpar)
   qtp.R <- (2 * (1 - skewpar) * log(QR)) * (pp > skewpar)

   qtp <- qtp.L + qtp.R
   qtp * scale + location
  }

#=========================================================================
#         rtwoplaplace
#=========================================================================

rtwoplaplace <- function(n, location = 0, scale = 1, skewpar = 0.5) {

  qtwoplaplace(p = runif(n), location = location, scale = scale,
               skewpar = skewpar)
}


     ### Two-piece Laplace (tpl) family ###

#=========================================================================
#         dtpl
#=========================================================================

dtpl <- function(x, location = 0, scale = 1, skewpar = 0.5,
                 log.arg = FALSE) {
 # Reference: Arash handnotes
   #if (any(
   #skewpar <= 0 |
   #      skewpar >= 1 |
    #      scale   <= 0 ,
    #       na.rm = TRUE)
    #       )
   # stop("some parameters out of bound")
 # Recycle the vectors to equal lengths
  LLL = max(length(x), length(location), length(scale),
            length(skewpar))
  if (length(x) != LLL) x = rep(x, length = LLL)
  if (length(location) != LLL) location = rep(location, length = LLL)
  if (length(scale) != LLL) scale = rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar = rep(skewpar, length = LLL)

  zedd <- (x - location) / scale
    s1 <-  exp(zedd / (2 * skewpar))
    #sech2(zedd / (4 * skewpar))
    s2 <- exp(-zedd/(2 * (1 - skewpar)))
    #sech2(zedd / (4 * (1 - skewpar)))
    dens <- s1
    dens[zedd > 0] <- s2[zedd > 0]
    dens <- 1 / (2 * scale) * dens

  if (log.arg) log(dens) else dens
}
#=========================================================================
#         ptpl
#=========================================================================
ptpl <- function(q, location = 0, scale = 1, skewpar = 0.5) {

     if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
          na.rm = TRUE))
    stop("some parameters out of bound")
   # Reference: Arash hand notes
    zedd <- (q - location) / scale
    p2 <- (1-skewpar)
  s1 <-  skewpar * (exp(zedd/ (2 * skewpar) - 0))
  s2 <- skewpar + p2 * (1 - exp(-zedd/ (2 * p2)))
 ans <- rep(0.0, length(zedd))
 ans[zedd <= 0] <- s1[zedd <= 0]
 ans[zedd > 0] <- s2[zedd > 0]
  ans
}
#=========================================================================
#         qtpl
#=========================================================================
qtpl <- function(p, location = 0, scale = 1, skewpar = 0.5){

  pp = p
  if (any(pp      <= 0 |
          pp      >= 1 |
          skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
             na.rm = TRUE))
    stop("some parameters out of bound")
    # Recycle the vectors to equal lengths
  LLL = max(length(pp), length(location), length(scale),
            length(skewpar))
  if (length(pp) != LLL) pp = rep(pp, length = LLL)
  if (length(location) != LLL) location = rep(location, length = LLL)
  if (length(scale) != LLL) scale = rep(scale, length = LLL)
  if (length(skewpar) != LLL) skewpar = rep(skewpar, length = LLL)

   QL = pp / skewpar; QL = ifelse(QL <=0, 0.1, QL)
   QR = (1-skewpar) / (1 - pp); QR = ifelse(QR <= 0, 0.1, QR)

   qtpL <- (2 * skewpar * log(QL)) * (pp <= skewpar)
   qtpR <- (2 * (1 - skewpar) * log(QR)) * (pp > skewpar)

   qtp <- qtpL + qtpR
   qtp * scale + location
  }

  qtpl((1:9)/10, skew = .8)
#=========================================================================
#         rtpl
#=========================================================================

rtpl <- function(n, location = 0, scale = 1, skewpar = 0.5) {

  qtpl(p = runif(n), location = location, scale = scale, skewpar = skewpar)
}


#=====================================================================
# The folowing programs estimate the MLE of TPL distribution.
#=====================================================================
#=====================================================================
# The "mu.tpl" program estimates the MLE of location parameter of tpl
#  distribution for given scale and shape(skewpar) parameters.
#=====================================================================
mu.tpl <- function(x, p){
    if (any(p <= 0 | p >= 1, na.rm = TRUE))
      stop("p must be between zero and one")
## Function for the number of digits
getndp <- function(x, tol=2*.Machine$double.eps)
        {
  ndp <- 0
  while(!isTRUE(all.equal(x, round(x, ndp), tol=tol))) ndp <- ndp+1
  if(ndp > -log10(tol)) warning("Tolerance reached, ndp possibly
underestimated.")
  ndp
      }
  x = sort(x)
  mu.all <- seq(min(x), max(x), length = 10^getndp(p))
  mu = mu.all[ceiling(length(mu.all) * p)]
  mu
}




#=====================================================================
# The "scale.tpl" program estimates the MLE of scale parameter of tpl
#  distribution for given location and shape parameters.
#=====================================================================
scale.tpl <- function(x, mu, p){

  if (any(p <= 0 | p >= 1, na.rm = TRUE))
      stop("p must be between zero and one")
  I1 <- function(x) ifelse(x <= 0, 1, 0)
  I2 <- function(x) ifelse(x > 0, 1, 0)
  n = length(x)
  zz = x - mu
  scale <- sum((-zz * I1(zz) / p) + (zz * I2(zz) / (1-p)))
  scale <- scale / (2 * n)
  scale
}

#====================================================================
# The "p.tpl" program estimates MLE of skewpar parameter of TPL
#  distribution for given location and scale parameters.
#====================================================================
p.tpl <- function(x, mu){
  I1 <- function(x) ifelse(x <= 0, 1, 0)
  I2 <- function(x) ifelse(x > 0, 1, 0)
  n = length(x)
  zz = x - mu
  temp1 <- sum(zz * I1(zz))
  temp2 <- sum(zz * I2(zz))
  p <- (temp1 +sqrt(-temp1 * temp2)) / (temp1 + temp2)
  p
}

#==================================================================
##########             log-Likelihood
#==================================================================
ell.tpl = function(location, scale, skewpar, x=dat)
	-sum(dtpl(x, location, scale, skewpar, log.arg = TRUE))

ell2.tpl = function(theta, x)
	ell.tpl(theta[1], theta[2], theta[3], x)
###################################################

#==================================================================
# The "mle.algp.tpl" estimates the exact MLE of three parameters
# of TPL distribution.
#==================================================================
mle.algp.tpl <- function(x, k = 1000, e.rate = 1e-10,
                         el.rate = length(x), max.iter=1000)
{
  n  = length(x)
  m  = max.iter
  pp = seq(.01, .99, length=k)
  M  = matrix (0.0, nr = k, nc = 4)
  for (j in 2 : length(pp) - 1){
    p   = pp[j]
    scale = sd(x)
    mu  = mean(x)

    t0 = c(mu = mu, scale = scale, p = p)
    t1 = t0;
    i = 1; e = 1; el = 1;
    while(i <= m && (e >= e.rate || el >= el.rate))
    {
    t1["p"]   = p
    t1["mu"]  = quantile(x, t1["p"])
    t1["scale"] = scale.tpl(x, t1["mu"], t1["p"])
    l0 = ell2.tpl(t0, x); l1 = ell2.tpl(t1, x)
    el = abs(l0 - l1); e = sum(abs(t1 - t0))
    t0 = t1;
    lt1 = ell2.tpl(t1, x)
    i = i + 1
    }
    M[j,] = c(t1, ell2.tpl(t1, x))
  }
  myname <- c("location", "scale", "skewpar", "-Log-liklihood")
  colnames(M) <- myname
  wh = which.min(M[-1000, 4])
  Ans = M[wh, ]
  Ans
}

#==================================================================
# The "mle.algp.tpl" estimates the exact MLE of three parameters
# of TPL distribution.
#==================================================================
mle.alg.tpl <- function(x, k = 1000, e.rate = 1e-10,
                        el.rate = length(x), max.iter=1000)
{
  n  = length(x)
  m  = max.iter
  pp = seq(.01, .99, length=k)
  M  = matrix (0.0, nr = k, nc = 4)
  for (j in 2 : length(pp) - 1){
    ppp   = pp[j]
    scale = sd(x)
    mu  = mean(x)

    t0 = c(mu = mu, scale = scale, p = ppp)
    t1 = t0;
    i = 1; e = 1; el = 1;
    while(i <= m && (e >= e.rate || el >= el.rate))
    {
    t1["p"]   = p.tpl(x, t1["mu"])
    t1["mu"]  = quantile(x, t1["p"])
    t1["scale"] = scale.tpl(x, t1["mu"], t1["p"])
    l0 = ell2.tpl(t0, x); l1 = ell2.tpl(t1, x)
    el = abs(l0 - l1); e = sum(abs(t1 - t0))
    t0 = t1;
    lt1 = ell2.tpl(t1, x)
    i = i + 1
    }
    M[j,] = c(t1, ell2.tpl(t1, x))
  }
  myname <- c("location", "scale", "skewpar", "-Log-liklihood")
  colnames(M) <- myname
  wh = which.min(M[-1000, 4])
  Ans = M[wh, ]
  Ans
}

 #==================================================================
# The "mle.algn.tpl" estimates the exact MLE of three parameters
# of TPL distribution according to paper.
#==================================================================
mle.algn.tpl <- function(x, e.rate = 1e-10, k = k,
                         el.rate = length(x), max.iter=1000)
{
  n  = length(x)
  m  = max.iter
  #pp = numeric(0)
  x = sort(x)
  M  = matrix (NA, nr = m, nc = 4)
  for (j in 2 : (length(x) - 1)){
    ppp   = j / n #pp[j]
    scale = sd(x)
    mu  = x[j]

    t0 = c(mu = mu, scale = scale, p = ppp)
    t1 = t0;
    i = 1; e = 1; el = 1;
    while(i <= m && (e >= e.rate || el >= el.rate))
    {
    t1["p"]   = p.tpl(x, t1["mu"])
    t1["mu"]  = quantile(x, t1["p"])
    t1["scale"] = scale.tpl(x, t1["mu"], t1["p"])
    l0 = ell2.tpl(t0, x); l1 = ell2.tpl(t1, x)
    el = abs(l0 - l1); e = sum(abs(t1 - t0))
    t0 = t1;
    lt1 = ell2.tpl(t1, x)
    i = i + 1
    }
    M[j,] = c(t1, ell2.tpl(t1, x))
  }
  myname <- c("location", "scale", "skewpar", "-Log-liklihood")
  colnames(M) <- myname
  wh = which.min(M[-1000, 4])
  Ans = M[wh, ]
  Ans
}

#==================================================================
# The "fit.tpl" fits tpl distribution on data.
#==================================================================

fit.tpl <- function(x, k = 1000, e.rate = 1e-10,
                    el.rate = length(x), max.iter=1000, exact.MLE = TRUE,
                    graph = TRUE,  main = NULL, xlab = NULL, col = NULL)
{
###################  check inputs ##################
  dat <- x
if (ncol(cbind(dat)) != 1)
    stop("x must be a vector or a one-column matrix")
if (is.null(xlab))
    xlab <- deparse(substitute(x))

if (is.null(col))
	col <- "light blue"

if (exact.MLE)
	if (k < 1000)
	stop("k musut be more than 1000 for exact MLE method")

if (graph != TRUE && any(mode(main) != "NULL" ||
    mode(col) != "NULL" || mode(xlab) != "NULL"))
warning("Some arguments use only for graph density")

  #### Initial Values #################
  n       <- length(dat)
  sigma   <- sd(dat)
  mu      <- mean(dat)
  pdat    <- mean(dat <= mu)
  min.dat <- min(dat)
  max.dat <- max(dat)

  den.dat <- density(dat)
  x       <- den.dat$x

 ######## Parameter estimation ############################
  theta <- matrix(0.0, ncol = 4, nrow = 1)

  if (exact.MLE) {fit.alg <- mle.algn.tpl(dat, k = k, e.rate = e.rate,
						el.rate = el.rate, max.iter = 1000)
                        theta[1, c(1:3)] = fit.alg[1:3] } else {
 NULL # Should work on it
	}

  theta[1, 4] <- ell2.tpl(theta[, 1:3], dat)
  colnames(theta) <- c("location", "scale", "skewpar", "-Log-likelihood")
  if (theta[1, 3] < 0.1 | theta[1, 3] > 0.9)
  	warning("The tpl distribution is not apropriate for this data set")

############### Making Y : The density function in points X = den.dat$x

  y <- matrix(ncol = length(x), nrow = 1)

  y[1, ] <- dtpl(x, theta[1, 1], theta[1, 2], theta[1, 3])

####################### ABE and SQE distance ############
  ab.error <- sum(abs(y[1, ] - den.dat$y))
  sq.error <- sum((y[1, ] - den.dat$y)^2)

################   KS-TEST  ###########################
  Fyi <- function(x)
  ptpl(x, theta[1,1], theta[1,2], theta[1,3])
  mykstest <- (ks.test(dat, "ptpl", theta[1,1], theta[1,2], theta[1,3] ))

  txt <- paste("Kolmogorov-Smirnov test", "\n",
             "D = ", round(mykstest$statistic,4), "\n",
             "P-value = ", round(mykstest$p.value,3), sep="")


######################### Graph ######################
  if (graph) {
    myxlab <- xlab
    max.y  <- c(max(y[1, ], den.dat$y))
    My     <- max(max.y[!(is.na(max.y))])
    par(mar = c(4.2, 4, 2, 1))
    hist(dat, freq = FALSE, ylim = c(0, My), col = col,
  	    main = main, xlab = xlab, las = 1)

    lines(x, y[1, ], type = "l", col = 3, lty = 1, lwd = 2)

    lines(den.dat, lty = 4)

    text(ifelse(median(dat) < mu, quantile(dat, 0.95), quantile(dat, 0.05)),
                quantile(den.dat$y, 0.75), txt, adj = c(0,1), pos = 3)

    legend(ifelse(median(dat) < mu, "topright", "topleft"),
          c("Two-piece Laplace", "Non-parametric"), col = c("green", 1),
          text.col = 1, lty = c(1, 4), pch = c(-1, -1), lwd = c(2, 1),
          bg = 'gray90',   merge = TRUE)

    rug(dat)
  }
##################### output ########################
 # if (exact.MLE) {
 # 	fit <- fit.alg} else {
 #	fit <- fit.vglm
  #}
  list( fit = theta, KSTest = mykstest )

}
