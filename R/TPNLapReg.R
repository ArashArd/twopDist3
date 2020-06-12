
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  From here there are some codes for quantile regression normal Laplace model
#  The codes have written by S3 class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  This is the derivative of location parameter of two-piece normal Laplace dist.
# Warning : The fisrst derivative needs to review for discontinuity
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


XUtptpnl = function(y, mylocat, myscale, myskew){

zedd  <- (y - mylocat) / myscale # Is a n x NOS matrix
cond2 <- (zedd > 0)

dl.dlocat        <- zedd / (8 * myskew^2) # cond1
dl.dlocat[cond2] <- (1 / (sqrt(2 * pi) * (1 - myskew)))
#dl.dlocat        <- dl.dlocat / (4 * myscale)

dl.dlocat
}


###################################################
I.location=function(location,x=dat)
c(sum((x-location)^2*(x<=location)),sum((x-location)*(x>location)))
###############################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Scale function for finding Sigma
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Scaletpnl = function(location,skewpar,x=dat)
{
n=length(x)
i=I.location(location,x)
i1=i[1];i2=i[2]
a=(4*n*skewpar^2)*(n*(1-skewpar)*sqrt(2*pi))
b=-i2*(4*n*skewpar^2)
c=-i1*(n*(1-skewpar)*sqrt(2*pi))
D=b^2-4*a*c
if(D<0)
{cat('Warning: Delta for scale is negative\n')
return(.0001)
}
if(a==0) return(0.0001)
return((-b + sqrt(D))/(2 * a))
}


Scale2tpnl = function(E,skewpar)
{
n=length(E)
A = sum((E[E<0])^2) / (8 * skewpar^2)
B = sum((E[E > 0]) ) / ((1-skewpar) * sqrt(2*pi))
sig = (B + sqrt(B^2 + 8 * n * A)) / (2*n)
sig}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  TPNL_REst
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TPNL_REst = function(y, x, myscale, myskew, Tol = 1e-5, maxit = 100){
qx0 <- qr(x)
betao = solve.qr(qx0, y)
Bo = betao
X = cbind(1, x)
muo = x %*% betao
mun = x %*% (betao +.1)
BWo = (betao)
BWn = (betao +.1)
X = x
EE0 = y - muo
sighat0 = Scale2tpnl(EE0, myskew)
myscale = sighat0
mymaxit = 0
while(any(abs(BWo - BWn) > Tol)){
mymaxit = mymaxit + 1
if (mymaxit == maxit) break
      BWo = BWn
      muo = mun
      XU  = XUtptpnl(y, muo, myscale, myskew)
      XtU = t(X) %*% XU
      EE = y - mun
      WW =  matrix(rep(as.vector(XU / (EE)), ncol(X)), nc = ncol(X)) * X
      temp15 = (t(y) %*% (XU / (EE)))
      #BWn = solve(t(WW) %*% (X)) %*% (t(WW) %*% y)
      cma <- chol(t(WW) %*% (X))

      BWn = chol2inv(cma) %*% (t(WW) %*% y)

      dat = EE
      #sighat = Scale(quantile(dat, myskew), myskew, scale(y)) #scale(y)
      sighat = Scale2tpnl(EE, myskew)
      Temp20 = myskew * (1 - myskew)
      mygamma = (4 * myskew + pi * (1 - myskew))/( 8 * Temp20 * pi )
      Imu = (1 + myskew * (2- pi) + pi)/(4 * Temp20 * pi)

      v = Imu/mygamma^2  * sighat^2

      VV  = v * solve(t(X) %*% (X))
      ## Must be modify by Arash
      myscale =  sighat

      mun = X %*% BWn
}

dn <- colnames(x)
if (is.null(dn))
        dn = paste0("x", 1L:p)
BetaHat = as.vector(t(BWn))
names(BetaHat) = dn
df <- nrow(X) - ncol(X)
 list(coefficients = BetaHat , vcov = VV,
        sigma = sighat,  df = df, p = myskew, MaxIt =  mymaxit)
        }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  The class of  two-piece normal Laplace regression (TPNLR)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 TPNLR <- function(x, ...) UseMethod("TPNLR")

#=========================================================================
#       TPLogR.default
#=========================================================================
 TPNLR <- function(formula, data, weights,  method = "qr",
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE,  contrasts = NULL,
    offset, myskew, myscale, ...) {

  x <- as.matrix(x)
  y <- as.numeric(y)
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m  <- match(c("formula", "data", #"subset", "weights", "na.action",
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    if (method == "model.frame") return(mf)

    mt <- attr(mf, "terms")
    y  <- model.response(mf, "numeric")
        offset <- as.vector(model.offset(mf))
        if (is.empty.model(mt)) {
        x <- NULL
        est <- list(coefficients = if(is.matrix(y)) matrix(, 0, 3) else numeric(),
                           residuals = y, fitted.values = 0 * y, rank = 0L )
        if (!is.null(offset)) {
            est$fitted.values <- offset
            est$residuals <- y - offset
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)
         est <- TPNL_REst(y, x, myscale, myskew, Tol = 1e-5)
        }

  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$p <- myskew
  #est$call <- match.call()
  class(est) <- "TPNLR"
    est$offset <- offset
    est$contrasts <- attr(x, "contrasts")
    est$xlevels <- .getXlevels(mt, mf)
    est$call <- cl
    est$terms <- mt
    if (model)
        est$model <- mf
    if (ret.x)
        est$x <- x
    if (ret.y)
        est$y <- y
    if (!qr)
        est$qr <- NULL

  est
}

#=========================================================================
#         print.TPLogR the print for the regression
#=========================================================================

print.TPNLR <- function(x, ...) {
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
}

#=========================================================================
#         Arash has to change the tests. The tests are not correct
#=========================================================================
summary.TPNLR <- function(object, ...) {
   myskew = object$p  ####### Checkkkkkkkkkkkkkkkk
  temp30 = myskew * (1 - myskew)

  #varepsilon = (1- 3*temp30) * (pi^2/3 - log(2)^2) + temp30 * log(2)^2
    se <- sqrt(diag(object$vcov)) # * varepsilon
  tval <- coef(object) / se
  TAB <- cbind(
    Estimate = coef(object), StdErr = se, t.value = tval,
    #if (object$p == 0.5) p.value = 2 * pt(-abs(tval), df = object$df) else
	p.value = 2 * pnorm(-abs(tval))
         )
  res <- list(call = object$call,  coefficients = TAB)
  class(res) <- "summary.TPNLR"
  res
}

#=========================================================================
#        print.summary
#=========================================================================

print.summary.TPNLR <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
}


#=========================================================================
#       confint.TPNLR
#=========================================================================
confint.TPNLR <- function(object, parm, level = 0.95, ...)
{
    cf <- coef(object)
    ses <- sqrt(diag(object$vcov))
    pnames <- names(ses)
    if (is.matrix(cf))
        cf <- setNames(as.vector(cf), pnames)
    if (missing(parm))
        parm <- pnames
    else if (is.numeric(parm))
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qnorm(a)
    pct <- stats:::format.perc(a, 3)
    ci <- array(NA_real_, dim = c(length(parm), 2L),
                              dimnames = list(parm,  pct))
    ci[] <- cf[parm] + ses[parm] %o% fac
    ci
}


#=========================================================================
#         TPNLR.formula TPNLR.formula
#=========================================================================
formula.TPNLR <- function(formula, data = list(), ...) {
  mf <- model.frame(formula = formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- TPNLR(x, y, myskew, myscale, ...)  # .default
  est$call <- match.call()
  est$formula <- formula
  est
}


#=========================================================================
#        predict.TPNLR
#=========================================================================
predict.TPNLR <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) y <- fitted(object) else{
   if(!is.null(object$formula)){
  ## model has been fitted using formula interface
  x <- model.matrix(object$formula, newdata)
   }      else{x <- newdata}
  y <- as.vector(as.matrix(cbind(1,x)) %*% as.matrix(coef(object), nc = 1))
   }
  y
}






