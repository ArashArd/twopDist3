#Last modiffied by Arash 26th May /2020
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  This is the derivative of location parameter of TPLogistic dist.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XUtplog = function(y, mylocat, myscale, myskew){
zedd  <- (y - mylocat) / myscale # Is a n x NOS matrix
cond2 <- (zedd > 0)

dl.dlocat        <- 2 * tanh(zedd / (4 * myskew)) / (myskew) # cond1
dl.dlocat[cond2] <- (2 * tanh(zedd / (4 * (1-myskew))) / ( 1-myskew))[cond2]
dl.dlocat        <- dl.dlocat / (4 * myscale)

dl.dlocat
}



Scale2tplogis = function(EE, skewpar)
{
nn = length(EE)
cond1 <- (EE <= 0)
cc = 0.324
sighat = cc * (sum(abs(EE[cond1]/skewpar)) +
                                  sum(abs(EE[!cond1]/(1-skewpar))))/(nn * 0.9)
sighat}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  TPL.logisticEst
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 logisticREst = function(y, x, myscale, myskew, Tol = 1e-5){
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
sighat0 = Scale2tplogis(EE0, myskew)
myscale = sighat0
while(any(abs(BWo - BWn) > Tol)){
      BWo = BWn
      muo = mun
      XU  = XUtplog(y, muo, myscale, myskew)
      XtU = t(X) %*% XU
      EE = y - mun
      WW =  matrix(rep(as.vector(XU / (EE)), ncol(X)), nc = ncol(X)) * X
      temp15 = (t(y) %*% (XU / (EE)))
      BWn = solve(t(WW) %*% (X)) %*% (t(WW) %*% y)
      sighat = Scale2tplogis(EE, myskew)
      Temp20 = myskew * (1 - myskew)
      Temp21 = 4 * ((1 - 3 * Temp20) * ((pi^2 / 3) - log(2)^2) + Temp20 * log(2)^2)
      v = Temp21 * sighat^2

      VV  = v * solve(t(X) %*% (X)) * (0.084/ Temp20)^3
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
        sigma = sighat,  df = df, p = myskew)
        }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  The class of TPLogR
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 TPLogR <- function(x, ...) UseMethod("TPLogR")

#=========================================================================
#       TPLogR.default
#=========================================================================
TPLogR  <- function(formula, data, weights,  method = "qr",
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE,  contrasts = NULL,
    offset, myskew, ...) {

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
         est <- logisticREst(y, x, myscale, myskew, Tol = 1e-5)
        }

  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$p <- myskew
  #est$call <- match.call()
  class(est) <- "TPLogR"
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

print.TPLogR <- function(x, ...) {
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
}

#=========================================================================
#         Arash has to change the tests. The tests are not correct
#=========================================================================
summary.TPLogR <- function(object, ...) {
   myskew = object$p
  temp30 = myskew * (1 - myskew)

  varepsilon = (1- 3*temp30) * (pi^2/3 - log(2)^2) + temp30 * log(2)^2
    se <- sqrt(diag(object$vcov)* varepsilon)
  tval <- coef(object) / se
  TAB <- cbind(
    Estimate = coef(object), StdErr = se, t.value = tval,
    #if (object$p == 0.5) p.value = 2 * pt(-abs(tval), df = object$df) else
	p.value = 2 * pnorm(-abs(tval))
         )
  res <- list(call = object$call,  coefficients = TAB)
  class(res) <- "summary.TPLogR"
  res
}

#=========================================================================
#        print.summary
#=========================================================================

print.summary.TPLogR <- function(x, ...) {
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
#         TPLogR.formula
#=========================================================================
TPLogR.formula <- function(formula, data = list(), ...) {
  mf <- model.frame(formula = formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- TPLogR(x, y, myskew, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#=========================================================================
#        predict.TPLogR
#=========================================================================
predict.TPLogR <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) y <- fitted(object) else{
   if(!is.null(object$formula)){
  ## model has been fitted using formula interface
  x <- model.matrix(object$formula, newdata)
   }      else{x <- newdata}
  y <- as.vector(as.matrix(cbind(1,x)) %*% as.matrix(coef(object), nc = 1))
   }
  y
}

