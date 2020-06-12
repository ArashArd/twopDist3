
#Last modified 26/05/2020 the vcov needs to check and also
# there are two algorithms which i need to add that.
#=========================================================================
#        The parameter estimation for a regression model using TPN distribution
#         The codes are objected orinted and S3
#=========================================================================


#=========================================================================
#         sig.tpn
#=========================================================================
 sig.tpn <- function(y, mu, p) {
    if (any(p <= 0 | p >= 1, na.rm = TRUE))
      stop("p must be between zero and one")
  I1 <- function(y) ifelse(y <= 0, 1, 0)
  I2 <- function(y) ifelse(y > 0, 1, 0)
  n = length(y)
  zz = y - mu
  zzL =  zz * I1(zz) / p
  zzR =  zz * I2(zz) /(1-p)
  sig <- t(zzL) %*% zzL + t(zzR) %*% zzR
  sig <- sqrt(sig / (4 * n))
  as.vector(sig)
}

#=========================================================================
#         tpnrnEst
#=========================================================================
tpnrnEst =
  function(x, y, p) {
  #x <- cbind(1, x)
  qx0 <- qr(x)
  ## compute (x'x)^(-1) x'y
  coef0 <- solve.qr(qx0, y)
  W0 = y - x %*% coef0
   if (p == 0.5) {Beta.hat = coef0; W = W0; qx = qx0} else {
     p  = p
     W0 = ifelse(W0 > 0, 1 / (8 * (1 - p)^2), 1 / (8 * p^2) )
     W0 = diag(as.vector(W0))
     qx <- qr(W0^0.5 %*% x)
     yt <- W0^0.5 %*% y
     Beta.hat <- solve.qr(qx, as.vector(yt))
     Beta.hat0 = coef0
      while (any(abs(Beta.hat0 - Beta.hat) > 1e-6)) {
        Beta.hat0 <- Beta.hat
        W = y - x %*% Beta.hat0
        W = ifelse(W > 0, 1/(8 * (1 - p)^2), 1/(8 * p^2) )
        W = diag(as.vector(W))
        qx <- qr(W^0.5 %*% x) # W^.5
        yt <- W^0.5 %*% y      # W^.5

        Beta.hat <- solve.qr(qx, as.vector(yt))
        #cat(Beta.hat, "\n")
        }
     }
  #The following quantities must be correct.
  df <- nrow(x) - ncol(x)
  sigma <- sig.tpn(y, x %*% Beta.hat, p)
  ## compute sigma^2 * (x'x)^-1
  #qx <- qr(W^0.5 %*% x)
 #######################
  ## Need to check
 ######################
  vcov <- sigma^2 * chol2inv(qx$qr)
  colnames(vcov) <- rownames(vcov) <- colnames(x)

  list(coefficients = Beta.hat, vcov = vcov,
        sigma = sigma,  df = df, p = p)
    }


#=========================================================================
#         tpn
#=========================================================================
 tpnr1Est =     function(x, y, p) {

    mod0 = lsfit(x, y)
     Beta.hat0 <- mod0$coef
     W0 = mod0$residuals
     p  = p
     W0 = ifelse(W0 > 0, 1 / (8 * (1 - p)^2), 1 / (8 * p^2) )
     W0 = as.vector(W0)
     mod1 = lsfit(x, y, wt = W0)
     Beta.hat <- mod1$coef
     W0 = mod1$residuals
     W0 = ifelse(W0 > 0, 1 / (8 * (1 - p)^2), 1 / (8 * p^2) )
     W0 = as.vector(W0)
      while (any(abs(Beta.hat0 - Beta.hat) > 1e-6)) {
        Beta.hat0 <- Beta.hat
        modn = lsfit(foodexp, income, wt = W0)
        W0 = modn$residuals
        W0 = ifelse(W0 > 0, 1 / (8 * (1 - p)^2), 1 / (8 * p^2) )
        W0 = as.vector(W0)
        Beta.hat = modn$coef
        cat(Beta.hat, "\n")
        }

  #The following quantities must be correct.
   x = cbind(1,x)
  df <- nrow(x) - ncol(x)
  sigma <- sig.tpn(y, x %*% Beta.hat, p)
  ## compute sigma^2 * (x'x)^-1
   qx <- qr(W0^0.5 %*% x)
  vcov <- sigma^2 #* chol2inv(qx$qr)
  #colnames(vcov) <- rownames(vcov) <- colnames(x)

  list(coefficients = Beta.hat, vcov = vcov,
        sigma = sigma,  df = df, p = p)
    }



#=========================================================================
#         tpnrn
#=========================================================================
TPNR <- function(x, ...) UseMethod("TPNR")
#=========================================================================
#         tpnrn.default
#=========================================================================
TPNR <- function(formula, data, weights,  method = "qr",
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE,  contrasts = NULL,
    offset, p, ...) {

  x <- as.matrix(x)
  y <- as.numeric(y)


    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", #"subset", "weights", "na.action",
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame")
        return(mf)

    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    offset <- as.vector(model.offset(mf))

    if (is.empty.model(mt)) {
        x <- NULL
        est <- list(coefficients = if (is.matrix(y)) matrix(, 0,
            3) else numeric(), residuals = y, fitted.values = 0 * y,
             rank = 0L )
        if (!is.null(offset)) {
            est$fitted.values <- offset
            est$residuals <- y - offset
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)
         est <- tpnrnEst(x, y, p)
        }
  ##############################################

  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$p <- p
  #est$call <- match.call()
  class(est) <- "TPNR"

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
#         print.tpnrn
#=========================================================================

print.TPNR <- function(x, ...) {
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
}

#=========================================================================
#         summary.tpnrn
#=========================================================================
summary.TPNR <- function(object, ...) {
   myskew = object$p
  temp30 = myskew * (1 - myskew)
  ## Need to check
  varepsilon = 4 * temp30
    se <- sqrt(diag(object$vcov)* varepsilon)
  tval <- coef(object) / se
  TAB <- cbind(
    Estimate = coef(object), StdErr = se, t.value = tval,
    #if (object$p == 0.5) p.value = 2 * pt(-abs(tval), df = object$df) else
	p.value = 2 * pnorm(-abs(tval))
         )
  res <- list(call = object$call,  coefficients = TAB)
  class(res) <- "summary.TPNR"
  res
}

#=========================================================================
#         print.summary.tpnrn
#=========================================================================

print.summary.TPNR <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE)
}

#=========================================================================
#       confint.TPNR
#=========================================================================
confint.TPNR <- function(object, parm, level = 0.95, ...)
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
#         tpnrn.formula
#=========================================================================
formula.TPNR <- function(formula, data = list(), ...) {
  mf <- model.frame(formula = formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- TPNR(x, y, p, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}

#=========================================================================
#         predict.tpnrn
#=========================================================================
predict.TPNR <- function(object, newdata = NULL, ...) {
  if(is.null(newdata)) y <- fitted(object) else{
   if(!is.null(object$formula)){
  ## model has been fitted using formula interface
  x <- model.matrix(object$formula, newdata)
   }      else{x <- newdata}
  y <- as.vector(x %*% coef(object))
   }
  y
}



