# agenstudentt.R

# Last modified:
# 20110429 by Thomas Yee 
# 20110502 by Thomas Yee 
# 20120519 by Arash

# Reference: zhu and galbraith, 2010, journal of econometrics, 157: 297--305.
# Also, Appendix C that goes with it.


############################# DDS.twopt #############################
DDS.twopt <- function (df) {
# Reference: p.301, eqn (29).
# Handles df = Inf.

  ans = digamma(0.5 * (df + 1)) - digamma(0.5 * df)
  ans[is.infinite(df)] = 0
  ans[df <= 0] = NA
  ans
}


###############################  DDSp.twopt #########################
DDSp.twopt <- function (df) {
# Reference: p.301, eqn (29).
# Handles df = Inf.

  ans = 0.5 * (trigamma(0.5 * (df + 1)) - trigamma(0.5 * df))
# ans[is.infinite(df)] = 0   # Not needed
  ans[df <= 0] = NA
  ans
}


############################# Kkn.twopt #############################
Kkn.twopt <- function(df, log.arg = FALSE) {
# Reference: p.298, between eqn (1) and eqn (2).
# The log-scale helps avoid overflow:
# Handles df = Inf.

  finite.df <- is.finite(df)
  df.finite <- df[finite.df]
  log.kkn <- df
  log.kkn[finite.df] <- lgamma((df.finite + 1) / 2) -
                        lgamma(df.finite / 2) - 0.5 * log(pi * df.finite)
# log.kkn <- lgamma((df + 1) / 2) - lgamma(df / 2) - 0.5 * log(pi * df)
  log.kkn[!finite.df] <- -0.5 * log(2 * pi)
  log.kkn[df <= 0] = NA
  if (log.arg) {
    log.kkn
  } else {
    exp(log.kkn)
  }
}

################      dtwopt       ##################################

dtwopt <- function(x, location = 0, scale = 1, skewpar = 0.5,
                 dfL = 2, dfR = 2, log.arg = FALSE) {

# Reference: equation (5)

  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 |
          dfL     <= 0 |
          dfR     <= 0, na.rm = TRUE))
    stop("some parameters out of bound")

# Recycle the vectors to equal lengths
  LLL = max(length(x), length(location), length(scale),
            length(skewpar), length(dfL), length(dfR))
  if (length(x)        != LLL)   x        = rep(x,        length = LLL)
  if (length(location) != LLL)   location = rep(location, length = LLL)
  if (length(scale)    != LLL)   scale    = rep(scale,    length = LLL)
  if (length(skewpar)  != LLL)   skewpar  = rep(skewpar,  length = LLL)
  if (length(dfL)      != LLL)   dfL      = rep(dfL,      length = LLL)
  if (length(dfR)      != LLL)   dfR      = rep(dfR,      length = LLL)


  zedd <- (x - location) / scale

  temp2 <- skewpar * Kkn.twopt(dfL)
  skewpar.s <- temp2 / (temp2 + (1 - skewpar) * Kkn.twopt(dfR)) # aka alphas

  log.s1 <- log(skewpar / skewpar.s) +
            Kkn.twopt(dfL, log.arg = TRUE) +
            (-0.5*(dfL + 1)) * log1p((zedd / (2*     skewpar.s ))^2 / dfL) +
            -log(scale)
  log.s2 <- log1p(-skewpar) - log1p(-skewpar.s) +
            Kkn.twopt(dfR, log.arg = TRUE) +
            (-0.5*(dfR + 1)) * log1p((zedd / (2*(1 - skewpar.s)))^2 / dfR) +
            -log(scale)
  logdensity <- log.s1
  logdensity[zedd > 0] <- log.s2[zedd > 0]


  indexL = is.infinite(dfL) & zedd <= 0
  indexR = is.infinite(dfR) & zedd >= 0
  if (any(indexL)) {
    logdensity[indexL] = dnorm(zedd[indexL] / (2 * skewpar.s[indexL]),
                               log = TRUE) - log(scale[indexL]) +
                         log(skewpar[indexL] / skewpar.s[indexL])
  }
  if (any(indexR)) {
    logdensity[indexR] = dnorm(zedd[indexR] / (2 * (1 - skewpar.s[indexR])),
                               log = TRUE) - log(scale[indexR]) +
                         log1p(-skewpar[indexR]) -
                         log1p(-skewpar.s[indexR])
  }

  if (log.arg) logdensity else exp(logdensity)
}

################      ptwopt         ################################
ptwopt <- function(q, location = 0, scale = 1, skewpar = 0.5,
                 dfL = 2, dfR = 2) {

  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 |
          dfL     <= 0 |
          dfR     <= 0,
          na.rm = TRUE))
    stop("some parameters out of bound")

# Reference: equation (6)
  # Recycle the vectors to equal lengths
  # Recycle the vectors to equal lengths
  LLL = max(length(q), length(location), length(scale),
            length(skewpar), length(dfL), length(dfR))
  if (length(q)        != LLL)    q        = rep(q,        length = LLL)
  if (length(location) != LLL)   location = rep(location, length = LLL)
  if (length(scale)    != LLL)   scale    = rep(scale,    length = LLL)
  if (length(skewpar)  != LLL)   skewpar  = rep(skewpar,  length = LLL)
  if (length(dfL)      != LLL)   dfL      = rep(dfL,      length = LLL)
  if (length(dfR)      != LLL)   dfR      = rep(dfR,      length = LLL)

  temp2 <- skewpar * Kkn.twopt(dfL)
  skewpar.s <- temp2 / (temp2 + (1 - skewpar) * Kkn.twopt(dfR))

  zedd <- (q - location) / scale
  prt <- 2 * skewpar * pt(pmin(zedd, 0) / (2 * skewpar.s), df = dfL) +
         2 * (1 - skewpar) *
        (pt(pmax(zedd, 0) / (2 * (1 - skewpar.s)), df = dfR) - 0.5)
  prt
}



#####################  ############################################
qtwopt <- function(p, location = 0, scale = 1, skewpar = 0.5,
                 dfL = 2, dfR = 2) {

  pp = p
  if (any(pp      <= 0 ||
          pp      >= 1 |
          skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 |
          dfL     <= 0 |
          dfR     <= 0, na.rm = TRUE))
    stop("some parameters out of bound")
    
    
  # Recycle the vectors to equal lengths
  LLL = max(length(pp), length(location), length(scale),
            length(skewpar), length(dfL), length(dfR))
  if (length(pp)       != LLL)   pp       = rep(pp,       length = LLL)
  if (length(location) != LLL)   location = rep(location, length = LLL)
  if (length(scale)    != LLL)   scale    = rep(scale,    length = LLL)
  if (length(skewpar)  != LLL)   skewpar  = rep(skewpar,  length = LLL)
  if (length(dfL)      != LLL)   dfL      = rep(dfL,      length = LLL)
  if (length(dfR)      != LLL)   dfR      = rep(dfR,      length = LLL)

# Reference: equation (7)
# It needs to check
  temp2 <- skewpar * Kkn.twopt(dfL)
  skewpar.s <- temp2 / (temp2 + (1 - skewpar) * Kkn.twopt(dfR))

  qrt <- 2 * skewpar.s * qt(pmin(pp, skewpar) / (2 * skewpar), df = dfL) +
         2 * (1 - skewpar.s) * qt((pmax(pp, skewpar) +
             1 - (2 * skewpar)) / (2 * (1 - skewpar)), df = dfR)
  qrt <-  qrt * scale + location
  qrt
}

###########     rtwopt     ##########################################

rtwopt <- function(n, location = 0, scale = 1, skewpar = 0.5,
                 dfL = 2, dfR = 2) {

  qtwopt(p = runif(n), location = location, scale = scale,
       skewpar = skewpar, dfL = dfL, dfR = dfR)
}


#######################  agstf family  ############################

twopt <- function(llocation = "identity",   elocation = list(),
                  lscale  = "loge",         escale    = list(),
                  ldfL    = "loglog",       edfL      = list(),
                  ldfR    = "loglog",       edfR      = list(),
                  lskew   = "logit",        eskew     = list(),
                  ilocation = NULL, iscale = NULL, idfL = NULL,
                  idfR      = NULL, iskew  = NULL, 
                  smallno = 1.0e-10,
                  imethod = 1,  zero = 2:5)
{
# Notes:
# 1. "sqw" is a 3-character acronym for "skew".
# 2. Yet to do: get ncol(y) > 1 working.
#


  lloc <- llocation; lsca <- lscale;  ldfL <- ldfL; ldfR = ldfR; lsqw = lskew
  eloc <- elocation; esca <- escale;  edfL <- edfL; edfR = edfR; esqw = eskew
  iloc <- ilocation; isca <- iscale;  idfL <- idfL; idfR = idfR; isqw = iskew

  if (mode(lloc) != "character" && mode(lloc) != "name")
    lloc <- as.character(substitute(lloc))
  if (!is.list(eloc)) eloc <- list()

  if (mode(lsca) != "character" && mode(lsca) != "name")
    lsca <- as.character(substitute(lsca))
  if (!is.list(esca)) esca <- list()

  if (mode(ldfL) != "character" && mode(ldfL) != "name")
    ldfL <- as.character(substitute(ldfL))
  if (!is.list(edfL)) edfL <- list()

  if (mode(ldfR) != "character" && mode(ldfR) != "name")
    ldfR <- as.character(substitute(ldfR))
  if (!is.list(edfR)) edfR <- list()

  if (mode(lsqw) != "character" && mode(lsqw) != "name")
    lsqw <- as.character(substitute(lsqw))
  if (!is.list(esqw)) esqw <- list()

  if (length(iloc))
    if (!is.Numeric(iloc))
      stop("bad input in argument 'ilocation'")
  if (length(isca))
    if (!is.Numeric(isca, positive = TRUE))
      stop("argument 'iscale' should be > 0")

  if (length(idfL))
    if (!is.Numeric(idfL) || any(idfL <= 1))
      stop("argument 'idfL' should be > 1")

  if (length(idfR))
    if (!is.Numeric(idfR) || any(idfR <= 1))
      stop("argument 'idfR' should be > 1")

  if (length(isqw))
    if (!is.Numeric(isqw, positive = TRUE) || any(isqw >= 1))
      stop("argument 'iskew' should be < 1 and > 0")

  new("vglmff",
  blurb = c("Asymmetric generalized Student t-distribution\n\n",
            "Link:     ",
            namesof("location", lloc, earg = eloc), ", ",
            namesof("scale",    lsca, earg = esca), ", ",
            namesof("dfL",      ldfL, earg = edfL), ", ",
            namesof("dfR",      ldfR, earg = edfR), ", ",
            namesof("skewness", lsqw, earg = esqw)),
  constraints = eval(substitute(expression({
    constraints <- cm.zero.vgam(constraints, x, .zero, M)
  }), list( .zero = zero ))),
  infos = eval(substitute(function(...) {
    list(Musual = 5,
         zero   = .zero)
  }, list( .zero = zero ))),

  initialize = eval(substitute(expression({
    Musual <- 5

    y <- as.matrix(y)
    extra$NOS <- NOS <- ncoly <- ncol(y)  # Number of species
    extra$Musual <- Musual
    M <- Musual * ncoly #

    if (ncoly > 1) 
      stop("cannot handle matrix response yet")

    mynames1 <- ptwopte("location", if (NOS > 1) 1:NOS else "", sep = "")
    mynames2 <- ptwopte("scale",    if (NOS > 1) 1:NOS else "", sep = "")
    mynames3 <- ptwopte("dfL",      if (NOS > 1) 1:NOS else "", sep = "")
    mynames4 <- ptwopte("dfR",      if (NOS > 1) 1:NOS else "", sep = "")
    mynames5 <- ptwopte("skewness", if (NOS > 1) 1:NOS else "", sep = "")
    predictors.names <-
        c(namesof(mynames1, .lloc, earg = .eloc, tag = FALSE),
          namesof(mynames2, .lsca, earg = .esca, tag = FALSE),
          namesof(mynames3, .ldfL, earg = .edfL, tag = FALSE),
          namesof(mynames4, .ldfR, earg = .edfR, tag = FALSE),
          namesof(mynames5, .lsqw, earg = .esqw, tag = FALSE))
    predictors.names <-
      predictors.names[interleave.VGAM(Musual * NOS, M = Musual)]

    if (!length(etastart)) {

# zz needs improvement overall

      init.loc <- if (length( .iloc )) {
        .iloc 
      } else if ( .imethod == 2) {
        apply(y, 2, median)
      } else if ( .imethod == 3) {
# Transpose it since byrow = TRUE below
         c(t(y))
      } else {
           colSums(w * y) / sum(w)
      }

      sdvec <- apply(y, 2, sd)
      init.sca <- if (length( .isca )) .isca else
                  sdvec / 2.3

      init.dfL <- if (length( .idfL )) .idfL else   
                   ((2 * sdvec^2 / init.sca^2) / (sdvec^2 / init.sca^2  - 1))
      if (!is.Numeric(init.dfL) || init.dfL <= 1)
          init.dfL <- rep(5.0, len = ncoly)
  
      init.dfR <- if (length( .idfR )) .idfR else   
                   ((2 * sdvec^2 / init.sca^2) / (sdvec^2 / init.sca^2  - 1))
      if (!is.Numeric(init.dfR) || init.dfR <= 1)
          init.dfR <- rep(5.0, len = ncoly)
  
      init.sqw <- if (length( .isqw )) .isqw else
                    mean(y < mean(y))
      if (!is.Numeric(init.sqw, posit = TRUE) || init.sqw >= 1)
          init.sqw <- 0.55
  
 print("summary(init.loc)")
 print( summary(init.loc) )
 print("summary(init.sca)")
 print( summary(init.sca) )
 print("summary(init.dfL)")
 print( summary(init.dfL) )
 print("summary(init.dfR)")
 print( summary(init.dfR) )
 print("summary(init.sqw)")
 print( summary(init.sqw) )
  
      mat1 <- matrix(theta2eta(init.loc, .lloc, earg = .eloc), n, NOS,
                     byrow = TRUE)
      mat2 <- matrix(theta2eta(init.sca, .lsca, earg = .esca), n, NOS,
                     byrow = TRUE)
      mat3 <- matrix(theta2eta(init.dfL, .ldfL, earg = .edfL), n, NOS,
                     byrow = TRUE)
      mat4 <- matrix(theta2eta(init.dfR, .ldfR, earg = .edfR), n, NOS,
                     byrow = TRUE)
      mat5 <- matrix(theta2eta(init.sqw, .lsqw, earg = .esqw), n, NOS,
                     byrow = TRUE)
      etastart <- cbind(mat1, mat2, mat3, mat4, mat5)
  
 print("head(etastart)")
 print( head(etastart) )
  
      etastart <- etastart[, interleave.VGAM(ncol(etastart), M = Musual)]
    }
  }), list( .lloc = lloc, .eloc = eloc, .iloc = iloc,
            .lsca = lsca, .esca = esca, .isca = isca,
            .ldfL = ldfL, .edfL = edfL, .idfL = idfL,
            .ldfR = ldfR, .edfR = edfR, .idfR = idfR,
            .lsqw = lsqw, .esqw = esqw, .isqw = isqw,
            .imethod = imethod ))),
  inverse = eval(substitute(function(eta, extra = NULL) {
    NOS <- extra$NOS
    Musual <- extra$Musual
    Loc <-  eta2theta(eta[, Musual*(1:NOS) - 4], .lloc, earg = .eloc)
    DfL <-  eta2theta(eta[, Musual*(1:NOS) - 2], .ldfL, earg = .edfL)
    DfR <-  eta2theta(eta[, Musual*(1:NOS) - 1], .ldfR, earg = .edfR)
    Loc[DfL < 1] <- NA
    Loc[DfR < 1] <- NA
    Loc
  }, list(
    .lloc = lloc, .lsca = lsca, .ldfL = ldfL, .ldfR = ldfR, .lsqw = lsqw,
    .eloc = eloc, .esca = esca, .edfL = edfL, .edfR = edfR, .esqw = esqw
     ))),
  last = eval(substitute(expression({
    misc$link <- c(rep( .lloc, length = NOS),
                   rep( .lsca, length = NOS),
                   rep( .ldfL, length = NOS),
                   rep( .ldfR, length = NOS),
                   rep( .lsqw, length = NOS))
    misc$link <- misc$link[interleave.VGAM(Musual * NOS, M = Musual)]
    temp.names <- c(mynames1, mynames2, mynames3, mynames4, mynames5)
    temp.names <- temp.names[interleave.VGAM(Musual * NOS, M = Musual)]
    names(misc$link) <- temp.names

    misc$earg <- vector("list", Musual * NOS)
    names(misc$earg) <- temp.names
    for(ii in 1:NOS) {
        misc$earg[[Musual*ii - 4]] <- .eloc
        misc$earg[[Musual*ii - 3]] <- .esca
        misc$earg[[Musual*ii - 2]] <- .edfL
        misc$earg[[Musual*ii - 1]] <- .edfR
        misc$earg[[Musual*ii - 0]] <- .esqw
    }

    misc$imethod <- .imethod
  }), list(
    .lloc = lloc, .lsca = lsca, .ldfL = ldfL, .ldfR = ldfR, .lsqw = lsqw,
    .eloc = eloc, .esca = esca, .edfL = edfL, .edfR = edfR, .esqw = esqw,
    .imethod = imethod ))),

  loglikelihood = eval(substitute(
    function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
    NOS <- extra$NOS
    Musual <- extra$Musual
    Loc <- eta2theta(eta[, Musual*(1:NOS)-4], .lloc, earg = .eloc)
    Sca <- eta2theta(eta[, Musual*(1:NOS)-3], .lsca, earg = .esca)
    DfL <- eta2theta(eta[, Musual*(1:NOS)-2], .ldfL, earg = .edfL)
    DfR <- eta2theta(eta[, Musual*(1:NOS)-1], .ldfR, earg = .edfR)
    Skw <- eta2theta(eta[, Musual*(1:NOS)-0], .lsqw, earg = .esqw)

 print("summary(DfL)")
 print( summary(DfL) )
 print("summary(DfR)")
 print( summary(DfR) )
 print("summary(Skw)")
 print( summary(Skw) )

 ## zedd <- (y - Loc) / Sca
    if (residuals) stop("loglikelihood residuals not ",
                        "implemented yet") else {
      sum(w * dtwopt(x = Loc, scale = Sca, skewpar = Skw,
                   dfL = DfL, dfR = DfR, log.arg = TRUE))
    }
  }, list( .lloc = lloc, .eloc = eloc,
           .lsca = lsca, .esca = esca,
           .ldfL = ldfL, .edfL = edfL,
           .ldfR = ldfR, .edfR = edfR,
           .lsqw = lsqw, .esqw = esqw ))),
  vfamily = c("agenstudentt"),
  deriv = eval(substitute(expression({
    Musual <- extra$Musual
    NOS  <- extra$NOS
    Loc <- eta2theta(eta[, Musual*(1:NOS) - 4], .lloc, earg = .eloc)
    Sca <- eta2theta(eta[, Musual*(1:NOS) - 3], .lsca, earg = .esca)
    DfL <- eta2theta(eta[, Musual*(1:NOS) - 2], .ldfL, earg = .edfL)
    DfR <- eta2theta(eta[, Musual*(1:NOS) - 1], .ldfR, earg = .edfR)
    Skw <- eta2theta(eta[, Musual*(1:NOS) - 0], .lsqw, earg = .esqw)

    dthetas.detas <- cbind(
      dloc.deta <- dtheta.deta(theta = Loc, .lloc, earg = .eloc),
      dsca.deta <- dtheta.deta(theta = Sca, .lsca, earg = .esca),
      ddfL.deta <- dtheta.deta(theta = DfL, .ldfL, earg = .edfL),
      ddfR.deta <- dtheta.deta(theta = DfR, .ldfR, earg = .edfR),
      dsqw.deta <- dtheta.deta(theta = Skw, .lsqw, earg = .esqw))
    dthetas.detas <-
      dthetas.detas[, interleave.VGAM(Musual * NOS, M = Musual)]

# Reference: first derivatives in Appendix C, pp.26--27
    zedd  <- (y - Loc) / Sca

    LL <- 1 + (zedd / (2 *      Skw  * Kkn.twopt(DfL)) )^2 / DfL
    LR <- 1 + (zedd / (2 * (1 - Skw) * Kkn.twopt(DfR)) )^2 / DfR
    tempL <- 1 - 1 / LL
    tempR <- 1 - 1 / LR
   
# Equation (56):
    dl.dsqw <-  (((DfL + 1) /      Skw ) * tempL) * (zedd <= 0) +
               -(((DfR + 1) / (1 - Skw)) * tempR) * (zedd >  0)
  
# Equation (57):
    dl.ddfL <- (-0.5 * log(LL) +
                 0.5 * (DfL + 1) * DDS.twopt(DfL) * tempL) * (zedd <= 0)
   
# Equation (58):
    dl.ddfR <- (-0.5 * log(LR) +
                 0.5 * (DfR + 1) * DDS.twopt(DfR) * tempR) * (zedd  > 0)
  
# Equation (59):
    dl.dloc <- 0.25 * ((DfL + 1) *  zedd / (LL * DfL *
               ((     Skw  * Kkn.twopt(DfL))^2 * Sca))) * (zedd <= 0) +
               0.25 * ((DfR + 1) *  zedd / (LR * DfR *
               (((1 - Skw) * Kkn.twopt(DfR))^2 * Sca))) * (zedd  > 0)
  
# Equation (60):
    dl.dsca <- ((DfL + 1) / Sca) * tempL * (zedd <= 0) +
               ((DfR + 1) / Sca) * tempR * (zedd  > 0) +
               -1 / Sca
    AAns <- c(w) * cbind(dl.dloc,
                         dl.dsca,
                         dl.ddfL,
                         dl.ddfR,
                         dl.dsqw) * dthetas.detas
 print("head(AAns) first deriv")
 print(head(AAns))
    c(w) * cbind(dl.dloc,
                 dl.dsca,
                 dl.ddfL,
                 dl.ddfR,
                 dl.dsqw) * dthetas.detas
   }), list(
    .lloc = lloc, .lsca = lsca, .ldfL = ldfL, .ldfR = ldfR, .lsqw = lsqw,
    .eloc = eloc, .esca = esca, .edfL = edfL, .edfR = edfR, .esqw = esqw ))),
  weight = eval(substitute(expression({
# Reference: Zhu & Galbraith (2010) page 301, eqn (29):

    wz   <- matrix(as.numeric(NA), n, dimm(M))

    edl.dsqw2   <-  3 * ((DfL + 1) / (     Skw  * (DfL + 3)) +
                         (DfR + 1) / ((1 - Skw) * (DfR + 3)))

    edl.ddfL2 <-    (DfL * DDS.twopt(DfL)^2 / (DfL + 3) -
                       2 * DDS.twopt(DfL)   / (DfL + 1) -
                     DDSp.twopt(DfL)) * Skw / 2

    edl.ddfR2 <-    (DfR * DDS.twopt(DfR)^2 / (DfR + 3) -
                       2 * DDS.twopt(DfR)   / (DfR + 1) -
                     DDSp.twopt(DfR)) * (1 - Skw) / 2

    edl.dloc2 <-    ((DfL + 1) / (     Skw  * (DfL + 3) * Kkn.twopt(DfL)^2 ) +
                     (DfR + 1) / ((1 - Skw) * (DfR + 3) * Kkn.twopt(DfR)^2 )
                    ) / (4 * Sca^2)

    edl.dsca2 <-  (     Skw  * DfL / (DfL + 3) +
                   (1 - Skw) * DfR / (DfR + 3)) * 2 / Sca^2

    edl.dsqwddfL <-  -1 / (DfL + 1) + DfL * DDS.twopt(DfL) / (DfL + 3)
    edl.dsqwddfR <-   1 / (DfR + 1) - DfR * DDS.twopt(DfR) / (DfR + 3)

    edl.dsqwdloc <-  -2 *  edl.dsqw2 / (3 * Sca)

    edl.dsqwdsca <- (DfL / (DfL + 3) - DfR / (DfR + 3)) * 2 /  Sca

    edl.ddfLdsca <- edl.dsqwddfL * Skw / Sca

    edl.ddfLdloc <-  (1 / (DfL + 1) -
                       (DfL + 1) * DDS.twopt(DfL) / (DfL + 3)) / Sca

    edl.ddfRdloc <- -(1 / (DfR + 1) -
                       (DfR + 1) * DDS.twopt(DfR) / (DfR + 3)) / Sca

    edl.ddfRdsca <-  -edl.dsqwddfR * (1 - Skw) / Sca


# zz; this may be in error (from equation (2))
    edl.dlocdsca <-  -edl.dsqwdsca * 2 / (3 * Sca)

# zz; this may be correct (from Appendix C)
    edl.dlocdsca <- 2 * (-(DfL + 1) / (DfL + 3) +
                          (DfR + 1) / (DfR + 3)) / Sca^2


    edl.ddfLddfR <- 0

    wz[, iam(1, 1, M)] <- edl.dloc2 # * dloc.deta^2
    wz[, iam(2, 2, M)] <- edl.dsca2 # * dsca.deta^2
    wz[, iam(3, 3, M)] <- edl.ddfL2 # * ddfL.deta^2
    wz[, iam(4, 4, M)] <- edl.ddfR2 # * ddfR.deta^2
    wz[, iam(5, 5, M)] <- edl.dsqw2 # * dsqw.deta^2

    wz[, iam(1, 2, M)] <-  edl.dlocdsca # * dsca.deta * dloc.deta
    wz[, iam(1, 3, M)] <-  edl.ddfLdloc # * ddfL.deta * dloc.deta
    wz[, iam(1, 4, M)] <-  edl.ddfRdloc # * ddfR.deta * dloc.deta
    wz[, iam(1, 5, M)] <-  edl.dsqwdloc # * dsqw.deta * dloc.deta

    wz[, iam(2, 3, M)] <-  edl.ddfLdsca # * ddfL.deta * dsca.deta
    wz[, iam(2, 4, M)] <-  edl.ddfRdsca # * ddfR.deta * dsca.deta
    wz[, iam(2, 5, M)] <-  edl.dsqwdsca # * dsqw.deta * dsca.deta

    wz[, iam(3, 4, M)] <-  edl.ddfLddfR
    wz[, iam(3, 5, M)] <-  edl.dsqwddfL # * dsqw.deta * ddfL.deta

    wz[, iam(4, 5, M)] <-  edl.dsqwddfR # * dsqw.deta * ddfR.deta


    myind = iam(NA, NA, M = Musual, both = TRUE, diag = TRUE)
    wz = wz * dthetas.detas[, myind$row.index] *
              dthetas.detas[, myind$col.index]


  if (ncol(cbind(y)) > 1)
    stop("code below might not work")



 if (FALSE) {
# Partition the working weight matrix according to LHS vs RHS.
    indexL = (zedd <= 0)
    indexR = (zedd >  0)
# Obsn on LHS means no information on dfR
    for (ii in 1:Musual)
      wz[indexL, iam(ii, 4, M = Musual)] <- 0
    wz[indexL, iam(4, 4, M = Musual)] <- .smallno

# Obsn on RHS means no information on dfL
    for (ii in 1:Musual)
      wz[indexR, iam(ii, 3, M = Musual)] <- 0
    wz[indexR, iam(3, 3, M = Musual)] <- .smallno
 }


#  For debugging:
  AAnsE <- c(w) * wz
 print("head(AAnsE) working weights")
 print(head(AAnsE))


    c(w) * wz
  }), list( .lloc = lloc, .eloc = eloc,
            .lsca = lsca, .esca = esca,
            .ldfL = ldfL, .edfL = edfL,
            .ldfR = ldfR, .edfR = edfR,
            .lsqw = lsqw, .esqw = esqw,
            .smallno = smallno ))))
}
