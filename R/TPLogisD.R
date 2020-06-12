#LOGISTIC FAMILY modiffied 04 2020 by arash


#=============================================================
#           dtwoplogis
#=============================================================

sech2 <- function(x) 1 / (cosh(x))^2

dtwoplogis <- function(x, location = 0, scale = 1, skewpar = 0.5,
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
    s1 <-  (2 * skewpar) * dlogis(zedd, scale = 2 * skewpar)
    #sech2(zedd / (4 * skewpar))
    s2 <- (2 * ( 1 - skewpar)) * dlogis(zedd, scale = 2 * (1-skewpar))
    #sech2(zedd / (4 * (1 - skewpar)))
    dens <- s1
    dens[zedd > 0] <- s2[zedd > 0]
    dens <- 1 / (1 * scale) * dens

  if (log.arg) log(dens) else dens
}



#=============================================================
#         ptwoplogis
#=============================================================


ptwoplogis <- function(q, location = 0, scale = 1, skewpar = 0.5) {
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
    p2 <- (1-skewpar)
  s1 <-  2 * skewpar * plogis(zedd, scale = 2 * skewpar) #/ scale
  s2 <- skewpar +
        2 * p2 *(plogis(zedd, scale = 2 * p2) - plogis(0, scale = 2 * p2))
 ans <- rep(0.0, length(zedd))
 ans[zedd <= 0] <- s1[zedd <= 0]
 ans[zedd > 0] <- s2[zedd > 0]
  ans
}


#============================================================
#         qtwoplogis
#============================================================

pos <- function(x) ifelse(x > 0, x, 0.0)
lone <- function(x) ifelse(x <= 1, x, 0.0)

qtwoplogis <- function(p, location = 0, scale = 1, skewpar = 0.5){

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

  #LLL_L <-  sum(pp <= skewpar)
  #qtwop <- rep(as.numeric(NA), LLL)
   QL = pp / (2 * skewpar); QL = ifelse(QL >= 1 | QL <=0, .1, QL)
   QR = (pp - skewpar) / (2 * (1-skewpar)) +
        plogis(0, scale = 2 * (1-skewpar))
        QR = ifelse(QR >= 1 | QR <=0, .1, QR)
  qtwopL <- qlogis(QL, location = 0, scale = 2 * skewpar) * (pp <= skewpar)
  qtwopR <- qlogis(QR , location = 0, scale = 2 * (1-skewpar)) * (pp > skewpar)

  qtwop <- qtwopL + qtwopR
   qtwop * scale + location
}


#=============================================================
#         rtwoplogis
#=============================================================

rtwoplogis <- function(n, location = 0, scale = 1, skewpar = 0.5) {
  qtwoplogis(p = runif(n), location = location, scale = scale,
          skewpar = skewpar)
}

#=============================================================
#          EIM  twoplogis
#=============================================================
EIMtplog = function(myscale, myskew){
EIM = diag(3)
temp10 <- myskew * (1 - myskew)
    EIM[1,1]      <- 1 / (12 * temp10 * myscale^2)
  EIM[2,2]         <-  (1/3 + (pi / 3)^2) /  myscale^2
  EIM[3,3]        <- (4/3 + (pi / 3)^2) / temp10
  EIM[1,3]   <- -((1/6) + (log(2)/3)) / (temp10 * myscale)

 EIM[3,1]  = EIM[1,3]
EIM}

#=============================================================
#         Inverse of EIM  twoplogis
#=============================================================
IEIMtplog = function(myscale, myskew){
EIM = diag(3)
p = myskew
sigma = myscale
Pi = pi

 EIM[1,1]  = 12*sigma^2*(4+Pi)*p*(-1+p)/(4*log(2)^2+4*log(2)-Pi-3)

 EIM[1,3]  = -(6+12*log(2))*p*sigma*(-1+p)/(4*log(2)^2+4*log(2)-Pi-3)
 EIM[3,1]  =  EIM[1,3]
 EIM[2,2]  = 9*sigma^2/(Pi^2+3)

 EIM[3,3]  = 3*p*(-1+p)/(4*log(2)^2+4*log(2)-Pi-3)
EIM}

#=============================================================
#          First derivative  twoplogis
#=============================================================

Utplog = function(y, mylocat, myscale, myskew){

 zedd <- (y - mylocat) / myscale # Is a n x NOS matrix

        #  cond1 <-    (zedd <= 0)
          cond2 <-    (zedd > 0)

      dl.dlocat        <- 2 * tanh(zedd / (4 * myskew)) / (myskew)  # cond1
      dl.dlocat[cond2] <- (2 * tanh(zedd / (4 * (1-myskew))) / ( 1-myskew))[cond2]
      dl.dlocat        <- dl.dlocat / (4 * myscale)

      dl.dscale        <-  2 * tanh(zedd / (4 * myskew)) * zedd / (myskew)
      dl.dscale[cond2] <- (2 * tanh(zedd / (4 * (1-myskew))) * zedd / (1 - myskew))[cond2]
      dl.dscale        <- (-1 + dl.dscale / 4) / myscale

      dl.dskewpar        <-  tanh(zedd / (4 * myskew)) * zedd / (2 * myskew^2)
      dl.dskewpar[cond2] <- (-tanh(zedd / (4 * (1-myskew))) *
                                zedd / (2 * (1-myskew)^2))[cond2]

U = c(sum(dl.dlocat), sum(dl.dscale), sum(dl.dskewpar))
U}



