\name{TPNormD}

\alias{dtwopnorm}
\alias{ptwopnorm}
\alias{qtwopnorm}
\alias{rtwopnorm}

\title{Two-piece Normal Distribution }

\description{ Density, distribution function, quantile function and
random generation for the 2 and 3-parameter two-piece normal distribution
(TPND). The 2-parameter TPND may be used for quantile/expectile
regression.
}

\usage{
dtwopnorm(x, location = 0, scale = 1, skewpar = 0.5, log.arg = FALSE)
ptwopnorm(q, location = 0, scale = 1, skewpar = 0.5)
qtwopnorm(p, location = 0, scale = 1, skewpar = 0.5)
rtwopnorm(n, location = 0, scale = 1, skewpar = 0.5)
}

\arguments{
  
  \item{x, q}{  vector of quantiles.}
  \item{p   }{  vector of probabilities.}
  \item{n   }{  number of observations.
    If \code{length(n) > 1} then the length is taken to be the number
    required. }
  %\item{q}{vector of quantiles.}
  \item{location}{ location parameter. }
  \item{scale}{ positive scale parameter. }
  \item{skewpar}{ skew/asymmetric parameter, between zero and one. }
  \item{log.arg}{ logical. If \code{log.arg = TRUE} then the logarithm
    of the density is returned. }

}

\value{
  \code{dtwopnorm} gives the density,
  \code{ptwopnorm} gives the distribution function,
  \code{qtwopnorm} gives the quantile function, and
  \code{rtwopnorm} generates random deviates.
 }

\author{Arash Ardalan}

\details{ The \pkg{twopDist} family function for estimating the parameters, for the
formula of the probability density function and other details.
}

\seealso{
  \code{\link{TPNR}}.
}

\examples{

location = 1; scale = 1.2; skewpar = 0.25; nn = 201
x = seq(-1.5 , 6.5, len = nn)
y =  dtwopnorm(x, location, scale, skewpar)

\dontrun{
plot(x, y, type = "l", las = 1, ylim = c(0, 1.2), ylab =
     paste("twopnorm(location = ", location, ", scale = ", scale,
                  ", skewpar = ", skewpar, ")"),
     col = "blue", cex.main = 0.8,
     main = "Blue is density, red is cumulative distribution function",
     sub = "Purple lines are the 10, 20, ..., 90 percentiles")
lines(x, ptwopnorm(x, location, scale, skewpar), col = "red")
probs = seq(0.1, 0.9, by = 0.1)
Q = qtwopnorm(probs, location, scale, skewpar)
lines(Q, dtwopnorm(Q, location, scale, skewpar), col = "purple",
      lty = 3, type = "h")
lines(Q, ptwopnorm(Q, location, scale, skewpar), col = "purple",
      lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
abline(h = 0, col = "gray")
max(abs(ptwopnorm(Q, location, scale, skewpar) - probs)) # Should be 0
}
}
\keyword{distribution}
