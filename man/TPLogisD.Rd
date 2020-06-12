\name{TPLogisD}

\alias{dtwoplogis}
\alias{ptwoplogis}
\alias{qtwoplogis}
\alias{rtwoplogis}

\title{Two-piece Logistic Distribution }

\description{ Density, distribution function, quantile function and
random generation for the 3-parameter two-piece logistic distribution
(TPLogisD). The 2-parameter TPLogisD may be used for quantile/expectile
type regression. This family also is more efficient and more robust in
comparison with classical quantile regression.
}

\usage{
dtwoplogis(x, location = 0, scale = 1, skewpar = 0.5, log.arg = FALSE)
ptwoplogis(q, location = 0, scale = 1, skewpar = 0.5)
qtwoplogis(p, location = 0, scale = 1, skewpar = 0.5)
rtwoplogis(n, location = 0, scale = 1, skewpar = 0.5)
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
  \item{log.arg}{ Logical. If \code{log.arg = TRUE} then the logarithm
    of the density is returned. }

}

\value{
  \code{dtwoplogis} gives the density,
  \code{ptwoplogis} gives the distribution function,
  \code{qtwoplogis} gives the quantile function, and
  \code{rtwoplogis} generates random deviates.
 }

\author{Arash Ardalan}

\details{The \pkg{twopDist} family function for estimating the parameters, for the
formula of the probability density function and other details.
}

\seealso{
  \code{\link{TPLogR}}.
}

\examples{

location = 1; scale = 1.2; skewpar = 0.25; nn = 201
  # location = 1, scale  = 1.2;
x = seq(-1.5 , 6.5, len=nn)
y =  dtwoplogis(x, location, scale, skewpar)

\dontrun{
plot(x, y, type="l", las=1, ylim = c(0, 1.2), ylab =
     paste("ftwoplogis(location=", location, ", scale=", scale,
                  ", skewpar=", skewpar, ")"),
     col="blue", cex.main = 0.8,
     main="Blue is density, red is cumulative distribution function",
     sub="Purple lines are the 10, 20, ..., 90 percentiles")
lines(x, ptwoplogis(x, location, scale, skewpar), col="red")
probs = seq(0.1, 0.9, by = 0.1)
Q = qtwoplogis(probs, location, scale, skewpar)
lines(Q, dtwoplogis(Q, location, scale, skewpar), col = "purple",
        lty = 3, type = "h")
lines(Q, ptwoplogis(Q, location, scale, skewpar), col = "purple",
         lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
abline(h = 0,col= "gray")
max(abs(ptwoplogis(Q, location, scale, skewpar) - probs)) # Should be 0
}
}
\keyword{distribution}
