\name{TPLapD}

\alias{dtwoplaplace}
\alias{ptwoplaplace}
\alias{qtwoplaplace}
\alias{rtwoplaplace}

\title{Two-piece Laplace Distribution }

\description{ Density, distribution function, quantile function and
random generation for the 3-parameter two-piece Laplace distribution
(TPLapD). The 2-parameter TPLapD may be used for classical quantile
regression. 
}

\usage{
dtwoplaplace(x, location = 0, scale = 1, skewpar = 0.5, log.arg = FALSE)
ptwoplaplace(q, location = 0, scale = 1, skewpar = 0.5)
qtwoplaplace(p, location = 0, scale = 1, skewpar = 0.5)
rtwoplaplace(n, location = 0, scale = 1, skewpar = 0.5)
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
  \item{skewpar}{ skew/asymmetric parameter, between zero and one.}
  \item{log.arg}{ Logical. If \code{log.arg = TRUE} then the logarithm
        of the density is returned. }

}

\value{
  \code{dtwoplaplace}  gives the density,
  \code{ptwoplaplace}  gives the distribution function,
  \code{qtwoplaplace}  gives the quantile function, and
  \code{rtwoplaplace}  generates random deviates.
 }

\author{Arash Ardalan}

%\details{ See \code{\link{twoplaplace3}} and \code{\link{twoplaplace}}, the
%\pkg{twopDist} family function for estimating the parameters, for the
%formula of the probability density function and other details.
%}

%\seealso{
%  \code{\link{twoplaplace3}}.
%}

\examples{
location = 0; scale = 1.5; skewpar = 0.25; nn = 201
x = seq(-3.5 ,6.5, len = nn)
y = dtwoplaplace(x, location, scale, skewpar)

\dontrun{
plot(x, y, type = "l", las = 1, ylim = c(0, 1.2), ylab =
     paste("ftwoplaplace(location = ", location, ", scale = ", scale,
                  ", skewpar = ", skewpar, ")"),
     col  = "blue", cex.main = 0.8,
     main = "Blue is density, red is cumulative distribution function",
     sub  = "Purple lines are the 10, 20, ..., 90 percentiles")
lines(x, ptwoplaplace(x, location, scale, skewpar), col = "red")
probs = seq(0.1, 0.9, by = 0.1)
Q = qtwoplaplace(probs, location, scale, skewpar)
lines(Q, dtwoplaplace(Q, location, scale, skewpar), col = "purple",
      lty = 3, type = "h")
lines(Q, ptwoplaplace(Q, location, scale, skewpar), col = "purple",
      lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
abline(h = 0, col = "gray")
max(abs(ptwoplaplace(Q, location, scale, skewpar) - probs)) # Should be 0
}
}
\keyword{distribution}
