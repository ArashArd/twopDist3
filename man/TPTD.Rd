\name{TPTD}

\alias{dtwopt}
\alias{ptwopt}
\alias{qtwopt}
\alias{rtwopt}

\title{Two-piece Student's t-distribution}

\description{ Density, distribution function, quantile function and
random generation for the 5-parameter two-piece Student's t-distribution
(TPTD). The 2-parameter TPTD (with given skewpar, dfL and dfR ) may be
used for quantile type regression. 
}

\usage{
dtwopt(x, location = 0, scale = 1, skewpar = 0.5, dfL = 2, dfR = 2,
       log.arg = FALSE)
ptwopt(q, location = 0, scale = 1, skewpar = 0.5, dfL = 2, dfR = 2)
qtwopt(p, location = 0, scale = 1, skewpar = 0.5, dfL = 2, dfR = 2)
rtwopt(n, location = 0, scale = 1, skewpar = 0.5, dfL = 2, dfR = 2)
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
  \item{dfL}{the left degree of freedom.}
  \item{dfR}{the right degree of freedom.}
  

}

\value{
  \code{dtwopt}  gives the density,
  \code{ptwopt}  gives the distribution function,
  \code{qtwopt}  gives the quantile function, and
  \code{rtwopt}  generates random deviates.
 }

\author{Arash Ardalan}

%\details{ See \code{\link{twopt3}} and \code{\link{twopt}}, the
%\pkg{twopDist} family function for estimating the parameters, for the
%formula of the probability density function and other details.
%}

%\seealso{
%  \code{\link{twopt3}}.
%}

\examples{

location = 0; scale = 1.5; skewpar = 0.25; nn = 201
x = seq(-3.5 ,6.5, len = nn)
y = dtwopt(x, location, scale, skewpar)

\dontrun{
plot(x, y, type = "l", las = 1, ylim = c(0, 1.2), ylab =
     paste("ftwopt(location = ", location, ", scale = ", scale,
           ", dfL = dfR = ", 2, ", skewpar = ", skewpar, ")"),
     col  = "blue", cex.main = 0.8,
     main = "Blue is density, red is cumulative distribution function",
     sub  = "Purple lines are the 10, 20, ..., 90 percentiles")
lines(x, ptwopt(x, location, scale, skewpar), col = "red")

probs = seq(0.1, 0.9, by = 0.1)
Q = qtwopt(probs, location, scale, skewpar)
lines(Q, dtwopt(Q, location, scale, skewpar), col = "purple",
      lty = 3, type = "h")
lines(Q, ptwopt(Q, location, scale, skewpar), col = "purple",
      lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
abline(h = 0, col = "gray")
max(abs(ptwopt(Q, location, scale, skewpar) - probs)) # Should be 0
}
}
\keyword{distribution}
