\name{TPNormLapD}

\alias{dtwopnormlap}
\alias{ptwopnormlap}
\alias{qtwopnormlap}
\alias{rtwopnormlap}

\title{Two-piece Laplace Distribution }

\description{ Density, distribution function, quantile function and
random generation for the 3-parameter two-piece normal Laplace distribution
(TPNormLapD). This distribution from one side is short tail(normal) and in
the other side(Laplace) is long tail, see Ardalan et. All (2012).
The 2-parameter TPNormLapD may be used for quantile type
regression model. 

}

\usage{
dtwopnormlap(x, location = 0, scale = 1, skewpar = 0.5, log.arg = FALSE)
ptwopnormlap(q, location = 0, scale = 1, skewpar = 0.5)
qtwopnormlap(p, location = 0, scale = 1, skewpar = 0.5)
rtwopnormlap(n, location = 0, scale = 1, skewpar = 0.5)
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
  \code{dtwopnormlap}  gives the density,
  \code{ptwopnormlap}  gives the distribution function,
  \code{qtwopnormlap}  gives the quantile function, and
  \code{rtwopnormlap}  generates random deviates.
 }

 \references{Arash Ardalan, S. M. Sadooghi-Alvandi and A.N. Nematollahi
  (2012), The Two-Piece Normal-Laplace Distribution, 35 pages (Comm. in
  Statist.: Theory & Methods. ).
}
\author{Arash Ardalan}

%\details{ See \code{\link{twopnormlap3}} and \code{\link{twopnormlap}}, the
%\pkg{twopDist} family function for estimating the parameters, for the
%formula of the probability density function and other details.
%}

%\seealso{
%  \code{\link{twopnormlap3}}.
%}

\examples{

location = 0; scale = 1.5; skewpar = 0.5; nn = 201
x = seq(-3.5 , 6.5, len = nn)
y = dtwopnormlap(x, location, scale, skewpar)

\dontrun{
plot(x, y, type = "l", las = 1, ylim = c(0, 1.2), ylab =
     paste("ftwopnormlap(location=", location, ", scale=", scale,
                  ", skewpar=", skewpar, ")"),
     col  = "blue", cex.main = 0.8,
     main = "Blue is density, red is cumulative distribution function",
     sub  = "Purple lines are the 10, 20, ..., 90 percentiles")
lines(x, ptwopnormlap(x, location, scale, skewpar), col = "red")
probs = seq(0.1, 0.9, by = 0.1)
Q = qtwopnormlap(probs, location, scale, skewpar)
lines(Q, dtwopnormlap(Q, location, scale, skewpar), col = "purple",
      lty = 3, type = "h")
lines(Q, ptwopnormlap(Q, location, scale, skewpar), col = "purple",
      lty = 3, type = "h")
abline(h = probs, col = "purple", lty = 3)
abline(h = 0,col = "gray")
max(abs(ptwopnormlap(Q, location, scale, skewpar) - probs)) # Should be 0
}
}
\keyword{distribution}
