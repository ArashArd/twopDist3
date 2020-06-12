\name{mle_twoplaplace}
\alias{mle_twoplaplace}

\title{ The maximum likelihood estimation of the 3-parameter two-piece
Laplace distribution}

\description{
Maximum likelihood distribution for the 3-parameter two-piece Laplace distribution
(TPLapD). The 2-parameter TPLapD may be used for classical quantile
regression.
}
\usage{
mle_twoplaplace(x, k = 1000, e.rate = 1e-10, el.rate = length(x), max.iter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of data.}

\item{k}{
 number of skewpar segments.
}
  \item{e.rate}{
the tolerance of log-likelihood.
}
  \item{el.rate}{
 the tolerance of log-likelihood.
}
  \item{max.iter}{
the maximum of iteration
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{Est }{The MLE of parameters, obj@Est}
 \item{logLike }{The log-likelihood}
 \item{var }{The variance of MLEs}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Arash Ardalan}

\details{ See \code{\link{twopDist-class}}, the
\pkg{twopDist} family function for estimating the parameters, for the
formula of the probability density function and other details.
}

%\seealso{
%  \code{\link{twoplaplace3}}.
%}

\examples{
location = 2; scale = 1.5; skewpar = 0.25; nn = 201; set.seed(13)
y = rtwoplaplace(nn, location, scale, skewpar)

\dontrun{

MLE = mle_twoplaplace(y, k = 1000, e.rate = 1e-10, el.rate = length(y), max.iter = 1000)
summary(MLE, level = 0.95)
dplot(MLE, col = 'gray60')
}
}