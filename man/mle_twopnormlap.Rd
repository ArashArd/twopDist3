\name{mle_twopnormlap}
\alias{mle_twopnormlap}

\title{
The maximum likelihood estimation of the 3-parameter two-piece normal
Laplace quantile distribution
}
\description{
Maximum likelihood distribution for the 3-parameter
 two-piece normal Laplace distribution.
}
\usage{
mle_twopnormlap(x, e = 1e-10, max_iter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of data.
}
  \item{max_iter}{
the maximum of iteration
}
  \item{e}{
The tolerance rate.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{Est }{The MLE of parameters, obj@Est}
 \item{logLike }{The log-likelihood}
 \item{var }{The variance of MLEs}
}
\references{
Arash Ardalan manuscripts.
}
\author{
Arash Ardalan
}
\note{
%%  ~~further notes~~
}





\examples{
location = 2; scale = 1.5; skewpar = 0.35; nn = 120; set.seed(13)
y = rtwopnormlap(nn, location, scale, skewpar)

\dontrun{
 InV = c(mean(y), sd(y), 0.3)
MLE = mle_twopnormlap(y, max_iter = 100)
summary(MLE, level = 0.95)
dplot(MLE, col = 'gray60')
}
}
