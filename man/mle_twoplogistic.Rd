\name{mle_twoplogistic}
\alias{mle_twoplogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The maximum likelihood estimation of the 3-parameter two-piece logistic quantile distribution
}
\description{
Maximum likelihood distribution for the 3-parameter two-piece logistic quantile distribution
(TPLapD). The 2-parameter two-piece logistic quantile distribution may be used for quantile type
regression.
}
\usage{
mle_twoplogistic(x, maxit = 1000, intialVal = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of data.
}
  \item{maxit}{
the maximum of iteration
}
  \item{intialVal}{
The intial value for starting the algorithm.
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
y = rtwoplogis(nn, location, scale, skewpar)

\dontrun{
 InV = c(mean(y), sd(y), 0.3)
MLE = mle_twoplogistic(y, maxit = 1000, intialVal = InV)
summary(MLE, level = 0.95)
dplot(MLE, col = 'gray60')
}
}
