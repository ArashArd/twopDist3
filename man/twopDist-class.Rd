\name{twopDist-class}
\Rdversion{1.1}
\docType{class}
\alias{twopDist-class}
\alias{twopDist}
\alias{ConfIint,twopDist-method}
\alias{plot,twopDist,ANY-method}
\alias{show,twopDist-method}
\alias{summary,twopDist-method}

\title{Class \code{"twopDist"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{twopDist(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Est}:}{Object of class \code{"vector"}. Estimated parameters.}
    \item{\code{logLike}:}{Object of class \code{"numeric"}. Extract maximized log-likelihood.  }
    \item{\code{var}:}{Object of class \code{"vector"}. Approximate variances. }
    \item{\code{data}:}{Object of class \code{"numeric"}. The observations. }
    \item{\code{Fitted}:}{Object of class \code{"numeric"}. The fitted values.}
    \item{\code{KST}:}{Object of class \code{"vector"}. The D statistics and the P-value of K--S test.  }
  }
}
\section{Methods}{
  \describe{
    \item{ConfIint}{\code{signature(object = "twopDist")}: Confidence intervals from likelihood profiles. }
    \item{plot}{\code{signature(x = "twopDist", y = "ANY")}: Plot fitted distribution with histogram.}
    \item{show}{\code{signature(object = "twopDist")}: Display object briefly.}
    \item{summary}{\code{signature(object = "twopDist")}: Generate object summary.}
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
Arash Ardalan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mle_twoplaplace}}
}
\examples{
showClass("twopDist")
}
\keyword{classes}
\usage{
twopDist(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
