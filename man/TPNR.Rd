\name{TPNR}
\alias{TPNR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
two-pice normal quantile model
}
\description{
TPNR is used to fit two-pice normal quantile models.
It can be used to carry out regression quuantile taype (meantile) at a specefic quantile.
}
\usage{
TPNR(formula, data, weights, method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE, contrasts = NULL, offset, p, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%
an object of class \code{formula}.
}
\item{data}{
  an optional data frame, list or environment (or object coercible by
\code{as.data.frame} to a data frame) containing the variables in the model.
If not found in data, the variables
are taken from environment(formula), typically the environment
from which lm is called.
 }

\item{method}{
a maximum likelihood method.
}

\item{myskew}{
  a specific quantile to be estimated, this is a number  between 0 and 1.
}
\item{myscale}{
an optional value, the scale parameter of the model.
}
  \item{\dots}{
...
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
TPNLR returns an object of class "TPNLR".

\item{coefficients }{a vector of coefficients}
\item{residuals }{the residuals, that is response minus fitted values.}
\item{fitted.values}{the fitted values.}

...
}
\references{
Arash Ardalan manuscripts.
}
\author{
Arash Ardalan
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(swiss)
summary(lm(Fertility ~ . , data = swiss))
TPNL_model_1 = TPNR(Fertility ~ . , data = swiss, p = 0.5)
TPNL_model_2 = TPNR(Fertility ~ . , data = swiss, p = 0.85)
summary(TPNL_model_1)
summary(TPNL_model_2)


}
