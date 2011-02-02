\name{score2proba}
\alias{score2proba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute the survival probabilities from a risk score }
\description{
  the function uses either a Cox model or the product-limit estimator to compute the survival probabilities from a risk score for a specific point in time.
}
\usage{
score2proba(data.tr, score, yr, method = c("cox", "prodlim"), conf.int = 0.95, which.est= c ("point", "lower", "upper"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.tr}{ the data frame for the training set. This data frame must contain three columns for the times, the event occurrence and the risk score. These columns are called "time", "event" and "score" respectively }
  \item{score}{ risk score for the test set }
  \item{yr}{ a point in time for which the survival probabilities must be computed }
  \item{method}{ method for survival probabilities estimation, either \code{cox} or \code{prodlim} for the Cox model or the product-limit estimator, respectively }
  \item{conf.int}{ value in [0,1]. Default at 0.95 }
  \item{which.est}{ which estimation to be returned? \code{point} for the point estimate, \code{lower} for the lower bound and \code{upper} for the upper bound }
}
%\details{}
\value{
  vector of predicted survival probabilities
}
\references{ Cox, D. R. (1972) "Regression Models and Life Tables", \emph{Journal of the Royal Statistical Society Series B}, \bold{34}, pages 187--220.

Andersen, P. K. and Borgan, O. and Gill, R. D. and Keiding, N. (1993) "Statistical Models Based on Counting Processes", \emph{Springer}. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link[survival]{coxph}}, \code{\link[prodlim]{prodlim}} }
\examples{
set.seed(12345)
age <- rnorm(30, 50, 10)
stime <- rexp(30)
cens <- runif(30,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
dd <- data.frame("time"=stime, "event"=sevent, "score"=age)

#Cox's model
score2proba(data.tr=dd, score=dd$score, yr=1, method="cox")
#product-limit estimator
score2proba(data.tr=dd, score=dd$score, yr=1, method="prodlim")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
