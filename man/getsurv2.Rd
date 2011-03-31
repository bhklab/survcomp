\name{getsurv2}
\alias{getsurv2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to retrieve the survival probabilities at a specific point in time }
\description{
  The function retrieves the survival probabilities from a survfit object, for a specific point in time.
}
\usage{
getsurv2(sf, time, which.est = c("point", "lower", "upper"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf}{ \code{survfit} object }
  \item{time}{ time at which the survival probabilities must be retrieved }
  \item{which.est}{ which estimation to be returned? \code{point} for the point estimate, \code{lower} for the lower bound and \code{upper} for the upper bound }
}
\details{ The survival probabilities are estimated through the \code{\link[survival]{survfit}} function.}
\value{
  vector of survival probabilities
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link[survival]{survfit}} }
\examples{
require(survival)
set.seed(12345)
age <- rnorm(30, 50, 10)
stime <- rexp(30)
cens <- runif(30,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
sf <- survfit(Surv(stime, sevent) ~ 1)
getsurv2(sf, time=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
