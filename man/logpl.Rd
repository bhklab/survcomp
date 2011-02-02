\name{logpl}
\alias{logpl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute the log partial likelihood of a Cox model }
\description{
  The function computes the log partial likelihood of a set of coefficients given some survival data.
}
\usage{
logpl(x, surv.time, surv.event, strata, beta, na.rm = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surv.time}{ vector of times to event occurrence }
  \item{surv.event}{ vector of indicators for event occurrence }
  \item{beta}{ vector of coefficients fitted by a Cox model for instance }
  \item{x}{ data matrix }
  \item{strata}{ stratification variable }
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise }
  \item{verbose}{ verbosity of the function }
}
%\details{}
\value{
vector of two elements: \code{logpl} and \code{event} for the estimation of the log partial likelihood and the number of events, respectively
}
\references{ Cox, D. R. (1972) "Regression Models and Life Tables", \emph{Journal of the Royal Statistical Society Series B}, \bold{34}, pages 187--220. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link[survival]{coxph}}, \code{\link{cvpl}} }
\examples{
require(survival)
set.seed(12345)
age <- rnorm(100, 50, 10)
stime <- rexp(100)
cens <- runif(100,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
##Cox model
coxm <- coxph(Surv(stime, sevent) ~ age)
##log partial likelihood of the null model
logpl(surv.time=stime, surv.event=sevent, beta=0, x=age)
##log partial likelihood of the Cox model
logpl(surv.time=stime, surv.event=sevent, beta=coxm$coefficients, x=age)
##equivalent to
coxm$loglik
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
