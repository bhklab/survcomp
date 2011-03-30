\name{logpl}
\alias{logpl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute the log partial likelihood of a Cox model }
\description{
  The function computes the log partial likelihood of a set of coefficients given some survival data.
}
\usage{
logpl(pred, surv.time, surv.event, strata, na.rm = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surv.time}{ vector of times to event occurrence }
  \item{surv.event}{ vector of indicators for event occurrence }
  \item{pred}{ linear predictors computed using the Cox model }
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
dd <- data.frame("stime"=stime, "sevent"=sevent, "age"=age)
##Cox model
coxm <- coxph(Surv(stime, sevent) ~ age, data=dd)
##log partial likelihood of the null model
logpl(pred=rep(0, nrow(dd)), surv.time=stime, surv.event=sevent)
##log partial likelihood of the Cox model
logpl(pred=predict(object=coxm, newdata=dd), surv.time=stime, surv.event=sevent)
##equivalent to
coxm$loglik
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
