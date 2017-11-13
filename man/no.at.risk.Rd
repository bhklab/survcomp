\name{no.at.risk}
\alias{no.at.risk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute the number of individuals at risk }
\description{
Function to compute the number of individuals at risk at certain time points, as used in the Kaplan-Meier estimator for instance, depending on stratification.
}
\usage{
no.at.risk(formula.s, data.s, sub.s = "all", t.step, t.end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.s}{formula composed of a \code{Surv} object and a strata variable (i.e. stratification).}
  \item{data.s}{data frame composed of the variables used in the formula.}
  \item{sub.s}{vector of booleans specifying if only a subset of the data should be considered.}
  \item{t.step}{time step at which the number of individuals at risk is computed.}
  \item{t.end}{maximum time to be considered.}
}
\details{
  The original version of this function was kindly provided by Dr Christos Hatzis (January, 17th 2006).
}
\value{
number of individuals at risk at each time step specified in \code{t.step} up to \code{t.end}.
}
%\references{}
\author{Christos Hatzis, Benjamin Haibe-Kains}
%\note{ ~~further notes~~ }
\seealso{\code{\link[survival]{survfit}}, \code{\link[survcomp]{km.coxph.plot}}}
\examples{
set.seed(12345)
stime <- rexp(100)
cens   <- runif(100,.5,2)
sevent  <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
strat <- sample(1:3, 100, replace=TRUE)
dd <- data.frame("surv.time"=stime, "surv.event"=sevent, "strat"=strat)
no.at.risk(formula.s=Surv(surv.time,surv.event) ~ strat, data.s=dd,
  sub.s="all", t.step=0.05, t.end=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{univar}% __ONLY ONE__ keyword per line
