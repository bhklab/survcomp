\name{hazard.ratio}
\alias{hazard.ratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to estimate the hazard ratio through Cox regression}
\description{
Function to compute the hazard ratio for a risk prediction.
}
\usage{
hazard.ratio(x, surv.time, surv.event, weights, strat, alpha = 0.05, na.rm = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of risk predictions.}
  \item{surv.time}{a vector of event times.}
  \item{surv.event}{a vector of event occurrence indicators.}
  \item{weights}{weight of each sample.}
  \item{strat}{stratification indicator.}
  \item{alpha}{apha level to compute confidence interval.}
  \item{na.rm}{\code{TRUE} if missing values should be removed.}
  \item{...}{additional parameters to be passed to the \code{\link[survival]{coxph}} function.}
}
\details{
  The hazard ratio is computed using the Cox model.
}
\value{
  \item{hazard.ratio }{hazard ratio estimate.}
  \item{se }{standard error of the estimate.}
  \item{lower }{lower bound for the confidence interval.}
  \item{upper }{upper bound for the confidence interval.}
  \item{p.value }{p-value for the statistical test if the estimate if different from 0.5.}
  \item{n }{number of samples used for the estimation.}
  \item{coxm }{\code{\link[survival]{coxph.object}} fitted on the survival data  and \code{x} (see below).}
  \item{data }{list of data used to compute the hazard ratio (\code{x}, \code{surv.time} and \code{surv.event}).}
}
\references{Cox, D. R. (1972) "Regression Models and Life Tables", \emph{Journal of the Royal Statistical Society Series B}, \bold{34}, pages 187--220. }
\author{Benjamin Haibe-Kains}
%\note{ ~~further notes~~ }
\seealso{\code{\link[survival]{coxph}}, \code{\link[survival]{coxph.object}}}
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
stime <- rexp(100)
cens   <- runif(100,.5,2)
sevent  <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
strat <- sample(1:3, 100, replace=TRUE)
weight <- runif(100, min=0, max=1)
hazard.ratio(x=age, surv.time=stime, surv.event=sevent, weights=weight, strat=strat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{univar}% __ONLY ONE__ keyword per line
