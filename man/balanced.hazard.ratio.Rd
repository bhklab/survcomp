\name{balanced.hazard.ratio}
\alias{balanced.hazard.ratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to estimate the balanced hazard ratio through Cox regression}
\description{
Function to compute the balanced hazard ratio for a risk group prediction.
}
\usage{
balanced.hazard.ratio(x, surv.time, surv.event, alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of risk group predictions.}
  \item{surv.time}{a vector of event times.}
  \item{surv.event}{a vector of event occurrence indicators.}
  \item{alpha}{apha level to compute confidence interval.}
  \item{...}{additional parameters to be passed to the \code{\link{hazard.ratio}} and the \code{\link[survival]{coxph}} functions.}
}
\details{
  The balanced hazard ratio is computed using the Cox model.
}
\value{
  \item{balanced.hazard.ratio }{balanced hazard ratio estimate.}
	\item{coef }{coefficient (beta) estimated in the cox regression model.}
  \item{se }{standard error of the coefficient (beta) estimate.}
  \item{lower }{lower bound for the confidence interval.}
  \item{upper }{upper bound for the confidence interval.}
  \item{p.value }{p-value computed using the score (logrank) test whether the balanced hazard ratio is different from 1.}
  \item{n }{number of samples used for the estimation.}
  \item{coxm }{\code{\link[survival]{coxph.object}} fitted on the survival data  and \code{x}.}
  \item{data }{list of data used to compute the balanced hazard ratio (\code{x}, \code{surv.time} and \code{surv.event}).}
}
\references{Branders, S. and Dupont, P. (2015) "A balanced hazard ratio for risk group evaluation from survival data", \emph{Statistics in Medicine}, \bold{34}(17), pages 2528--2543. }
\author{Samuel Branders}
%\note{ ~~further notes~~ }
\seealso{\code{\link{hazard.ratio}}, \code{\link[survival]{coxph}}, \code{\link[survival]{coxph.object}}}
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
stime <- rexp(100)
cens   <- runif(100,.5,2)
sevent  <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
strat <- sample(1:3, 100, replace=TRUE)
weight <- runif(100, min=0, max=1)
balanced.hazard.ratio(x=age, surv.time=stime, surv.event=sevent, weights=weight,
  strat=strat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{univar}% __ONLY ONE__ keyword per line
