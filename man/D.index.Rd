\name{D.index}
\alias{D.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute the D index}
\description{
Function to compute the D index for a risk prediction, i.e. an estimate of the log hazard ratio comparing two equal-sized prognostic groups. This is a natural measure of separation between two independent survival distributions under the proportional hazards assumption.
}
\usage{
D.index(x, surv.time, surv.event, weights, strat, alpha = 0.05, method.test = c("logrank", "likelihood.ratio", "wald"), na.rm = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of risk predictions.}
  \item{surv.time}{a vector of event times.}
  \item{surv.event}{a vector of event occurrence indicators.}
  \item{weights}{weight of each sample.}
  \item{strat}{stratification indicator.}
  \item{alpha}{apha level to compute confidence interval.}
  \item{method.test}{Statistical test to use in order to compute the p-values related to a D. index, see \link[survival]{summary.coxph} for more details.}
  \item{na.rm}{\code{TRUE} if missing values should be removed.}
  \item{...}{additional parameters to be passed to the \code{\link[survival]{coxph}} function.}
}
\details{
  The D index is computed using the Cox model fitted on the scaled rankits of the risk scores instead of the risk scores themselves. The scaled rankits are the expected standard Normal order statistics scaled by \code{kappa = sqrt(8/pi)}. See (Royston and Sauerbrei, 2004) for details.

Note that the value D reported in (Royston and Sauerbrei, 2004) is given 
}
\value{
  \item{d.index }{D index estimate.}
  \item{se }{standard error of the estimate.}
  \item{lower }{lower bound for the confidence interval.}
  \item{upper }{upper bound for the confidence interval.}
  \item{p.value }{p-value for the statistical test if the estimate if different from 0.5.}
  \item{n }{number of samples used for the estimation.}
  \item{coxm }{\code{\link[survival]{coxph.object}} fitted on the survival data  and \code{z} (see below).}
  \item{data }{list of data used to compute the index (\code{x}, \code{z}, \code{surv.time} and \code{surv.event}). The item \code{z} contains the scaled rankits which are the expected standard Normal order statistics scaled by \code{kappa}.}
}
\references{Royston, P. and Sauerbrei, W. (2004) "A new measure of prognostic separation in survival data", Statistics in Medicine, \bold{23}, pages 723--748.}
\author{Benjamin Haibe-Kains}
%\note{ ~~further notes~~ }
\seealso{\code{\link[survival]{coxph}}, \code{\link[survival]{coxph.object}}, \code{\link[SuppDists]{normOrder}}}
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
stime <- rexp(100)
cens   <- runif(100,.5,2)
sevent  <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
strat <- sample(1:3, 100, replace=TRUE)
weight <- runif(100, min=0, max=1)
D.index(x=age, surv.time=stime, surv.event=sevent, weights=weight, strat=strat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{univar}% __ONLY ONE__ keyword per line
