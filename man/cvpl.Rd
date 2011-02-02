\name{cvpl}
\alias{cvpl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute the CVPL }
\description{
  The function computes the cross-validated partial likelihood (CVPL) for the Cox model.
}
\usage{
cvpl(x, surv.time, surv.event, strata, nfold = 1, setseed, na.rm = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ data matrix }
  \item{surv.time}{ vector of times to event occurrence }
  \item{surv.event}{ vector of indicators for event occurrence }
  \item{strata}{ stratification variable }
  \item{nfold}{ number of folds for the cross-validation }
  \item{setseed}{ seed for the random generator }
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise }
  \item{verbose}{ verbosity of the function }
}
%%\details{
%%}
\value{
  \item{cvpl }{mean cross-validated partial likelihood}
  \item{pl }{vector of cross-validated partial likelihoods}
  \item{convergence }{vector of booleans reporting the convergence of the Cox model in each fold}
  \item{n }{number of observations used to estimate the cross-validated partial likelihood}
}
\references{
Verweij PJM. and van Houwelingen H (1993) "Cross-validation in survival analysis", \emph{Statistics in Medicine}, \bold{12}, pages 2305--2314

van Houwelingen H, Bruinsma T, Hart AA, van't Veer LJ, and Wessels LFA (2006) "Cross-validated Cox regression on microarray gene expression data", \emph{Statistics in Medicine}, \bold{25}, pages 3201--3216.
}
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link{logpl}}, \code{\link[survival]{coxph}} }
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
stime <- rexp(100)
cens   <- runif(100,.5,2)
sevent  <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
strat <- sample(1:3, 100, replace=TRUE)
cvpl(x=age, surv.time=stime, surv.event=sevent, strata=strat, nfold=10, setseed=54321)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
