\name{sbrier.score2proba}
\alias{sbrier.score2proba}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute the BSCs from a risk score, for all the times of event occurrence }
\description{
  The function computes all the Brier scores (BSC) and the corresponding integrated Briser score (IBSC) from a risk score, for all the times of event occurrence. The risk score is first transformed in survival probabilities using either a Cox model or the product-limit estimator.
}
\usage{
sbrier.score2proba(data.tr, data.ts, method = c("cox", "prodlim"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.tr}{ the data frame for the training set. This data frame must contain three columns for the times, the event occurrence and the risk score. These columns are called "time", "event" and "score" respectively. }
  \item{data.ts}{ the data frame for the test set. This data frame must contain three columns for the times, the event occurrence and the risk score. These columns are called "time", "event" and "score" respectively. }
  \item{method}{ method for survival probabilities estimation using either a Cox model or the product-limit estimator }
}
%\details{}
\value{
  \item{time }{ vector of points in time}
  \item{bsc }{vector of Brier scores (BSC) at ome points in time}
  \item{bsc.integrated }{value of the integrated Brier score (IBSC)}
}
\references{ Brier, G. W. (1950) "Verification of forecasts expressed in terms of probabilities", \emph{Monthly Weather Review}, \bold{78}, pages 1--3.
	
	Graf, E. and Schmoor, C. and Sauerbrei, W. and Schumacher, M. (1999) "Assessment and comparison of prognostic classification schemes for survival data ", \emph{Statistics in Medicine}, \bold{18}, pages 2529--2545.
	
	Cox, D. R. (1972) "Regression Models and Life Tables", \emph{Journal of the Royal Statistical Society Series B}, \bold{34}, pages 187--220.
	
	Andersen, P. K. and Borgan, O. and Gill, R. D. and Keiding, N. (1993) "Statistical Models Based on Counting Processes", \emph{Springer}. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link[ipred]{sbrier}}, \code{\link[survival]{coxph}}, \code{\link[prodlim]{prodlim}} }
\examples{
set.seed(12345)
age <- rnorm(30, 50, 10)
stime <- rexp(30)
cens <- runif(30,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
dd <- data.frame("time"=stime, "event"=sevent, "score"=age)

#Cox's model
sbrier.score2proba(data.tr=dd, data.ts=dd, method="cox")
#product-limit estimator
sbrier.score2proba(data.tr=dd, data.ts=dd, method="prodlim")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ univar }% __ONLY ONE__ keyword per line
