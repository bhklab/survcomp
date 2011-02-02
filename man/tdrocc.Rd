\name{tdrocc}
\alias{tdrocc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute time-dependent ROC curves }
\description{
  The function is a wrapper for the \code{\link[survivalROC]{survivalROC}} function in order to compute the time-dependent ROC curves.
}
\usage{
tdrocc(x, surv.time, surv.event, surv.entry = NULL, time, cutpts = NA, na.rm = FALSE, verbose = FALSE, span = 0, lambda = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of risk scores. }
  \item{surv.time}{ vector of times to event occurrence. }
  \item{surv.event}{ vector of event occurrence indicators. }
  \item{surv.entry}{ entry time for the subjects. }
  \item{time}{ time point for the ROC curve. }
  \item{cutpts}{ cut points for the risk score. }
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise. }
  \item{verbose}{ verbosity of the function. }
  \item{span}{ Span for the NNE, need either lambda or span for NNE. }
  \item{lambda}{ smoothing parameter for NNE. }
  \item{\dots}{ additional arguments to be passed to the \code{\link[survivalROC]{survivalROC}} function.  }
}
%\details{}
\value{
  \item{spec}{specificity estimates}
  \item{sens}{sensitivity estimates}
  \item{rule}{rule to compute the predictions at each cutoff}
  \item{cuts}{cutoffs}
  \item{time}{time point at which the time-dependent ROC is computed}
  \item{survival}{overall survival at the time point}
  \item{AUC}{Area Under the Curve (AUC) of teh time-dependent ROC curve}
  \item{data}{survival data and risk score used to compute the time-dependent ROC curve}
}
\references{ Heagerty, P. J. and Lumley, T. L. and Pepe, M. S. (2000) "Time-Dependent ROC Curves for Censored Survival Data and a Diagnostic Marker", \emph{Biometrics}, \bold{56}, pages 337--344. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link[survivalROC]{survivalROC}} }
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
stime <- rexp(100)
cens <- runif(100,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
tdroc <- tdrocc(x=age, surv.time=stime, surv.event=sevent, time=1, na.rm=TRUE, verbose=FALSE)
##plot the time-dependent ROC curve
plot(x=1-tdroc$spec, y=tdroc$sens, type="l", xlab="1 - specificity", ylab="sensitivity", xlim=c(0, 1), ylim=c(0, 1))
lines(x=c(0,1), y=c(0,1), lty=3, col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
