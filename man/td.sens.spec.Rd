\name{td.sens.spec}
\alias{td.sens.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute sensitivity and specificity for a binary classification of survival data }
\description{
  The function is a wrapper for the \code{\link[survivalROC]{survivalROC.C}} function in order  to compute sensitivity and specificity for a binary classification of survival data.
}
\usage{
td.sens.spec(cl, surv.time, surv.event, time, span = 0, sampling = FALSE,
  na.rm = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{ vector of binary classes. }
  \item{surv.time}{ vector of times to event occurrence. }
  \item{surv.event}{ vector of event occurrence indicators. }
  \item{time}{ time point for sensitivity and specificity estimations. }
  \item{span}{ Span for the NNE. Default value is 0. }
  \item{sampling}{ jackknife procedure to estimate the standard error of sensitivity and specificity estimations.}
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise. }
  \item{\dots}{ additional arguments to be passed to the \code{\link[survivalROC]{survivalROC}} function.  }
}
\details{Only NNE method is used to estimate sensitivity and specificity (see \code{\link[survivalROC]{survivalROC.C}}). The standard error for sensitivity and specificity is estimated through jackknife procedure (see \code{\link[bootstrap]{jackknife}}).}
\value{
  \item{sens}{sensitivity estimate}
  \item{sens.se}{standard error for sensitivity estimate}
  \item{spec}{specificity estimate}
  \item{spec.se}{standard error for specificity estimate}
}
\references{ Heagerty, P. J. and Lumley, T. L. and Pepe, M. S. (2000) "Time-Dependent ROC Curves for Censored Survival Data and a Diagnostic Marker", \emph{Biometrics}, \bold{56}, pages 337--344. 

Efron, B. and   Tibshirani, R. (1986). "The Bootstrap Method for standard errors, confidence intervals, and other measures of statistical accuracy", \emph{Statistical Science}, \bold{1} (1), pages 1--35.
}
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link[survivalROC]{survivalROC}} }
\examples{
set.seed(12345)
gender <- sample(c(0,1), 100, replace=TRUE)
stime <- rexp(100)
cens <- runif(100,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
mysenspec <- td.sens.spec(cl=gender, surv.time=stime, surv.event=sevent,
  time=1, span=0, na.rm=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
