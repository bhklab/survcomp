\name{iauc.comp}
\alias{iauc.comp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compare two IAUCs through time-dependent ROC curves }
\description{
  This function compares two integrated areas under the curves (IAUC) through the results of time-dependent ROC curves at some points in time. The statistical test is a Wilcoxon rank sum test for dependent samples.
}
\usage{
iauc.comp(auc1, auc2, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{auc1}{ vector of AUCs computed from the first time-dependent ROC curves for some points in time }
  \item{auc2}{ vector of AUCs computed from the second time-dependent ROC curves for some points in time }
  \item{time}{ vector of points in time for which the AUCs are computed }
}
\details{
  The two vectors of AUCs must be computed from the same samples (and corresponding survival data) and for the same points in time. The function uses a Wilcoxon rank sum test for dependent samples.
}
\value{
  \item{p.value }{p-value from the Wilcoxon rank sum test for the comparison iauc1 > iauc2}
  \item{iauc1 }{value of the IAUC for the first set of time-depdent ROC curves}
  \item{iauc2 }{value of the IAUC for the second set of time-depdent ROC curves}
}
\references{ Wilcoxon, F. (1945) "Individual comparisons by ranking methods", \emph{Biometrics Bulletin}, \bold{1}, pages 80--83.

Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A comparative study of survival models for breast cancer prognostication based on microarray  data: does a single gene beat them all?", \emph{Bioinformatics}, \bold{24}, 19, pages 2200--2208. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link{tdrocc}}, \code{\link{wilcox.test}} }
\examples{
set.seed(12345)
age <- rnorm(30, 50, 10)
size <- rexp(30,1)
stime <- rexp(30)
cens <- runif(30,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
##time-dependent ROC curves
tt <- unique(sort(stime[sevent == 1]))
##size
mytdroc1 <- NULL
for(i in 1:length(tt)) {
	rr <- tdrocc(x=size, surv.time=stime, surv.event=sevent, time=tt[i], na.rm=TRUE, verbose=FALSE)
	mytdroc1 <- c(mytdroc1, list(rr))
}
auc1 <- unlist(lapply(mytdroc1, function(x) { return(x$AUC) }))
##age
mytdroc2 <- NULL
for(i in 1:length(tt)) {
	rr <- tdrocc(x=age, surv.time=stime, surv.event=sevent, time=tt[i], na.rm=TRUE, verbose=FALSE)
	mytdroc2 <- c(mytdroc2, list(rr))
}
auc2 <- unlist(lapply(mytdroc2, function(x) { return(x$AUC) }))
iauc.comp(auc1=auc1, auc2=auc2, time=tt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ htest }% __ONLY ONE__ keyword per line
