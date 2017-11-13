\name{ibsc.comp}
\alias{ibsc.comp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compare two IBSCs }
\description{
  This function compares two integrated Briers scores (IBSC) through the estimation of the Brier scores (BSC) at some points in time. The statistical test is a Wilcoxon rank sum test for dependent samples.
}
\usage{
ibsc.comp(bsc1, bsc2, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bsc1}{ vector of BSCs computed from the first predicted probabilities for some points in time }
  \item{bsc2}{ vector of BSCs computed from the second predicted probabilities for some points in time }
  \item{time}{ vector of points in time for which the BSCs are computed }
}
\details{
  The two vectors of BSCs must be computed from the same samples (and corresponding survival data) and for the same points in time. The function uses a Wilcoxon rank sum test for dependent samples.
}
\value{
  \item{p.value }{p-value from the Wilcoxon rank sum test for the comparison ibsc1 < ibsc2}
  \item{ibsc1 }{value of the IBSC for the first set of BSCs}
  \item{ibsc2 }{value of the IBSC for the second set of BSCs}
}
\references{ Wilcoxon, F. (1945) "Individual comparisons by ranking methods", \emph{Biometrics Bulletin}, \bold{1}, pages 80--83.
	
Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A comparative study of survival models for breast cancer prognostication based on microarray  data: does a single gene beat them all?", \emph{Bioinformatics}, \bold{24}, 19, pages 2200--2208. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link{sbrier.score2proba}}, \code{\link[ipred]{sbrier}} }
\examples{
set.seed(12345)
age <- rnorm(30, 50, 10)
size <- rexp(30,1)
stime <- rexp(30)
cens <- runif(30,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
##Brier scores
##size
dd <- data.frame("time"=stime, "event"=sevent, "score"=size)
bsc1 <- sbrier.score2proba(data.tr=dd, data.ts=dd, method="cox")
##size
dd <- data.frame("time"=stime, "event"=sevent, "score"=age)
bsc2 <- sbrier.score2proba(data.tr=dd, data.ts=dd, method="cox")
if(!all(bsc1$time == bsc2$time)) {
  stop("the two vector of BSCs must be computed for the same points in time!") }
ibsc.comp(bsc1=bsc1$bsc, bsc2=bsc2$bsc, time=bsc1$time)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ htest }% __ONLY ONE__ keyword per line
