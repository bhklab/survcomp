\name{hr.comp.meta}
\alias{hr.comp.meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compare two concordance indices}
\description{
This function compares two lists of hazard ratios computed from the same survival data by using the function \code{hazard.ratio}. The statistical test is a Student t test for dependent samples.
}
\usage{
hr.comp.meta(list.hr1, list.hr2, hetero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list.hr1}{first list of D indices as returned by the \code{hazard.ratio} function.}
  \item{list.hr2}{second list of D indices as returned by the \code{hazard.ratio} function.}
  \item{hetero}{if TRUE, a random effect model is use to compute the meta-estimators. Otherwise a fixed effect model is used.}
}
\details{
  In meta-analysis, we estimate the statistic of interest in several independent datasets. It results a list of estimates such as list of hazard ratios. The two lists of hazrd ratios must be computed from the same samples (and corresponding survival data). The function computes a meta-estimator for the correlations between the two scores and uses a Student t test for dependent samples.
}
\value{
  \item{p.value }{p-value from the Student t test for the comparison hr1 > hr2.}
  \item{hr1 }{meta-estimator of the first D index.}
  \item{hr2 }{meta-estimator of the second D index.}
}
\references{Cochrane, W. G. (1954) "The combination of estimates from different experiments", \emph{Biometrics}, \bold{10}, pages 101--129.

Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A comparative study of survival models for breast cancer prognostication based on microarray  data: does a single gene beat them all?", \emph{Bioinformatics}, \bold{24}, 19, pages 2200--2208.}
\author{ Benjamin Haibe-Kains}
%\note{ ~~further notes~~ }
\seealso{\code{\link{concordance.index}}.}
\examples{
#first dataset
set.seed(12345)
age <- rnorm(100, 50, 10)
size <- rexp(100,1)
stime <- rexp(100)
cens <- runif(100,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
h1.1 <- hazard.ratio(x=age, surv.time=stime, surv.event=sevent)
h2.1 <- hazard.ratio(x=size, surv.time=stime, surv.event=sevent)
#second dataset
set.seed(54321)
age <- rnorm(110, 53, 10)
size <- rexp(110,1.1)
stime <- rexp(110)
cens <- runif(110,.55,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
h1.2 <- hazard.ratio(x=age, surv.time=stime, surv.event=sevent)
h2.2 <- hazard.ratio(x=size, surv.time=stime, surv.event=sevent)
hr.comp.meta(list.hr1=list("hr.age1"=h1.1, "hr.age2"=h1.2), list.hr2=list("hr.size1"=h2.1, "hr.size2"=h2.2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{htest}% __ONLY ONE__ keyword per line
