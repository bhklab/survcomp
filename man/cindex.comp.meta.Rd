\name{cindex.comp.meta}
\alias{cindex.comp.meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compare two concordance indices}
\description{
This function compares two lists of concordance indices computed from the same survival data by using the function \code{concordance.index}. The statistical test is a Student t test for dependent samples.
}
\usage{
cindex.comp.meta(list.cindex1, list.cindex2, hetero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list.cindex1}{first list of concordance indices as returned by the \code{concordance.index} function.}
  \item{list.cindex2}{second list of concordance indices as returned by the \code{concordance.index} function.}
  \item{hetero}{if TRUE, a random effect model is use to compute the meta-estimators. Otherwise a fixed effect model is used.}
}
\details{
  In meta-analysis, we estimate the statistic of interest in several independent datasets. It results a list of estimates such as list of concordance indices. The two lists of concordance indices must be computed from the same samples (and corresponding survival data). The function computes a meta-estimator for the correlations between the two scores and uses a Student t test for dependent samples.
}
%%\note{
%%Since the comparison between the two concordance indices is done with a paired t-test that relies on the computation of the Pearson correlation coefficient between the original predictors \code{x}, any non-linear transformation of these predictors will affect the computation of the p-value.
%%We are now using spearman correlation to circumvent this difficulty.
%%}
\value{
  \item{p.value }{p-value from the Student t test for the comparison cindex1 > cindex2.}
  \item{cindex1 }{meta-estimator of the first concordance index.}
  \item{cindex2 }{meta-estimator of the second concordance index.}
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
c1.1 <- concordance.index(x=age, surv.time=stime, surv.event=sevent, method="noether")
c2.1 <- concordance.index(x=size, surv.time=stime, surv.event=sevent, method="noether")
#second dataset
set.seed(54321)
age <- rnorm(110, 53, 10)
size <- rexp(110,1.1)
stime <- rexp(110)
cens <- runif(110,.55,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
c1.2 <- concordance.index(x=age, surv.time=stime, surv.event=sevent, method="noether")
c2.2 <- concordance.index(x=size, surv.time=stime, surv.event=sevent, method="noether")
cindex.comp.meta(list.cindex1=list("cindex.age1"=c1.1, "cindex.age2"=c1.2), list.cindex2=list("cindex.size1"=c2.1, "cindex.size2"=c2.2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{htest}% __ONLY ONE__ keyword per line
