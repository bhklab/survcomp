\name{cindex.comp}
\alias{cindex.comp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compare two concordance indices}
\description{
This function compares two concordance indices computed from the same  data by using the function \code{concordance.index}. The statistical test is a Student t test for dependent samples.
}
\usage{
cindex.comp(cindex1, cindex2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cindex1}{first concordance index as returned by the \code{concordance.index} function.}
  \item{cindex2}{second concordance index as returned by the \code{concordance.index} function.}
}
\details{
  The two concordance indices must be computed from the same samples (and corresponding survival data). The function uses a Student t test for dependent samples.
}
%%\note{
%%Since the comparison between the two concordance indices is done with a paired t-test that relies on the computation of the Pearson correlation coefficient between the original predictors \code{x}, any non-linear transformation of these predictors will affect the computation of the p-value.
%%We are now using spearman correlation to circumvent this difficulty.
%%}
\value{
  \item{p.value }{p-value from the Student t test for the comparison cindex1 > cindex2.}
  \item{cindex1 }{value of the first concordance index.}
  \item{cindex2 }{value of the second concordance index.}
}
\references{Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A comparative study of survival models for breast cancer prognostication based on microarray  data: does a single gene beat them all?", \emph{Bioinformatics}, \bold{24}(19), pages 2200--2208.}
\author{ Benjamin Haibe-Kains}
%\note{ ~~further notes~~ }
\seealso{\code{\link{concordance.index}}.}
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
size <- rexp(100,1)
stime <- rexp(100)
cens <- runif(100,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
c1 <- concordance.index(x=age, surv.time=stime, surv.event=sevent, method="noether")
c2 <- concordance.index(x=size, surv.time=stime, surv.event=sevent, method="noether")
cindex.comp(c1, c2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{htest}% __ONLY ONE__ keyword per line
