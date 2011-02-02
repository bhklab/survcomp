\name{test.hetero.test}
\alias{test.hetero.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to test the heterogeneity of set of probabilities }
\description{
  The function tests whether a set of p-values are heterogeneous.
}
\usage{
test.hetero.test(p, weight, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of p-values }
  \item{weight}{ vector of weights (e.g. sample size of each study) }
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise }
}
\details{ The p-values should be one-sided and computed from the same null hypothesis. }
\value{
  \item{Q}{ Q statistic }
  \item{p.value}{ p-value of the heterogeneity test}
}
\references{ Cochrane, W. G. (1954) "The combination of estimates from different experiments", \emph{Biometrics}, \bold{10}, pages 101--129.

Whitlock, M. C. (2005) "Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach", \emph{J. Evol. Biol.}, \bold{18}, pages 1368--1373. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{combine.test} }
\examples{
p <- c(0.01, 0.13, 0.07, 0.2)
w <- c(100, 50, 200, 30)

#with equal weights
test.hetero.test(p=p)
#with p-values weighted by the sample size of the studies
test.hetero.test(p=p, weight=w)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
