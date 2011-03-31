\name{test.hetero.est}
\alias{test.hetero.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to test the heterogeneity of set of probabilities }
\description{
  The function tests whether a set of p-values are heterogeneous.
}
\usage{
test.hetero.est(x, x.se, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of estimates }
  \item{x.se}{ vector of standard errors of the corresponding estimates }
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise }
}
\details{ The heterogeneity test is known to be very conservative. Consider a p-value < 0.1 as significant. }
\value{
  \item{Q}{ Q statistic }
  \item{p.value}{ p-value of the heterogeneity test}
}
\references{ Cochrane, W. G. (1954) "The combination of estimates from different experiments", \emph{Biometrics}, \bold{10}, pages 101--129. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{combine.test} }
\examples{
set.seed(12345)
x1 <- rnorm(100, 50, 10) + rnorm(100, 0, 2)
m1 <- mean(x1)
se1 <- sqrt(var(x1))
x2 <- rnorm(100, 75, 15) + rnorm(100, 0, 5)
m2 <- mean(x2)
se2 <- sqrt(var(x2))

test.hetero.est(x=c(m1, m2), x.se=c(se1, se2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
