\name{combine.est}
\alias{combine.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to combine estimates }
\description{
  The function combines several estimators using meta-analytical formula to compute a meta-estimate.
}
\usage{
combine.est(x, x.se, hetero = FALSE, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of estimates }
  \item{x.se}{ vector of standard errors of the corresponding estimates }
  \item{hetero}{ \code{TRUE} is the heterogeneity should be taken into account (random effect model), \code{FALSE} otherwise (fixed effect model) }
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise }
}
%\details{}
\value{
  \item{estimate}{ meta-estimate }
  \item{se}{ standard error of the meta-estimate}
}
\references{ Cochrane, W. G. (1954) "The combination of estimates from different experiments", \emph{Biometrics}, \bold{10}, pages 101--129. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{test.hetero.est} }
\examples{
set.seed(12345)
x1 <- rnorm(100, 50, 10) + rnorm(100, 0, 2)
m1 <- mean(x1)
se1 <- sqrt(var(x1))
x2 <- rnorm(100, 75, 15) + rnorm(100, 0, 5)
m2 <- mean(x2)
se2 <- sqrt(var(x2))

#fixed effect model
combine.est(x=c(m1, m2), x.se=c(se1, se2), hetero=FALSE)
#random effect model
combine.est(x=c(m1, m2), x.se=c(se1, se2), hetero=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
