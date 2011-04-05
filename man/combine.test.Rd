\name{combine.test}
\alias{combine.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to combine probabilities }
\description{
  The function combines several p-value estimated from the same null hypothesis in different studies involving independent data.
}
\usage{
combine.test(p, weight, method = c("fisher", "z.transform", "logit"),
  hetero = FALSE, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of p-values }
  \item{weight}{ vector of weights (e.g. sample size of each study) }
  \item{method}{ \code{fisher} for the Fisher's combined probability test, \code{z.transform} for the Z-transformed test, \code{logit} for the weighted Z-method }
  \item{hetero}{ \code{TRUE} is the heterogeneity should be taken into account, \code{FALSE} otherwise }
  \item{na.rm}{ \code{TRUE} if the missing values should be removed from the data, \code{FALSE} otherwise }
}
\details{ The p-values must be one-sided and computed from the same null hypothesis. }
\value{
  p-value
}
\references{ Whitlock, M. C. (2005) "Combining probability from independent tests: the weighted Z-method is superior to Fisher's approach", \emph{J. Evol. Biol.}, \bold{18}, pages 1368--1373. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{test.hetero.test} }
\examples{
p <- c(0.01, 0.13, 0.07, 0.2)
w <- c(100, 50, 200, 30)

#with equal weights
combine.test(p=p, method="z.transform")
#with p-values weighted by the sample size of the studies
combine.test(p=p, weight=w, method="z.transform")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
