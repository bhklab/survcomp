\name{fisherz}
\alias{fisherz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute Fisher z transformation }
\description{
  The function computes the Fisher z transformation useful to calculate the confidence interval of Pearson's correlation coefficient.
}
\usage{
fisherz(x, inv = FALSE, eps = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ value, e.g. Pearson's correlation coefficient }
  \item{inv}{ \code{TRUE} for inverse Fisher z transformation, \code{FALSE} otherwise }
  \item{eps}{ tolerance for extreme cases, i.e. \deqn{latex}{|x| \approx 1} when inv = \code{FALSE} and \deqn{latex}{|x| \approx Inf} when inv = \code{TRUE} }
}
\details{
The sampling distribution of Pearson's \eqn{latex}{\rho} is not normally distributed. R. A. Fisher developed a transformation now called \dQuote{Fisher's z transformation} that converts Pearson's \eqn{latex}{\rho} to the normally distributed variable z. The formula for the transformation is
\deqn{latex}{z = 1 / 2 [ \log(1 + \rho) - \log(1 - \rho) ]}

Two attributes of the distribution of the z statistic: (1) It is normally distributed and (2) it has a known standard error of
\deqn{latex}{\sigma_z = 1 / \sqrt{N - 3}}
where \eqn{latex}{N} is the number of samples.

Fisher's z is used for computing confidence intervals on Pearson's correlation and for confidence intervals on the difference between correlations.
}
\value{
Fisher's z statistic
}
\references{ R. A. Fisher (1915) "Frequency distribution of the values of the correlation coefficient in samples of an indefinitely large population". \emph{Biometrika}, \bold{10},pages 507--521. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link[stats]{cor}} }
\examples{
set.seed(12345)
x1 <- rnorm(100, 50, 10)
x2 <- runif(100,.5,2)
cc <- cor(x1, x2)
z <- fisherz(x=cc, inv=FALSE)
z.se <- 1 / sqrt(100 - 3)
fisherz(z, inv=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
