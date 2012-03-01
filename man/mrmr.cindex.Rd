\name{mrmr.cindex}
\alias{mrmr.cindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute the concordance index for survival or binary class prediction}
\description{
Function to compute the minimum redundancy - maximum relevance (mRMR) ranking for a risk prediction or a binary classification task based on the concordance index. The mRMR feature selection has been adapted to use the concordance index to estimate the correlation between a variable and the output (binary or survival) data. 
}
\usage{
mrmr.cindex(x, surv.time, surv.event, cl, weights, comppairs=10, strat, alpha = 0.05, outx = TRUE, method = c("conservative", "noether", "nam"), alternative = c("two.sided", "less", "greater"), na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of risk predictions.}
  \item{surv.time}{a vector of event times.}
  \item{surv.event}{a vector of event occurence indicators.}
  \item{cl}{a vector of binary class indicators.}
  \item{weights}{weight of each sample.}
  \item{comppairs}{threshold for compairable patients.}
  \item{strat}{stratification indicator.}
  \item{alpha}{apha level to compute confidence interval.}
  \item{outx}{set to \code{TRUE} to not count pairs of observations tied on \code{x} as a relevant pair. This results in a Goodman-Kruskal gamma type rank correlation.}
  \item{method}{can take the value \code{conservative}, \code{noether} or \code{name} (see paper Pencina et al. for details).}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} (concordance index is greater than 0.5) or \code{"less"} (concordance index is less than 0.5). You can specify just the initial letter.}
  \item{na.rm}{\code{TRUE} if missing values should be removed.}
}
%%\details{
%%}
\note{
The "direction" of the concordance index (< 0.5 or > 0.5) is the opposite than the \link[Hmisc]{rcorr.cens} function in the \code{Hmisc} package. So you can easily get the same results than \link[Hmisc]{rcorr.cens} by  changing the sign of \code{x}.
}
\value{
  A mRMR ranking

}
\references{Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996) "Tutorial in biostatistics: multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing error", \emph{Statistics in Medicine}, \bold{15}, pages 361--387.

Pencina, M. J. and D'Agostino, R. B. (2004) "Overall C as a measure of discrimination in survival analysis: model specific population value and confidence interval estimation", \emph{Statistics in Medicine}, \bold{23}, pages 2109--2123, 2004.}
\author{Benjamin Haibe-Kains, Markus Schroeder}
%\note{ ~~further notes~~ }
\seealso{\code{\link[Hmisc]{rcorr.cens}}, \code{\link[CPE]{phcpe}}, \code{\link[clinfun]{coxphCPE}}}
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
sex <- sample(0:1, 100, replace=TRUE)
stime <- rexp(100)
cens   <- runif(100,.5,2)
sevent  <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
strat <- sample(1:3, 100, replace=TRUE)
weight <- runif(100, min=0, max=1)
comppairs <- 10
xx <- data.frame("age"=age, "sex"=sex)
cat("survival prediction:\n")
mrmr.cindex(x=xx, surv.time=stime, surv.event=sevent, strat=strat, weights=weight, method="noether", comppairs=comppairs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{univar}% __ONLY ONE__ keyword per line
