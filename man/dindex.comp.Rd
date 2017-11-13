\name{dindex.comp}
\alias{dindex.comp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compare two D indices  }
\description{
 This function compares two D indices from their betas and standard errors as computed by a Cox model for instance. The statistical test is a Student t test for dependent samples. The two D indices must be computed using the \code{\link{D.index}} function from the same survival data.
}
\usage{
dindex.comp(dindex1, dindex2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dindex1}{ first D index }
  \item{dindex2}{ second D index }
}
\details{
  The two D indices must be computed using the \code{\link{D.index}} function from the same samples (and corresponding survival data). The function uses a Student t test for dependent samples.
}
\value{
  \item{p.value }{p-value from the Wilcoxon rank sum test for the comparison dindex1 > dindex2}
  \item{dindex1 }{value of the first D index}
  \item{dindex2 }{value of the second D index}
}
\references{ Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A comparative study of survival models for breast cancer prognostication based on microarray  data: does a single gene beat them all?", \emph{Bioinformatics}, \bold{24}, 19, pages 2200--2208. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link{D.index}}, \code{\link{t.test}} }
\examples{
set.seed(12345)
age <- rnorm(100, 50, 10)
size <- rexp(100,1)
stime <- rexp(100)
cens <- runif(100,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
d1 <- D.index(x=age, surv.time=stime, surv.event=sevent)
d2 <- D.index(x=size, surv.time=stime, surv.event=sevent)
dindex.comp(d1, d2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ htest }% __ONLY ONE__ keyword per line
