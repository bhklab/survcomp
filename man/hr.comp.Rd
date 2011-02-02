\name{hr.comp}
\alias{hr.comp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to statistically compare two hazard ratios }
\description{
 This function compares two hazard ratios from their betas and standard errors as computed by a Cox model for instance. The statistical test is a Student t test for dependent samples. The two hazard ratios must be computed from the same survival data.
}
\usage{
hr.comp(hr1, hr2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hr1}{ first hazard ratio. }
  \item{hr2}{ second hazard ratio. }
}
\details{
  The two hazard ratios must be computed from the same samples (and corresponding survival data). The function uses a Student t test for dependent samples.
}
\value{
  \item{p.value }{p-value from the Student t test for the comparison beta1 > beta2 (equivalently hr1 > hr2)}
  \item{hr1 }{value of the first hazard ratio}
  \item{hr2 }{value of the second hazard ratio}
}
\references{ Student 1908) "The Probable Error of a Mean", \emph{Biometrika}, \bold{6}, 1, pages 1--25.
	
	Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A comparative study of survival models for breast cancer prognostication based on microarray  data: does a single gene beat them all?", \emph{Bioinformatics}, \bold{24}, 19, pages 2200--2208. }
\author{ Benjamin Haibe-Kains }
%\note{}
\seealso{ \code{\link[survival]{coxph}}, \code{\link{t.test}} }
\examples{
require(survival)
set.seed(12345)
age <- as.numeric(rnorm(100, 50, 10) >= 50)
size <- as.numeric(rexp(100,1) > 1)
stime <- rexp(100)
cens <- runif(100,.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
hr1 <- hazard.ratio(x=age, surv.time=stime, surv.event=sevent)
hr2 <- hazard.ratio(x=size, surv.time=stime, surv.event=sevent)
hr.comp(hr1=hr1, hr2=hr2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
\keyword{ htest }% __ONLY ONE__ keyword per line
