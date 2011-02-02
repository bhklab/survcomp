\name{censor.time}
\alias{censor.time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to artificially censor survival data }
\description{
  The function censors the survival data at a specific point in time. This is is useful if you used datasets having different follow-up periods.
}
\usage{
censor.time(surv.time, surv.event, time.cens = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surv.time}{ vector of times to event occurrence }
  \item{surv.event}{ vector of indicators for event occurrence }
  \item{time.cens}{ point in time at which the survival data must be censored }
}
%\details{}
\value{
  \item{surv.time.cens }{vector of censored times to event occurrence}
  \item{surv.event.cens }{vector of censored indicators for event occurrence}
}
%\references{}
\author{ Benjamin Haibe-Kains }
%\note{}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
set.seed(12345)
stime <- rexp(30)
cens <- runif(30,0.5,2)
sevent <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
censor.time(surv.time=stime, surv.event=sevent, time.cens=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
