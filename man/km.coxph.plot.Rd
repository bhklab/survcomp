\name{km.coxph.plot}
\alias{km.coxph.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to plot several Kaplan-Meier survival curves }
\description{
Function to plot several Kaplan-Meier survival curves with number of individuals at risk at some time points.
}
\usage{
km.coxph.plot(formula.s, data.s, x.label, y.label, main.title, sub.title, leg.text, leg.pos = "bottomright", leg.bty = "o", leg.inset = 0.05, o.text, v.line, h.line, .col = 1:4, .lty = 1, .lwd = 1, show.n.risk = FALSE, n.risk.step, n.risk.cex = 0.85, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula.s}{formula composed of a \code{Surv} object and a strata variable (i.e. stratification).}
  \item{data.s}{data frame composed of the variables used in the formula.}
  \item{x.label}{label for the y-axis.}
  \item{y.label}{label for the x-axis.}
  \item{main.title}{main title at the top of the plot.}
  \item{sub.title}{subtitle at the bottom of the plot.}
  \item{leg.text}{text in the legend.}
  \item{leg.pos}{the location may also be specified by setting 'x' to a single keyword from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}. This places the legend on the inside of the plot frame at the given location.}
  \item{leg.bty}{the type of box to be drawn around the legend. The allowed values are \code{"o"} (the default) and \code{"n"}.}
  \item{leg.inset}{inset distance from the margins as a fraction of the plot region. Default value is 0.05.}
  \item{o.text}{ plot the logrank p-value.}
  \item{v.line}{x coordinate(s) for vertical line(s).}
  \item{h.line}{y coordinate(s) for horizontal line(s).}
  \item{.col}{vector of colors for the different survival curves.}
  \item{.lty}{vector of line types for the different survival curves}
  \item{.lwd}{vector of line widths for the different survival curves.}
  \item{show.n.risk}{if \code{TRUE}, show the numbers of samples at risk for each time step.}
  \item{n.risk.step}{vector specifying the time to be the steps for displaying the number of individuals at risk.}
  \item{n.risk.cex}{size of the number of individuals at risk. Default value is 0.85.}
  \item{verbose}{ verbosity level (\code{TRUE} or \code{FALSE}). Default value is \code{TRUE}.}
  \item{\dots}{additional parameters to be passed to the \code{\link[graphics]{plot}} function.}
}
\details{
  The original version of this function was kindly provided by Dr Christos Hatzis (January, 17th 2006).
}
%\value{}
%\references{}
\author{Christos Hatzis, Benjamin Haibe-Kains}
%\note{ ~~further notes~~ }
\seealso{\code{\link[survival]{survfit}}, \code{\link[survival]{coxph}}}
\examples{
set.seed(12345)
stime <- rexp(100) * 10
cens   <- runif(100,.5,2) * 10
sevent  <- as.numeric(stime <= cens)
stime <- pmin(stime, cens)
strat <- sample(1:3, 100, replace=TRUE)
dd <- data.frame("surv.time"=stime, "surv.event"=sevent, "strat"=strat)

km.coxph.plot(formula.s=Surv(surv.time, surv.event) ~ strat, data.s=dd,
  sub.s="all", x.label="Time (years)", y.label="Probability of survival",
  main.title="", leg.text=paste(c("Low", "Intermediate", "High"), "   ", sep=""),
  leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen", "darkred"),
  .lty=c(1,1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.85, verbose=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
\keyword{univar}% __ONLY ONE__ keyword per line
