\name{forestplot.surv}
\alias{forestplot.surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Forest plots enables to display performance estimates of survival models}
\description{
Draw a forest plot together with a table of text.
}
\usage{
forestplot.surv(labeltext, mean, lower, upper, align = NULL,
  is.summary = FALSE, clip = c(-Inf, Inf), xlab = "", zero = 0,
  graphwidth = unit(2, "inches"), col = meta.colors(), xlog = FALSE,
  box.size = NULL, x.ticks = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{labeltext}{Matrix of strings or \code{NA}s for blank spaces }
  \item{mean}{Vector of centers of confidence intervals (or \code{NA}s for blank space) }
  \item{lower}{Vector of lower ends of confidence intervals }
  \item{upper}{Vector of upper ends of confidence intervals  }
  \item{align}{Vector giving alignment (\code{l},\code{r},\code{c}) for columns of table }
  \item{is.summary}{Vector of logicals. Summary lines have bold text and diamond confidence intervals.}
  \item{clip}{Lower and upper limits for clipping confidence intervals to arrows }
  \item{xlab}{x-axis label }
  \item{zero}{x-axis coordinate for zero line }
  \item{graphwidth}{Width of confidence interval graph }
  \item{col}{See \code{\link[rmeta]{meta.colors}}}
  \item{xlog}{If \code{TRUE}, x-axis tick marks are exponentiated}
  \item{box.size}{Override the default box size based on precision}
  \item{x.ticks}{Optional user-specified x-axis tick marks. Specify \code{NULL} to use the defaults, \code{numeric(0)} to omit the x-axis.}
  \item{\dots}{Not used.}
}
\details{
This function is more flexible than \code{\link[rmeta]{metaplot}} and the
\code{plot} methods for meta-analysis objects, but requires more work by
the user.

In particular, it allows for a table of text, and clips confidence
intervals to arrows when they exceed specified limits.
}
\value{
None
}
\references{
rmeta package, CRAN, Thomas Lumley <tlumley@u.washington.edu>. Functions for simple fixed and random effects meta-analysis for two-sample comparisons and cumulative meta-analyses. Draws standard summary plots, funnel plots, and computes summaries and tests for association and heterogeneity.
}
\seealso{\code{metaplot}, \code{forestplot}}
\examples{
##
}
\keyword{hplot}% at least one, from doc/KEYWORDS


