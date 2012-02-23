\name{metaplot.surv}
\alias{metaplot.surv}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Meta-analysis plot (forest plot)}
\description{
  Plot confidence intervals with boxes indicating the sample
  size/precision and optionally a diamond indicating a summary
  confidence interval. This function is usually called by \code{plot}
  methods for meta-analysis objects. Additional, you can specifiy your
  own lower and upper boarder from the confidence interval.
}
\usage{
metaplot.surv(mn, se=NULL, lower=NULL, upper=NULL, nn=NULL,
    labels=NULL, conf.level = .95, xlab = "", ylab = "", xlim = NULL,
    summn = NULL, sumse = NULL, sumlower = NULL, sumupper = NULL,
    sumnn = NULL, summlabel = "Summary", logeffect = FALSE,
    lwd = 2, boxsize = 1, zero = as.numeric(logeffect),
    colors, xaxt="s", logticks=TRUE, ... )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mn}{point estimates from studies}
  \item{se}{standard errors of \code{mn}}
  \item{lower}{Vector of lower ends of confidence intervals}
  \item{upper}{Vector of upper ends of confidence intervals}
  \item{nn}{precision: box ares is proportional to this. \code{1/se^2}
    is the default}
  \item{labels}{labels for each interval}
  \item{conf.level}{Confidence level for confidence intervals}
  \item{xlab}{label for the point estimate axis}
  \item{ylab}{label for the axis indexing the different studies}
  \item{xlim}{the range for the x axis.}
  \item{summn}{summary estimate}
  \item{sumse}{standard error of summary estimate}
  \item{sumlower}{lower end of confidence intervals of summary estimate}
  \item{sumupper}{upper end of confidence intervals of summary estimate}
  \item{sumnn}{precision of summary estimate}
  \item{summlabel}{label for summary estimate }
  \item{logeffect}{\code{TRUE} to display on a log scale }
  \item{lwd}{line width}
  \item{boxsize}{Scale factor for box size}
  \item{zero}{"Null" effect value}
  \item{xaxt}{use \code{"n"} for no x-axis (to add a customised one)}
  \item{logticks}{if \code{TRUE} and \code{logscale}, have tick values
    approximately equally spaced on a log scale}.
  \item{colors}{see \code{\link{meta.colors}}}
  \item{\dots}{Other graphical parameters }
}

\value{
  This function is used for its side-effect.
}

\seealso{
  \code{\link{forestplot.surv}} for more flexible plots
}

\examples{
metaplot.surv(mn=c(0.4,0.5,0.6), lower=c(0.35,0.4,0.57), upper=c(0.45,0.6,0.63), labels=c("A","B","C"), xlim=c(0.2,0.8), boxsize=0.5, zero=0.5, col=rmeta::meta.colors(box="royalblue",line="darkblue",zero="firebrick"))
}
\keyword{ hplot }

