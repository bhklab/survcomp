#' @title Function to plot several Kaplan-Meier survival curves
#'
#' @description
#' Function to plot several Kaplan-Meier survival curves with number of
#'   individuals at risk at some time points.
#'
#' @usage
#' km.coxph.plot(formula.s, data.s, weight.s, x.label, y.label, main.title, sub.title,
#' leg.text, leg.pos = "bottomright", leg.bty = "o", leg.inset = 0.05, o.text, v.line,
#' h.line, .col = 1:4, .lty = 1, .lwd = 1, show.n.risk = FALSE, n.risk.step,
#' n.risk.cex = 0.85, verbose = TRUE, ...)
#'
#' @param formula.s formula composed of a Surv object and a strata variable
#'   (i.e. stratification).
#' @param data.s data frame composed of the variables used in the formula.
#' @param weight.s vector of weights of length nrow(data.s).
#' @param x.label label for the y-axis.
#' @param y.label label for the x-axis.
#' @param main.title main title at the top of the plot.
#' @param sub.title subtitle at the bottom of the plot.
#' @param leg.text text in the legend.
#' @param leg.pos the location may also be specified by setting 'x' to a single
#'   keyword from the list "bottomright", "bottom", "bottomleft", "left",
#'   "topleft", "top", "topright", "right" and "center". This places the
#'   legend on the inside of the plot frame at the given location.
#' @param leg.bty the type of box to be drawn around the legend. The allowed
#'   values are "o" (the default) and "n".
#' @param leg.inset inset distance from the margins as a fraction of the plot
#'   region. Default value is 0.05.
#' @param o.text plot the logrank p-value.
#' @param v.line x coordinate(s) for vertical line(s).
#' @param h.line y coordinate(s) for horizontal line(s).
#' @param .col vector of colors for the different survival curves.
#' @param .lty vector of line types for the different survival curves
#' @param .lwd vector of line widths for the different survival curves.
#' @param show.n.risk if `TRUE`, show the numbers of samples at risk for each time
#'   step.
#' @param n.risk.step vector specifying the time to be the steps for displaying
#'   the number of individuals at risk.
#' @param n.risk.cex size of the number of individuals at risk. Default value is
#'   0.85.
#' @param verbose verbosity level (TRUE or FALSE). Default value is TRUE.
#' @param ... additional parameters to be passed to the plot function.
#'
#' @details
#' The original version of this function was kindly provided by Dr Christos
#'   Hatzis (January, 17th 2006).
#'
#' @return
#' Several Kaplan-Meier survival curves with number of individuals at risk at
#'   some time points.
#'
#' @seealso
#' [survival::survfit], [survival::coxph]
#'
#' @examples
#' set.seed(12345)
#' stime <- rexp(100) * 10
#' cens   <- runif(100,.5,2) * 10
#' sevent  <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' strat <- sample(1:3, 100, replace=TRUE)
#' dd <- data.frame("surv.time"=stime, "surv.event"=sevent, "strat"=strat)
#' ddweights <- array(1, dim=nrow(dd))
#'
#' if (interactive()) {
#' km.coxph.plot(formula.s=Surv(surv.time, surv.event) ~ strat, data.s=dd,
#'   weight.s=ddweights, x.label="Time (years)", y.label="Probability of survival",
#'   main.title="", leg.text=paste(c("Low", "Intermediate", "High"), "   ", sep=""),
#'   leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen", "darkred"),
#'   .lty=c(1,1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.85, verbose=FALSE)
#' }
#'
#' @md
#' @export
km.coxph.plot <-
function(formula.s, data.s, weight.s, x.label, y.label, main.title, sub.title, leg.text, leg.pos="bottomright", leg.bty="o", leg.inset=0.05, o.text, v.line, h.line, .col=1:4, .lty=1, .lwd=1, show.n.risk=FALSE, n.risk.step, n.risk.cex=0.85, verbose=TRUE, ...) {

	if (missing(sub.title)) { sub.title <- NULL }
	if (missing(leg.text)) { leg.text <- NULL }
  if (missing(weight.s)) { weight.s <- array(1, dim=nrow(data.s), dimnames=list(rownames(data.s))) }
  ## weights should be > 0
  data.s <- data.s[!is.na(weight.s) & weight.s > 0, , drop=FALSE]
  weight.s <- weight.s[!is.na(weight.s) & weight.s > 0]
  pos <- 1
  envir = as.environment(pos)
  assign("weight.s", weight.s, envir = envir)
  weighted <- length(sort(unique(weight.s))) > 1

	ng <- length(leg.text)
    old.mar <- par("mar")
    on.exit( par( mar = old.mar ) )
    .xaxt="s"
    .xlab=x.label
    if (show.n.risk) {
        par(mar = old.mar + c(ng,8,3,0))
        .xaxt="n"
        .xlab = ""
    }

    plot(survfit(formula.s, data=data.s, weights=weight.s), xaxt=.xaxt, col=.col, lty=.lty, lwd=.lwd, xlab=.xlab, ylab=y.label, ... )
    title(main.title)

    if (!missing(v.line) && !is.null(v.line)) { abline(v=v.line, lty=3, col="purple") }
    if (!missing(h.line) && !is.null(h.line)) { abline(h=h.line, lty=3, col="purple") }

    if (!is.null(leg.text)) { legend(x=leg.pos, xjust=0, yjust=1, legend=leg.text, col=.col, lty=.lty, lwd=.lwd, cex=0.9, bg="white", inset=leg.inset, bty=leg.bty) }
    if (!is.null(sub.title)) { mtext(sub.title, line=-4, outer=TRUE) }
    if (missing(o.text)) {
		  sdf <- summary(survival::coxph(formula.s, data=data.s, weights=weight.s))
	    if(verbose) { print(sdf) }
        p.val <- sdf$sctest["pvalue"]
        o.text <- sprintf("Logrank P = %.1E", p.val)
    }
    if (is.null(o.text)) { o.text <- FALSE }
    text(0,0, o.text, cex=0.85, pos=4)

    if (show.n.risk) {
        usr.xy <- par( "usr" )
        nrisk <- no.at.risk(formula.s=formula.s, data.s=data.s, sub.s="all", t.step=n.risk.step, t.end=floor(usr.xy[2]) )
        at.loc <- seq(0, usr.xy[2], n.risk.step)
        axis(1, at=at.loc)
        mtext(x.label, side=1, line=2)
        mtext("No. At Risk", side=1, line=3, at=-0.5*n.risk.step, adj=1, cex=n.risk.cex, font=2)
        #nrsk.lbs <- sapply( strsplit(levels(nrisk[,1]),"="), FUN=function(x) x[2] )
        #if( any(is.na(nrsk.lbs)) ) nrsk.lbs <- leg.text
        for( i in 1:nrow(nrisk) ) {
            mtext(leg.text[i], side=1, line=3+i, at=-0.5*n.risk.step, adj=1, cex=n.risk.cex)
            mtext(nrisk[i,-1], side=1, at=at.loc, line=3+i, adj=1, cex=n.risk.cex)
       }
    }

    if( exists("weight.s", envir=.GlobalEnv) ) remove("weight.s", envir=.GlobalEnv)
}
