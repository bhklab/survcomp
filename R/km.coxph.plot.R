'km.coxph.plot' <-
function(formula.s, data.s, weight.s, sub.s="all", x.label, y.label, main.title, sub.title, leg.text, leg.pos="bottomright", leg.bty="o", leg.inset=0.05, o.text, v.line, h.line, .col=1:4, .lty=1, .lwd=1, show.n.risk=FALSE, n.risk.step, n.risk.cex=0.85, verbose=TRUE, ...) {
	
  if(length(sub.s) == 1 && sub.s=="all") { sub.s <- rep(TRUE, nrow(data.s)) }
  assign("sub.s", sub.s, envir=.GlobalEnv)
	if(missing(sub.title)) { sub.title <- NULL }
	if(missing(leg.text)) { leg.text <- NULL }
  if(missing(weight.s)) { weight.s <- array(1, dim=nrow(data.s), dimnames=list(rownames(data.s))) }
  assign("weight.s", weight.s, envir=.GlobalEnv)
	
	ng <- length(leg.text)
    old.mar <- par("mar")
    on.exit( par( mar = old.mar ) )
    .xaxt="s"
    .xlab=x.label
    if( show.n.risk ) {
        par(mar = old.mar + c(ng,8,3,0))
        .xaxt="n"
        .xlab = ""
    }

    plot(survfit(formula.s, data=data.s, weights=weight.s, subset=sub.s), xaxt=.xaxt, col=.col, lty=.lty, lwd=.lwd, xlab=.xlab, ylab=y.label, ... )
    title(main.title)
	
    if(!missing(v.line) && !is.null(v.line)) { abline(v=v.line, lty=3, col="purple") }
	if(!missing(h.line) && !is.null(h.line)) { abline(h=h.line, lty=3, col="purple") }

    if(!is.null(leg.text)) { legend(x=leg.pos, xjust=0, yjust=1, legend=leg.text, col=.col, lty=.lty, lwd=.lwd, cex=0.9, bg="white", inset=leg.inset, bty=leg.bty) }
    if(!is.null(sub.title)) { mtext(sub.title, line=-4, outer=TRUE) }
    if(missing(o.text) ) {
		sdf <- survdiff(formula.s, data=data.s, subset=sub.s)
	    if(verbose) { print(sdf) }
        p.val <- 1-pchisq(sdf$chisq,length(sdf$n)-1)
        #if( p.val < 0.001 ) o.text <- "P < 0.001" else o.text <- paste("P =", signif(p.val,3)) #, "(log-rank test)")
        o.text <- sprintf("P = %.1E", p.val)
    }
    if(is.null(o.text)) { o.text <- FALSE }
    text(0,0, o.text, cex=0.85, pos=4)

    if( show.n.risk ) {
        usr.xy <- par( "usr" )
        nrisk <- no.at.risk( formula.s, data.s, sub.s, n.risk.step, floor(usr.xy[2]) )
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

    if( exists("sub.s", envir=.GlobalEnv) ) remove("sub.s", envir=.GlobalEnv)
}