'no.at.risk' <-
function(formula.s, data.s, sub.s="all", t.step, t.end) {

   require(survival)

   if( length(sub.s)==1 && sub.s=="all" ) sub.s <- rep(TRUE, nrow(data.s))
   assign("sub.s", sub.s, envir=.GlobalEnv)

   sf <- survfit( formula.s, data=data.s, subset=sub.s )
   if(is.null(sf$strata)) {
		n.strata <- 1
		nnn <- "all"
		tms <- list("all"=sf$time)
   		rsk <- list("all"=sf$n.risk)
	} else {
		n.strata <- length(sf$strata)
		nnn <- names(sf$strata)
		#nnnt <- sf$ntimes.strata
		nnnt <- sf$strata
   		tms <- split(sf$time, rep(1:n.strata, times=nnnt))
   		rsk <- split(sf$n.risk, rep(1:n.strata, times=nnnt))
	}

   find.idx <- function( x, y ) {
		myx <- x > y
		if(all(!myx)) { ix <- which.max(x) } else { ix <- match(min(x[myx]), x ) - 1 }
		if( ix == 0 ) ix <- 1
		return(ix)
   }

   t.pts <- seq(0, t.end, t.step)
   nar <- data.frame("risk class"=nnn)
   for( i in 1:length(t.pts) ) {
       idx <- sapply(tms, FUN=find.idx, y=t.pts[i])
       if(is.null(sf$strata)) { n.r <- rsk[[c(1, idx[1])]] } else { n.r <- sapply( 1:n.strata, FUN=function(j) rsk[[c(j,idx[j])]] ) }
       nar <- cbind(nar, n.r)
       names(nar)[i+1] <- as.character(t.pts[i])
   }

   remove("sub.s", envir=.GlobalEnv)
   return(nar)
}
