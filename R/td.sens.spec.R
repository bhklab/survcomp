td.sens.spec <-
function(cl, surv.time, surv.event, time, span=0, sampling=FALSE, na.rm=FALSE, ...) {
	#require(survivalROC)
	
	if((length(cl) + length(surv.time) + length(surv.event)) != (3 * length(cl))) { stop("paramaters cl, surv.time and surv.event must have the same length!") }
	if(is.null(names(cl))) { names(cl) <- names(surv.time) <- names(surv.event) <- paste("X", 1:length(cl), sep=".") }
	cc.ix <- complete.cases(cl, surv.time, surv.event)
	if(all(!cc.ix) & !na.rm) { stop("missing values are present!") }
	cl2 <- cl[cc.ix]
	ucl2 <- sort(unique(cl2))
	if(length(ucl2) != 2) { stop("cl must be binary!") }
	oo <- order(cl2, decreasing=FALSE)
	cl2 <- cl2[oo]
	st <- surv.time[cc.ix][oo]
	se <- surv.event[cc.ix][oo]
	marker.fake <- 1:length(cl2)
	names(marker.fake) <- names(cl2)
	mycutoff <- sum(cl2 == ucl2[1])
	
	##using the survival.C function
	rr <- survivalROC::survivalROC.C(Stime=st, status=se, marker=marker.fake, predict.time=time, span=span,  ...)
	#rr <- survivalROC::survivalROC(Stime=st, status=se, marker=marker.fake, cut.values=mycutoff, predict.time=time, span=span, lambda=lambda, ...)

	sens.se <- spec.se <- NA
	if(sampling) {
		#require(bootstrap)
		
		theta.foo1 <- function(x, cl, surv.time, surv.event, time, ...) {
			cl <- cl[x]
			oo <- order(cl, decreasing=FALSE)
			cl <- cl[oo]
			surv.time <- surv.time[x][oo]
			surv.event <- surv.event[x][oo]
			ucl <- sort(unique(cl))
			marker.fake <- 1:length(cl)
			names(marker.fake) <- names(cl)
			mycutoff <- sum(cl == ucl[1])
			rr <- survivalROC::survivalROC.C(Stime=surv.time, status=surv.event, marker=marker.fake,  predict.time=time, span=span, ...)
			#rr <- survivalROC::survivalROC(Stime=surv.time, status=surv.event, marker=marker.fake, cut.values=mycutoff, predict.time=time, span=span, lambda=lambda, ...)
			return(list("sens"=rr$TP[which(rr$cut.values == mycutoff)], "spec"=1-rr$FP[which(rr$cut.values == mycutoff)]))
		}
		myx <- 1:length(cl2)
		myx2 <- myx
		sens.values <- spec.values <- rep(NA, length(cl2))
		names(sens.values) <- names(spec.values) <- names(cl2)
		for(i in length(myx2):1) {
			tt <- theta.foo1(x=myx[-myx2[i]], cl=cl2, surv.time=st, surv.event=se, time=time, ...)
			sens.values[is.na(sens.values) & myx >= myx2[i]] <- tt$sens
			spec.values[is.na(spec.values) & myx >= myx2[i]] <- tt$spec
		}
		
		theta.foo2 <- function(x, stat=c("sens", "spec"), xdata) {
			stat <- match.arg(stat)
			xdata <- unlist(xdata[stat])
			return(xdata[-x])	
		}
		
		sens.se <- bootstrap::jackknife(x=1:length(cl2), theta=theta.foo2, stat="sens", xdata=list("sens"=sens.values, "spec"=spec.values))$jack.se
		spec.se <- bootstrap::jackknife(x=1:length(cl2), theta=theta.foo2, stat="spec", xdata=list("sens"=sens.values, "spec"=spec.values))$jack.se	
	}
	
	return(list("sens"=rr$TP[which(rr$cut.values == mycutoff)], "sens.se"=sens.se, "spec"=1-rr$FP[which(rr$cut.values == mycutoff)], "spec.se"=spec.se))	
}