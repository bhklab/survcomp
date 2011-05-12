`hazard.ratio` <-
function(x, surv.time, surv.event, weights, strat, alpha=0.05, na.rm=FALSE, ...) {
	require(survival)
	if(!missing(weights)) {
		if(length(weights) != length(x)) { stop("bad length for parameter weights!") }
	} else { weights <- rep(1,  length(x)) }
	if(!missing(strat)) {
		if(length(strat) != length(x)) { stop("bad length for parameter strat!") }
	} else { strat <- rep(1,  length(x)) }
	cc.ix <- complete.cases(x, surv.time, surv.event, weights, strat)
	if(any(!cc.ix) & !na.rm) { stop("NA values are present!") }
	sx <- x[cc.ix]
	oo <- order(sx, decreasing=FALSE)
	sx <- sx[oo]
	stime <- surv.time[cc.ix][oo]
	sevent <- surv.event[cc.ix][oo]
	sweights <- weights[cc.ix][oo]
	sstrat <- strat[cc.ix][oo]	
	data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event)
	#fit the cox model		
	options(warn=2)
	rr <- try(coxph(Surv(stime, sevent) ~ strata(sstrat) + sx, weights=sweights, ...))
	options(warn=0)
	if(class(rr) == "try-error") {
		res <- list("hazard.ratio"=NA, "coef"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=sum(cc.ix), "coxm"=NA, "data"=data)
	} else {
		hrcoef <- rr$coefficients
		hrse <- sqrt(drop(rr$var))
		names(hrcoef) <- names(hrse) <- NULL
		mypp <- pchisq(2 * (rr$loglik[2] - rr$loglik[1]), df=1, lower.tail=FALSE)
		## mypp <- pchisq((hrcoef / hrse)^2, df=1, lower.tail=FALSE)
		res <- list("hazard.ratio"=exp(hrcoef), "coef"=hrcoef, "se"=hrse, "lower"=exp(hrcoef - qnorm(alpha / 2, lower.tail=FALSE) * hrse), "upper"=exp(hrcoef + qnorm(alpha / 2, lower.tail=FALSE) * hrse), "p.value"=mypp, "n"=rr$n, "coxm"=rr, "data"=data)
	}
	
	return(res)
}