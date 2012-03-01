`hazard.ratio` <-
function(x, surv.time, surv.event, weights, strat, alpha=0.05, method.test=c("logrank", "likelihood.ratio", "wald"), na.rm=FALSE, ...) {
	require(survival)
	method.test <- match.arg(method.test)
	if(!missing(weights)) {
		if(length(weights) != length(x)) { stop("bad length for parameter weights!") }
	} else { weights <- rep(1,  length(x)) }
	if(!missing(strat)) {
		if(length(strat) != length(x)) { stop("bad length for parameter strat!") }
	} else { strat <- rep(1,  length(x)) }
	cc.ix <- complete.cases(x, surv.time, surv.event, weights, strat)
	if(sum(cc.ix) < 3) {
	## not enough observations	
		data <- list("x"=x, "z"=rep(NA, length(x)), "surv.time"=surv.time, "surv.event"=surv.event, "weights"=weights, "strat"=strat)
		return(list("hazard.ratio"=NA, "coef"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=sum(cc.ix), "coxm"=NA, "data"=data))	
	}
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
		mystat <- NA
		switch(method.test, 
		"logrank"={
			mystat <- rr$score
		},
		"likelihood.ratio"={
			mystat <- 2 * (rr$loglik[2] - rr$loglik[1])
		},
		"wald"={
			mystat <- rr$wald.test
			##(hrcoef / hrse)^2
		}) 
		mypp <- pchisq(mystat, df=1, lower.tail=FALSE)
		res <- list("hazard.ratio"=exp(hrcoef), "coef"=hrcoef, "se"=hrse, "lower"=exp(hrcoef - qnorm(alpha / 2, lower.tail=FALSE) * hrse), "upper"=exp(hrcoef + qnorm(alpha / 2, lower.tail=FALSE) * hrse), "p.value"=mypp, "n"=rr$n, "coxm"=rr, "data"=data)
	}
	
	return(res)
}