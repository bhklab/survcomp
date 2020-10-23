#' @name D.index
#'
#' @title Function to compute the D index
#'
#' @description
#' Function to compute the D index for a risk prediction, i.e. an estimate of
#'   the log hazard ratio comparing two equal-sized prognostic groups. This is
#'   a natural measure of separation between two independent survival
#'   distributions under the proportional hazards assumption.
#'
#' @usage
#' D.index(x, surv.time, surv.event, weights, strat, alpha = 0.05,
#'   method.test = c("logrank", "likelihood.ratio", "wald"), na.rm = FALSE, ...)
#'
#' @param x a vector of risk predictions.
#' @param surv.time a vector of event times.
#' @param surv.event a vector of event occurrence indicators.
#' @param weights weight of each sample.
#' @param strat stratification indicator.
#' @param alpha apha level to compute confidence interval.
#' @param method.test Statistical test to use in order to compute the p-values
#'   related to a D. index, see [survival::summary.coxph] for more details.
#' @param na.rm `TRUE` if missing values should be removed.
#' @param ... additional parameters to be passed to the [survival::coxph]
#'   function.
#'
#' @details
#' The D index is computed using the Cox model fitted on the scaled rankits of
#'   the risk scores instead of the risk scores themselves. The scaled rankits
#'   are the expected standard Normal order statistics scaled by
#'   kappa = sqrt(8/pi). See (Royston and Sauerbrei, 2004) for details.
#'   Note that the value D reported in (Royston and Sauerbrei, 2004) is given.
#'
#' @value
#' A list with items:
#' - d.index: D index (exponentiated, aka hazard ratio).
#' - coef: D index estimate (coefficient) fitted in the cox regression model.
#' - se: standard error of the estimate.
#' - lower: lower bound for the confidence interval.
#' - upper: upper bound for the confidence interval.
#' - p.value: p-value for the statistical test if the estimate if different from 0.5.
#' - n: number of samples used for the estimation.
#' - coxm: [survival::coxph.object] fitted on the survival data and z (see below).
#' - data: list of data used to compute the index (x, z, surv.time and
#' surv.event). The item z contains the scaled rankits which are the expected
#' standard Normal order statistics scaled by kappa.
#'
#' @authors
#' Benjamin Haibe-Kains
#'
#' @references
#' Royston, P. and Sauerbrei, W. (2004) "A new measure of prognostic separation
#'   in survival data", Statistics in Medicine, 23, pages 723–748.
#'
#' @seealso
#' [survival::coxph], [survival::coxph.object], [SuppDists::normOrder]
#'
#' @md
#' @export
D.index <-
function(x, surv.time, surv.event, weights, strat, alpha=0.05, method.test=c("logrank", "likelihood.ratio", "wald"), na.rm=FALSE, ...) {
	#require(SuppDists)
	method.test <- match.arg(method.test)
	if(!missing(weights)) {
		if(length(weights) != length(x)) { stop("bad length for parameter weights!") }
		## remove weights=0 because the coxph function does not deal with them properly
		iix <- weights <= 0
		if(any(iix)) { warning("samples with weight<=0 are discarded") }
		weights[iix] <- NA
	} else { weights <- rep(1,  length(x)) }
	if(!missing(strat)) {
		if(length(strat) != length(x)) { stop("bad length for parameter strat!") }
	} else { strat <- rep(1,  length(x)) }
	cc.ix <- complete.cases(x, surv.time, surv.event, weights, strat)
	if(sum(cc.ix) < 3) {
	## not enough observations
		data <- list("x"=x, "z"=rep(NA, length(x)), "surv.time"=surv.time, "surv.event"=surv.event, "weights"=weights, "strat"=strat)
		return(list("d.index"=NA, "coef"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=sum(cc.ix), "coxm"=NA, "data"=data))
	}
	if(any(!cc.ix) & !na.rm) { stop("NA values are present!") }
	sx <- x[cc.ix]
	oo <- order(sx, decreasing=FALSE)
	sx <- sx[oo]
	stime <- surv.time[cc.ix][oo]
	sevent <- surv.event[cc.ix][oo]
	sweights <- weights[cc.ix][oo]
	sstrat <- strat[cc.ix][oo]
	kap <- sqrt(8/pi)
	z <- kap^-1 * SuppDists::normOrder(N=length(sx))
	#ties?
	dup <- duplicated(sx)
	if(any(dup)) {
		udup <- unique(sx[dup])
		for(i in 1:length(udup)) { z[sx == udup[i]] <- mean(z[sx == udup[i]]) }
	}
	z2 <- x
	z2[!cc.ix] <- NA
	z2[cc.ix] <- z[match(1:length(sx), oo)]
	data <- list("x"=x, "z"=z2, "surv.time"=surv.time, "surv.event"=surv.event, "weights"=weights, "strat"=strat)
	#fit the cox model
	options(warn=2)
	rr <- try(coxph(Surv(stime, sevent) ~ strata(sstrat) + z, weights=sweights, ...))
	options(warn=0)
	if(class(rr) == "try-error") {
		res <- list("d.index"=NA, "coef"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=sum(cc.ix), "coxm"=NA, "data"=data)
	} else {
		dicoef <- rr$coefficients
		dise <- sqrt(drop(rr$var))
		names(dicoef) <- names(dise) <- NULL
		mystat <- NA
		switch(method.test,
		"logrank"={
			mystat <- rr$score
		},
		"likelihood.ratio"={
			mysat <- 2 * (rr$loglik[2] - rr$loglik[1])
		},
		"wald"={
			mystats <- rr$wald.test
			##(hrcoef / hrse)^2
		})
		mypp <- pchisq(mystat, df=1, lower.tail=FALSE)
		res <- list("d.index"=exp(dicoef), "coef"=dicoef, "se"=dise, "lower"=exp(dicoef - qnorm(alpha / 2, lower.tail=FALSE) * dise), "upper"=exp(dicoef + qnorm(alpha / 2, lower.tail=FALSE) * dise), "p.value"=mypp, "n"=rr$n, "coxm"=rr, "data"=data)
	}

	return(res)
}
