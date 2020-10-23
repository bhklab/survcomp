#' @name hazard.ratio
#'
#' @title Function to estimate the hazard ratio through Cox regression
#'
#' @description
#' Function to compute the hazard ratio for a risk prediction.
#'
#' @usage
#' hazard.ratio(x, surv.time, surv.event, weights, strat, alpha = 0.05,
#' method.test = c("logrank", "likelihood.ratio", "wald"), na.rm = FALSE, ...)
#'
#' @param x a vector of risk predictions.
#' @param surv.time a vector of event times.
#' @param surv.event a vector of event occurrence indicators.
#' @param weights weight of each sample.
#' @param strat stratification indicator.
#' @param alpha alpha level to compute confidence interval.
#' @param method.test Statistical test to use in order to compute the p-values
#'   related to a D.index, see [survival::summary.coxph] for more details.
#' @param na.rm `TRUE` if missing values should be removed.
#' @param ... additional parameters to be passed to the [survival::coxph]
#'   function.
#'
#' @details
#' The hazard ratio is computed using the Cox model.
#'
#' @return
#' A list with items:
#' - hazard.ratio: hazard ratio estimate.
#' - coef: coefficient (beta) estimated in the cox regression model.
#' - se: standard error of the coefficient (beta) estimate.
#' - lower: upper bound for the confidence interval.
#' - upper: upper bound for the confidence interval.
#' - p.value: p-value computed using the likelihood ratio test whether the hazard ratio is different from 1.
#' - n: number of samples used for the estimation.
#' - data: list of data used to compute the hazard ratio (x, surv.time and
#' surv.event).
#'
#' @references
#' Cox, D. R. (1972) "Regression Models and Life Tables", Journal of the Royal
#'   Statistical Society Series B, 34, pages 187–220.
#'
#' @seealso
#' [survival::coxph], [survival::coxph.object]
#'
#' @examples
#' set.seed(12345)
#' age <- rnorm(100, 50, 10)
#' stime <- rexp(100)
#' cens   <- runif(100,.5,2)
#' sevent  <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' strat <- sample(1:3, 100, replace=TRUE)
#' weight <- runif(100, min=0, max=1)
#' hazard.ratio(x=age, surv.time=stime, surv.event=sevent, weights=weight,
#'   strat=strat)
#'
#' @md
#' @export
hazard.ratio <-
function(x, surv.time, surv.event, weights, strat, alpha=0.05, method.test=c("logrank", "likelihood.ratio", "wald"), na.rm=FALSE, ...) {
	method.test <- match.arg(method.test)
	if(!missing(weights)) {
		if(length(weights) != length(x)) { stop("bad length for parameter weights!") }
	} else { weights <- rep(1,  length(x)) }
	if(!missing(strat)) {
		if(length(strat) != length(x)) { stop("bad length for parameter strat!") }
		## remove weights=0 because the coxph function does not deal with them properly
		iix <- weights <= 0
		if(any(iix)) { warning("samples with weight<=0 are discarded") }
		weights[iix] <- NA
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
