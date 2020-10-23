#' @name hr.comp.meta
#'
#' @title Function to estimate the hazard ratio through Cox regression
#'
#' @description
#' Function to compute the hazard ratio for a risk prediction.
#'
#' @usage
#' hazard.ratio(x, surv.time, surv.event, weights, strat, alpha = 0.05,
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
#' @param ... `r print("additional parameters to be passed to [survival::coxph] function")`
#'
#' @details
#' The hazard ratio is computed using the Cox model.
#'
#' @return
#' A list with items:
#' - hazard.ratio: hazard ratio estimate.
#' - coef: coefficient (beta) estimated in the cox regression model.
#' - se: standard error of the coefficient (beta) estimate.
#' - lower: lower bound for the confidence interval.
#' - upper: upper bound for the confidence interval.
#' - p.value: p-value computed using the likelihood ratio test whether the
#' hazard ratio is different from 1.
#' - n: number of samples used for the estimation.
#' - coxm: [survival::coxph.object] fitted on the survival data and x (see below).
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
hr.comp.meta <-
function(list.hr1, list.hr2, hetero=FALSE) {

	if(length(list.hr1) != length(list.hr2)) { stop("the concordance indices are computed from different number of samples!") }

	n <- 0
	x1 <- x1.se <- x2 <- x2.se <- corz <- corz.se <- NULL
	for(i in 1:length(list.hr1)) {
		nn <- list.hr1[[i]]$n
		if(nn != list.hr2[[i]]$n) { stop("the number of samples to compute the concordance indices is not the same!") }
		n <- n + nn
		x1 <- c(x1, list.hr1[[i]]$coef)
		x1.se <- c(x1.se, list.hr1[[i]]$se)
		x2 <- c(x2, list.hr2[[i]]$coef)
		x2.se <- c(x2.se, list.hr2[[i]]$se)
		cort <- cor(list.hr1[[i]]$data$x, list.hr2[[i]]$data$x, method="spearman", use="complete.obs")
		## since r is the spearman correlation coefficient and not the Pearson's one, we should apply a correction factor (see http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient for details)
		corz <- c(corz, sqrt((nn - 3) / 1.06) * fisherz(cort, inv=FALSE))
		if(nn > 3) { corz.se <- c(corz.se, 1 / sqrt(nn - 3)) } else { corz.se <- c(corz.se, NA) }
	}
	x1.meta <- combine.est(x=x1, x.se=x1.se, hetero=hetero, na.rm=TRUE)
	x2.meta <- combine.est(x=x2, x.se=x2.se, hetero=hetero, na.rm=TRUE) 
	if(x1.meta$estimate == x2.meta$estimate && x1.meta$se == x2.meta$se) {
	## same hazard ratios
		return(list("p.value"=1, "cindex1"=x1.meta$estimate, "cindex2"=x2.meta$estimate))
	}
	rz <- combine.est(x=corz, x.se=corz.se, na.rm=TRUE, hetero=hetero)$estimate
	## since r is the spearman correlation coefficient and not the Pearson's one, we should apply a correction factor (see http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient for details)
	rz <- rz / (sqrt((n - 3) / 1.06))
	r <- fisherz(rz, inv=TRUE)

	if(abs(r) < 1) {
		t.stat <- (x1.meta$estimate - x2.meta$estimate) / sqrt(x1.meta$se^2 + x2.meta$se^2 - 2 * r * x1.meta$se * x2.meta$se)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "hr1"=exp(x1.meta$estimate), "hr2"=exp(x2.meta$estimate)))
}
