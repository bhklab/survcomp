#' @name score2proba
#'
#' @title Function to compute the survival probabilities from a risk score
#'
#' @description
#' the function uses either a Cox model or the product-limit estimator to
#'   compute the survival probabilities from a risk score for a specific
#'   point in time.
#'
#' @param data.tr
#' @param score
#' @param yr
#' @param method
#' @param conf.int
#' @param which.est
#'
#' @return vector of predicted survival probabilities
#'
#' @seealso
#' [survival::coxph], [prodlim::prodlim]
#'
#' @references
#' Cox, D. R. (1972) "Regression Models and Life Tables", Journal of the Royal
#'   Statistical Society Series B, 34, pages 187–220.
#'
#' Andersen, P. K. and Borgan, O. and Gill, R. D. and Keiding, N. (1993)
#'   "Statistical Models Based on Counting Processes", Springer.
#'
#' @examples
#' set.seed(12345)
#' age <- rnorm(30, 50, 10)
#' stime <- rexp(30)
#' cens <- runif(30,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' dd <- data.frame("time"=stime, "event"=sevent, "score"=age)
#'
#' #Cox's model
#' score2proba(data.tr=dd, score=dd$score, yr=1, method="cox")
#' #product-limit estimator
#' score2proba(data.tr=dd, score=dd$score, yr=1, method="prodlim")
#'
#' @md
#' @export
score2proba <-
function(data.tr, score, yr, method=c("cox", "prodlim"), conf.int=0.95, which.est=c("point", "lower", "upper")) {
	method <- match.arg(method)
	which.est <- match.arg(which.est)
	cc.ix <- complete.cases(score)
	score2 <- score[cc.ix]
	pred <- rep(NA, length(score))
	names(pred) <- names(score)
	switch(method,
	"cox"={
		predm <- coxph(Surv(time, event) ~ score, data=data.tr)
		sf <- survfit(predm, newdata=data.frame("score"=score2), conf.int=conf.int)
		pred[cc.ix] <- getsurv2(sf=sf, time=yr, which.est=which.est)
	},
	"prodlim"={
		#require(prodlim)
		#require(KernSmooth)
		if(which.est != "point") { stop("not implemented yet!") }
		predm <- prodlim::prodlim(Surv(time, event) ~ score, data=data.tr, conf.int=conf.int)
		pred[cc.ix] <- unlist(predict(predm, newdata=data.frame("score"=score2), times=yr))
	})
	return(pred)
}

