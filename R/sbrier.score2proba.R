#' @title Function to compute the BSCs from a risk score, for all the times of
#'   event occurrence
#'
#' @description
#' The function computes all the Brier scores (BSC) and the corresponding
#'   integrated Briser score (IBSC) from a risk score, for all the times of
#'   event occurrence. The risk score is first transformed in survival
#'   probabilities using either a Cox model or the product-limit estimator.
#'
#' @param data.tr the data frame for the training set. This data frame must
#'   contain three columns for the times, the event occurrence and the risk
#'   score. These columns are called "time", "event" and "score" respectively.
#' @param data.ts the data frame for the test set. This data frame must contain
#'   three columns for the times, the event occurrence and the risk score.
#'   These columns are called "time", "event" and "score" respectively.
#' @param method method for survival probabilities estimation using either a
#'   Cox model or the product-limit estimator
#'
#' @return
#' A list containing items:
#' - time: vector of points in time
#' - bsc: vector of Brier scores (BSC) at ome points in time
#' - bsc.integrated: value of the integrated Brier score (IBSC)
#'
#' @references
#' Brier, G. W. (1950) "Verification of forecasts expressed in terms of
#'   probabilities", Monthly Weather Review, 78, pages 1–3.
#'
#' Graf, E. and Schmoor, C. and Sauerbrei, W. and Schumacher, M. (1999)
#'   "Assessment and comparison of prognostic classification schemes for
#'   survival data ", Statistics in Medicine, 18, pages 2529–2545.
#'
#' Cox, D. R. (1972) "Regression Models and Life Tables", Journal of the Royal
#'   Statistical Society Series B, 34, pages 187–220.
#'
#' Andersen, P. K. and Borgan, O. and Gill, R. D. and Keiding, N. (1993)
#'   "Statistical Models Based on Counting Processes", Springer.
#'
#' @seealso
#' [ipred::sbrier], [survival::coxph], [prodlim::prodlim]
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
#' sbrier.score2proba(data.tr=dd, data.ts=dd, method="cox")
#' #product-limit estimator
#' sbrier.score2proba(data.tr=dd, data.ts=dd, method="prodlim")
#'
#' @md
#' @export
sbrier.score2proba <-
function(data.tr, data.ts, method=c("cox", "prodlim")) {
	## require(ipred)
	method <- match.arg(method)
	## remove missing values and sort the data for the test set
	cc.ix <- complete.cases(data.ts)
	ot <- order(data.ts$time)[1:(length(cc.ix)-sum(!cc.ix))]
	data.ts <- data.ts[ot, ,drop=FALSE]
	surv.time.ts <- data.ts$time
	surv.event.ts <- data.ts$event
	score.ts <- data.ts$score
	btime <- surv.time.ts[surv.time.ts >= 0 & surv.time.ts <= max(surv.time.ts, na.rm=TRUE)]
	utime <- unique(surv.time.ts[surv.event.ts == 1])
	bsc <- rep(NA, length(btime))
	switch(method,
	"cox"={
		##require(survival)
		## fit the cox model for the training set
		coxm <- survival::coxph(Surv(time, event) ~ score, data=data.tr)
		## compute survival probabilities using the cox model fitted on the training set and the score from the test set
		#sf <- survfit(coxm, newdata=data.ts)
		dd <- data.frame("score"=score.ts)
		sf <- survfit(coxm, newdata=dd)
		for(i in 1:length(utime)) {
			mypred <- getsurv2(sf=sf, time=utime[i])
			bsc[is.na(bsc) & btime <= utime[i]] <- ipred::sbrier(obj=Surv(surv.time.ts, surv.event.ts), pred=mypred, btime=utime[i])
		}	
	},
	"prodlim"={
		#require(KernSmooth)
		prodlim.m <- prodlim::prodlim(Surv(time, event) ~ score, data=data.tr)
		lpred <- predict(prodlim.m, newdata=data.ts, times=utime)
		names(lpred) <- dimnames(data.ts)[[1]]
		bsc <- rep(NA, length(btime))
		for(i in 1:length(utime)) {
			mypred <- unlist(lapply(lpred, function(x, ix) { return(x[[ix]]) }, ix=i))
			bsc[is.na(bsc) & btime <= utime[i]] <- ipred::sbrier(obj=Surv(surv.time.ts, surv.event.ts), pred=mypred, btime=utime[i])
		}
	})
	if(sum(is.na(bsc)) > 0) { bsc[is.na(bsc)] <- bsc[ min(which(is.na(bsc)))-1] } 
	diffs <- c(btime[1], btime[2:length(btime)] - btime[1:(length(btime) - 1)])
	bsc.int <- sum(diffs * bsc)/max(btime)
	return(list("time"=btime, "bsc"=bsc, "bsc.integrated"=bsc.int))
}

