#' @title Function to compute the log partial likelihood of a Cox model
#'
#' @usage
#' logpl(pred, surv.time, surv.event, strata, na.rm = FALSE, verbose = FALSE)
#'
#' @param surv.time vector of times to event occurrence
#' @param surv.event vector of indicators for event occurrence
#' @param pred linear predictors computed using the Cox model
#' @param strata stratification variable
#' @param na.rm TRUE if the missing values should be removed from the data,
#'   FALSE otherwise
#' @param verbose verbosity of the function
#'
#' @return
#' vector of two elements: `logpl` and event for the estimation of the log
#'   partial likelihood and the number of events, respectively
#'
#' @references
#' Cox, D. R. (1972) "Regression Models and Life Tables", Journal of the Royal
#'   Statistical Society Series B, 34, pages 187–220.
#'
#' @seealso
#' [survival::coxph], [survcomp::cvpl]
#'
#' @examples
#' set.seed(12345)
#' age <- rnorm(100, 50, 10)
#' stime <- rexp(100)
#' cens <- runif(100,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' dd <- data.frame("stime"=stime, "sevent"=sevent, "age"=age)
#' ##Cox model
#' coxm <- coxph(Surv(stime, sevent) ~ age, data=dd)
#' ##log partial likelihood of the null model
#' logpl(pred=rep(0, nrow(dd)), surv.time=stime, surv.event=sevent)
#' ##log partial likelihood of the Cox model
#' logpl(pred=predict(object=coxm, newdata=dd), surv.time=stime, surv.event=sevent)
#' ##equivalent to
#' coxm$loglik
#'
#' @md
#' @export
logpl <-
function(pred, surv.time, surv.event, strata, na.rm=FALSE, verbose=FALSE) {

	##############
	#internal function
	##############
	
	logpl1 <- function(pred, surv.time, surv.event, verbose=FALSE) {	
	
		n <- length(pred)
		r <- rank(surv.time)
		ita <- pred
		epita <- exp(ita)
		d <- rep(0, n)
		dono <- rep(0, n)
		for(i in 1:n) {
			d[i] <- sum(surv.event[r == r[i]])
			dono[i] <- sum(epita[r >= r[i]])
		}
		risk <- d/dono
		risk1 <- d/dono^{	2}
		culrisk1 <- culrisk <- rep(0, n)
		for(i in 1:n) {
			culrisk[i] <- sum(unique(risk[r <= r[i]]))
			culrisk1[i] <- sum(unique(risk1[r <= r[i]]))
		}
		lik <- sum((ita - log(dono)) * surv.event)
		res <- c(lik, sum(surv.event))
		names(res) <- c("logpl", "event")
		return(res)
	}
	
	##############
	
	## remove NA values
	if(missing(strata)) { strata <- rep(1, length(pred)) } 
	cc.ix <- complete.cases(surv.time, surv.event, pred, strata)
	surv.time <- surv.time[cc.ix]
	surv.event <- surv.event[cc.ix]
	pred <- pred[cc.ix]
	strata <- strata[cc.ix]
    n <- sum(cc.ix)
    if (!all(cc.ix) && !na.rm) { stop("NA values are present!") }
    if(verbose) { message(sprintf("%i cases are removed due to NA values",as.integer(sum(!cc.ix)))) }
    
    ss <- unique(strata)
    if(length(ss) < 2) {
    	res <- logpl1(surv.time=surv.time, surv.event=surv.event, pred=pred, verbose=verbose)
    }
    else {
    	res1 <- 0
    	res2 <- 0
    	for(i in 1:length(ss)) {
    		myx <- strata == ss[i]
    		rr <- logpl1(surv.time=surv.time[myx], surv.event=surv.event[myx], pred=pred[myx], verbose=verbose)
    		res1 <- res1 + rr[1]
    		res2 <- res2 + rr[2]
    	}
    	res <- c(res1, res2)
    	names(res) <- c("logpl", "event")
    }
    return(res)
}
