#' @name td.sens.spec
#'
#' @title Function to compute sensitivity and specificity for a binary
#'   classification of survival data
#'
#' @details
#' The function is a wrapper for the survivalROC.C function in order to compute
#'   sensitivity and specificity for a binary classification of survival data.
#'
#' @usage
#' td.sens.spec(cl, surv.time, surv.event, time, span = 0, sampling = FALSE,
#'   na.rm = FALSE, ...)
#'
#' @param cl vector of binary classes.
#' @param surv.time vector of times to event occurrence.
#' @param surv.event vector of event occurrence indicators.
#' @param time time point for sensitivity and specificity estimations.
#' @param span Span for the NNE. Default value is 0.
#' @param sampling jackknife procedure to estimate the standard error of
#'   sensitivity and specificity estimations.
#' @param na.rm TRUE if the missing values should be removed from the data,
#'   FALSE otherwise.
#' @param ... additional arguments to be passed to the
#'   [survivalROC::survivalROC] function.
#'
#' @details
#' Only NNE method is used to estimate sensitivity and specificity
#'   (see [survivalROC::survivalROC.C]). The standard error for sensitivity and
#'   specificity is estimated through jackknife procedure (see
#'   [bootstrap::jackknife]).
#'
#' @return
#' A list of items:
#' - sens: sensitivity estimate
#' - sens.se: standard error for sensitivity estimate
#' - spec: specificity estimate
#' - spec.se: standard error for specificity estimate
#'
#' @references
#' Heagerty, P. J. and Lumley, T. L. and Pepe, M. S. (2000) "Time-Dependent ROC
#'   Curves for Censored Survival Data and a Diagnostic Marker", Biometrics,
#'   56, pages 337–344.
#'
#' Efron, B. and Tibshirani, R. (1986). "The Bootstrap Method for standard
#'   errors, confidence intervals, and other measures of statistical accuracy",
#'   Statistical Science, 1 (1), pages 1–35.
#'
#' @seealso
#' [survivalROC::survivalROC]
#'
#' @examples
#' set.seed(12345)
#' gender <- sample(c(0,1), 100, replace=TRUE)
#' stime <- rexp(100)
#' cens <- runif(100,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' mysenspec <- td.sens.spec(cl=gender, surv.time=stime, surv.event=sevent,
#'   time=1, span=0, na.rm=FALSE)
#'
#' @md
#' @export
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