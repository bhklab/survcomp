#' @title Function to compute time-dependent ROC curves
#' 
#' @description
#' The function is a wrapper for the survivalROC.C function in order to compute
#'   sensitivity and specificity for a binary classification of survival data.
#'
#' @param cl vector of binary classes.
#' @param surv.time vector of times to event occurrence.
#' @param surv.event vector of event occurrence indicators.
#' @param time time point for sensitivity and specificity estimations.
#' @param span Span for the NNE. Default value is 0.
#' @param sampling jackknife procedure to estimate the standard error of
#'   sensitivity and specificity estimations.
#' @param na.rm `TRUE` if the missing values should be removed from the data,
#'   `FALSE` otherwise.
#' @param ... additional arguments to be passed to the
#'   [survivalROC::survivalROC] function.
#'
#' @details
#' Only NNE method is used to estimate sensitivity and specificity
#'   (see [survivalROC::survivalROC.C]). The standard error for sensitivity and
#'   specificity is estimated through jackknife procedure
#'   (see [bootstrap::jackknife]).
#'
#' @return
#' A list with items:
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
tdrocc <-
function(x, surv.time, surv.event, surv.entry=NULL, time, cutpts=NA, na.rm=FALSE, verbose=FALSE, span=0, lambda=0, ...) {
	#require(survivalROC)	
	data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event)
	cc.ix <- complete.cases(x, surv.time, surv.event, surv.entry)
   if (!all(cc.ix) && !na.rm) { stop("NA values are present!") }
   if(verbose) { message(sprintf("%i cases are removed due to NA values",as.integer(sum(!cc.ix)))) }
    
   x2 <- x[cc.ix]
   surv.time2 <- surv.time[cc.ix]
   surv.event2 <- surv.event[cc.ix]
    
   if(!all(sort(unique(surv.event2)) == c(0, 1))) { stop("survival event variable must take values 0 or 1") }
   if(is.na(cutpts)) {
       ux <- unique(sort(x2))
       delta <- min(diff(ux))/2
       cutpts <- c(ux - delta, ux[length(ux)] + delta)
   }
   myrule <- function(x, thresh) { return(ifelse(x > thresh, 1, 0)) }
   
	if(all(time < surv.time2[surv.event2 == 1])) { return(list("spec"=NA, "sens"=NA, "rule"=myrule, "cuts"=cutpts, "time"=time, "survival"=NA, "AUC"=NA, "data"=data)) }
    
	rocco <- survivalROC::survivalROC(Stime=surv.time2, status=surv.event2, marker=x2, predict.time=time, cut.values=cutpts, entry=surv.entry, span=span, lambda=lambda, ...)
	res <- list("spec"=1-rocco$FP, "sens"=rocco$TP, "rule"=myrule, "cuts"=cutpts)
	
	return(list("spec"=1-rocco$FP, "sens"=rocco$TP, "rule"=myrule, "cuts"=cutpts, "time"=time, "survival"=rocco$Survival, "AUC"=rocco$AUC, "data"=data))
}

