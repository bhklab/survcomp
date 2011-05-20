`tdrocc` <-
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
    
	rocco <- survivalROC:::survivalROC(Stime=surv.time2, status=surv.event2, marker=x2, predict.time=time, cut.values=cutpts, entry=surv.entry, span=span, lambda=lambda, ...)
	res <- list("spec"=1-rocco$FP, "sens"=rocco$TP, "rule"=myrule, "cuts"=cutpts)
	
	return(list("spec"=1-rocco$FP, "sens"=rocco$TP, "rule"=myrule, "cuts"=cutpts, "time"=time, "survival"=rocco$Survival, "AUC"=rocco$AUC, "data"=data))
}

