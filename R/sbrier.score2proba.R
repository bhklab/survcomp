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

