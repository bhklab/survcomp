`score2proba` <-
function(data.tr, score, yr, method=c("cox", "prodlim"), conf.int=0.95, which.est=c("point", "lower", "upper")) {
	method <- match.arg(method)
	which.est <- match.arg(which.est)
	cc.ix <- complete.cases(score)
	score2 <- score[cc.ix]
	pred <- rep(NA, length(score))
	names(pred) <- names(score)
	switch(method,
	"cox"={
		require(survival)
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

