getsurv2 <-
function(sf, time, which.est=c("point", "lower", "upper")) {
	which.est <- match.arg(which.est)
	switch(which.est,
	"lower"={
		survres <- as.matrix(sf$lower)	
	},
	"point"={
		survres <- as.matrix(sf$surv)
	},
	"upper"={
		survres <- as.matrix(sf$upper)
	})
	if(time >= max(sf$time)) {
		time.ix <- length(sf$time)
	} else { time.ix <- order(sf$time <= time)[1] - 1 }
	if(time.ix == 0) {
		res <- rep(1, ncol(sf$surv[time.ix, ]))
		names(res) <- dimnames(sf$surv)[[2]]
	} else { res <- survres[time.ix, ] }	
	return(res)
}

