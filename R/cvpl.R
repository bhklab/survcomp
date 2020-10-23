cvpl <-
function(x, surv.time, surv.event, strata, nfold=1, setseed, na.rm=FALSE, verbose=FALSE) {
	x <- as.data.frame(x)
	if(is.null(dimnames(x))) { dimnames(x) <- list(names(surv.time), "x") }
	if(missing(strata)) { strata <- rep(1, length(surv.time)) }
	## remove NA values
	cc.ix <- complete.cases(x, surv.time, surv.event, strata)
	surv.time <- surv.time[cc.ix]
	surv.event <- surv.event[cc.ix]
	x <- x[cc.ix, ,drop=FALSE]
	strata <- strata[cc.ix]
	nr <- sum(cc.ix)
	if (!all(cc.ix) && !na.rm) { stop("NA values are present!") }
	if(verbose) { message(sprintf("%i cases (%i cases are removed due to NA values)", nr, sum(!cc.ix))) }

	## k-fold cross-validation
	if(nfold == 1) {
		k <- 1
		nfold <- nr
	} else { k <- floor(nr / nfold) }

	## set the random seed to use the same data partition
	## for the cross-validation
	if (!missing(setseed)) {
		set.seed(setseed)
	}
	smpl <- sample(nr)
	res.cvpl <- 0
	conv <- pl <- NULL
	dd <- data.frame("stime"=surv.time, "sevent"=surv.event, "strat"=strata, x)
	for (i in 1:nfold) {
		#index of samples to hold out
		if(i == nfold) { s.ix <- smpl[c(((i - 1) * k + 1):nr)] } else { s.ix <- smpl[c(((i - 1) * k + 1):(i * k))] }
		## convergence ?
		#lwa <- options("warn")$warn
		#options("warn"=2)
		ff <- sprintf("Surv(stime, sevent) ~ strata(strat) + %s", paste(dimnames(dd)[[2]][4:ncol(dd)], collapse=" + "))
		try(m <- coxph(formula=formula(ff), data=dd[-s.ix, , drop=FALSE]))
		#options("warn"=lwa)
		if (class(m) != "try-error") {
			conv <- c(conv, TRUE)
			li <- m$loglik[2]
			mypred <- predict(object=m, newdata=dd)
			l <- logpl(surv.time=dd[ , "stime"], surv.event=dd[ , "sevent"], pred=mypred, strata=dd[ , "strat"])[1]
		} else {
			conv <- c(conv, FALSE)
			l <- NA
			li <- NA
		}

		res.cvpl <- res.cvpl + (li - l)
		pl <- c(pl, (li - l) / length(s.ix)) # dividing by the number of events instead?
	}
	res.cvpl <- res.cvpl / nr # dividing by the number of events instead?
	names(conv) <- names(pl) <- paste(rep("split", nfold), 1:nfold, sep=".")
	return (list("cvpl"=res.cvpl, "pl"=pl, "convergence"=conv, "n"=nr))
}

