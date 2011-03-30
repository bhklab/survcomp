`concordance.index` <-
function(x, surv.time, surv.event, cl, weights, strat, alpha=0.05, outx=TRUE, method=c("conservative", "noether", "nam"), na.rm=FALSE) {
	method <- match.arg(method)
	if(!missing(weights)) {
		if(length(weights) != length(x)) { stop("bad length for parameter weights!") }
	} else { weights <- rep(1, length(x)) }
	if(!missing(strat)) {
		if(length(strat) != length(x)) { stop("bad length for parameter strat!") }
	} else { strat <- rep(1, length(x)) }
	
	if(missing(cl) && (missing(surv.time) || missing(surv.event))) { stop("binary classes and survival data are missing!") }
	if(!missing(cl) && (!missing(surv.time) || !missing(surv.event))) { stop("choose binary classes or survival data but not both!") }
	msurv <- FALSE
	if(missing(cl)) { #survival data
		msurv <- TRUE
		cl <- rep(0, length(x))
	} else { surv.time <- surv.event <- rep(0, length(x)) } #binary classes
	
	cc.ix <- complete.cases(x, surv.time, surv.event, cl, weights, strat)
	if(all(!cc.ix)) {
		if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
		return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=0, "data"=data))	
	}
	if(any(!cc.ix) & !na.rm) { stop("NA values are present!") }
	#remove samples whose the weight is equal to 0
	cc.ix <- cc.ix & weights != 0
	x2 <- x[cc.ix]
	cl2 <- cl[cc.ix]
	st <- surv.time[cc.ix]
	se <- surv.event[cc.ix]
	weights <- weights[cc.ix]
	strat <- strat[cc.ix]
	strat <- as.numeric(as.factor(strat))
	ustrat <- sort(unique(strat)) ## to check later
	N <- length(x2)
	
	ch <- dh <- uh <- rph <- rep(0, times=length(strat))
	lenS <- length(strat)
	lenU <- length(ustrat)
	out <- .C("concordanceIndexC", as.integer(as.logical(msurv)), as.integer(ustrat), as.double(x2),
			as.integer(cl2), as.double(st), as.integer(se), as.double(weights), as.integer(strat),
			as.integer(N), as.integer(as.logical(outx)), ch = as.integer(ch),dh = as.integer(dh),
			uh = as.integer(uh), rph = as.integer(rph), as.integer(lenS), as.integer(lenU), PACKAGE="survcomp")
     ch <- out$ch
     dh <- out$dh
     uh <- out$uh
     rph <- out$rph
	pc <- (1 / (N * (N - 1))) * sum(ch)
	pd  <- (1 / (N * (N - 1))) * sum(dh)
	cindex <- pc / (pc + pd)
	
	switch(method,
	"noether"={
	pcc <- (1 / (N * (N - 1) * (N - 2))) * sum(ch * (ch - 1))
	pdd <- (1 / (N * (N - 1) * (N - 2))) * sum(dh * (dh - 1))
	pcd <- (1 / (N * (N - 1) * (N - 2))) * sum(ch * dh)
	varp <- (4 / (pc + pd)^4) * (pd^2 * pcc - 2 * pc * pd * pcd + pc^2 * pdd)
	ci <- qnorm(p=alpha / 2, lower.tail=FALSE) * sqrt(varp / N)
	lower <- cindex - ci
	upper <- cindex + ci
	p <- pnorm((cindex - 0.5) / sqrt(varp / N), lower.tail=cindex < 0.5)
	},
	"conservative"={
	C <- cindex
	sum.ch <- sum(ch)
	sum.dh <- sum(dh)
	pc <- (1 / (N * (N - 1))) * sum.ch
	pd  <- (1 / (N * (N - 1))) * sum.dh
	w <- (2 * qnorm(p=alpha / 2, lower.tail=FALSE)^2) / (N * (pc + pd))
	ci <- sqrt(w^2 + 4 * w * C * (1 - C)) / (2 * (1 + w))
	point <- (w + 2 * C) / (2 * (1 + w))
	lower <- point - ci
	upper <- point + ci
	cindex <- C
	p <- NA
	varp <- NA
	},
	"name"={
		stop("method not implemented!")	
	})
	#bound the confidence interval
	lower <- ifelse(lower < 0, 0, lower)
	lower <- ifelse(lower > 1, 1, lower)
	upper <- ifelse(upper < 0, 0, upper)
	upper <- ifelse(upper > 1, 1, upper)
	if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
	return(list("c.index"=cindex, "se"=sqrt(varp / N), "lower"=lower, "upper"=upper, "p.value"=p, "n"=N, "data"=data))
}
