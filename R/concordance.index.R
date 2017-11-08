`concordance.index` <-
function(x, surv.time, surv.event, cl, weights, comppairs=10, strat, alpha=0.05, outx=TRUE, method=c("conservative", "noether", "nam"), alternative=c("two.sided", "less", "greater"), na.rm=FALSE) {
	method <- match.arg(method)
	alternative <- match.arg(alternative)
	if(!missing(weights)) {
		if(length(weights) != length(x)) { stop("bad length for parameter weights!") }
		if(min(weights, na.rm=TRUE) < 0 && max(weights, na.rm=TRUE) > 1) { stop("weights must be a number between 0 and 1!") }
	} else { weights <- rep(1, length(x)) }
	if(!missing(strat)) {
		if(length(strat) != length(x)) { stop("bad length for parameter strat!") }
	} else { strat <- rep(1, length(x)) }

	if(missing(cl) && (missing(surv.time) || missing(surv.event))) { stop("binary classes and survival data are missing!") }
	if(!missing(cl) && (!missing(surv.time) || !missing(surv.event))) { stop("choose binary classes or survival data but not both!") }
	msurv <- FALSE
	if(missing(cl)) { ## survival data
		msurv <- TRUE
		cl <- rep(0, length(x))
	} else { surv.time <- surv.event <- rep(0, length(x)) } ## binary classes

	cc.ix <- complete.cases(x, surv.time, surv.event, cl, weights, strat)
	if(sum(cc.ix) < 3) {
	## not enough observations
		if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
		return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=sum(cc.ix), "data"=data, "comppairs"=NA))
	}
	if(any(!cc.ix) & !na.rm) { stop("NA values are present!") }
  ## remove samples whose the weight is equal to 0 to speed up the computation of the concordance index
	cc.ix <- cc.ix & weights != 0
	x2 <- x[cc.ix]
	cl2 <- cl[cc.ix]
	st <- surv.time[cc.ix]
	se <- surv.event[cc.ix]
	if(msurv && sum(se) == 0) {
		warning("\nno events, the concordance index cannot be computed!")
		data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event)
		return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=0, "data"=data, "comppairs"=NA))
	}
	if(!msurv && length(unique(cl2)) == 1) {
		warning("\nonly one class, the concordance index cannot be computed!")
		data  <- list("x"=x, "cl"=cl)
		return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=0, "data"=data, "comppairs"=NA))
	}
	weights <- weights[cc.ix]
	strat <- strat[cc.ix]
	strat <- as.numeric(as.factor(strat))
	ustrat <- sort(unique(strat)) ## to check later
	N <- sum(weights) ##length(x2)
	if(N <= 1) {
    warning("\nweights of observations are too small (sum should be > 1), the concordance index cannot be computed!")
    if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
    return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=length(x2), "data"=data, "comppairs"=NA))
  }

	ch <- dh <- uh <- rph <- rep(0, times=length(strat))
	lenS <- length(strat)
	lenU <- length(ustrat)
	out <- .C(.C_concordanceIndexC, as.integer(as.logical(msurv)), as.integer(ustrat), as.double(x2),
			as.integer(cl2), as.double(st), as.integer(se), as.double(weights), as.integer(strat),
			as.integer(N), as.integer(as.logical(outx)), ch = as.numeric(ch), dh = as.numeric(dh),
			uh = as.numeric(uh), rph = as.numeric(rph), as.integer(lenS), as.integer(lenU), PACKAGE="survcomp")
  ch <- out$ch
  dh <- out$dh
  uh <- out$uh
  rph <- out$rph
  cscount <- sum(ch + dh) ## comparable pairs
  if(sum(ch)==0 || sum(dh)==0 || sum(ch * (ch - 1))==0 || sum(dh * (dh - 1))==0 || sum(ch * dh)==0 || cscount < comppairs){
    if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
    return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=length(x2), "data"=data, "comppairs"=cscount))
  }

	pc <- (1 / (N * (N - 1))) * sum(ch)
	pd  <- (1 / (N * (N - 1))) * sum(dh)
	cindex <- pc / (pc + pd)

	switch(method,
	"noether"={
    pcc <- (1 / (N * (N - 1) * (N - 2))) * sum(ch * (ch - 1))
    pdd <- (1 / (N * (N - 1) * (N - 2))) * sum(dh * (dh - 1))
    pcd <- (1 / (N * (N - 1) * (N - 2))) * sum(ch * dh)
    varp <- (4 / (pc + pd)^4) * (pd^2 * pcc - 2 * pc * pd * pcd + pc^2 * pdd)
    if((varp / N) > 0) {
      ci <- qnorm(p=alpha / 2, lower.tail=FALSE) * sqrt(varp / N)
      lower <- cindex - ci
      upper <- cindex + ci
      switch(alternative,
		   "two.sided"={ p <- pnorm((cindex - 0.5) / sqrt(varp / N), lower.tail=cindex < 0.5) * 2 },
		   "less"={ p <- pnorm((cindex - 0.5) / sqrt(varp / N), lower.tail=TRUE) },
		   "greater"={  p <- pnorm((cindex - 0.5) / sqrt(varp / N), lower.tail=FALSE) }
	    )
    } else { ci <- lower <- upper <- p <- NA }
  },
	"conservative"={
    C <- cindex
    ## pc and pd have been computed previously
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
  if(!is.na(varp) && (varp / N) > 0) { se <- sqrt(varp / N) } else { se <- NA }
  return(list("c.index"=cindex, "se"=se, "lower"=lower, "upper"=upper, "p.value"=p, "n"=length(x2), "data"=data, "comppairs"=cscount))
}
