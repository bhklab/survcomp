'combine.est' <-
function(x, x.se, hetero=FALSE, na.rm=FALSE) {
	cc.ix <- complete.cases(x, x.se)
	if(!all(cc.ix) && !na.rm) { stop("missing values are present!") }
	if(all(!cc.ix)) { return(list("estimate"=NA, "se"=NA)) } ## all estimates / standard errors are missing
	x <- x[cc.ix]
	x.se <- x.se[cc.ix]
	k <- length(x)
	if(any(x.se == 0)) {
		warning("standard deviation of zero is present!")
		x.se[x.se == 0] <- 10^-16	
	}
	if(k == 1) { return(list("estimate"=x, "se"=x.se)) }
	wi <- x.se^-2
	if(hetero) {
		w.bar <- sum(wi / k)
		s2w <- (sum(wi^2) - k * w.bar^2) / (k - 1)
		U <- (k - 1) * (w.bar - s2w / (k * w.bar))
		Q <- test.hetero.est(x=x, x.se=x.se)$Q
		tau2 <- ifelse(Q <= (k - 1), 0, (Q - (k - 1)) / U)
		wi <- 1 / ((1 / wi) + tau2)
	}
	ce <- c(sum(wi * x) / sum(wi), sqrt(1/sum(wi)))
	return(list("estimate"=ce[1], "se"=ce[2]))
}
