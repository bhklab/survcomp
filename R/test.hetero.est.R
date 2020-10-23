test.hetero.est <-
function(x, x.se, na.rm=FALSE) {
	cc.ix <- complete.cases(x, x.se);
	if(!all(cc.ix) && !na.rm) { stop("missing values are present!"); }
	x <- x[cc.ix];
	k <- length(x);
	if(k == 1) {
		Q <- NA;
		qpv <- NA;	
	}
	else {
		x.se <- x.se[cc.ix];
		wi <- x.se^-2;
		Q <- sum(wi * x^2) - (sum(wi * x))^2 / sum(wi);
		qpv <- pchisq(Q, df=k-1, lower.tail=FALSE);
	}
	return(list("Q"=Q, "p.value"=qpv));
}