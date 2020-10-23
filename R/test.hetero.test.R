#' @name
#'
#'
#'
#'
#'
#' @md
#' @export
test.hetero.test <-
function(p, weight, na.rm=FALSE) {
	k <- length(p);
	if(missing(weight)) { weight <- rep(1, k); }
	cc.ix <- !is.na(p);
	if(!all(cc.ix) && !na.rm) { stop("missing values are present!"); }
	p <- p[cc.ix];
	weight <- weight[cc.ix];
	z <- qnorm(p, lower.tail=FALSE);
	Q <- sum(weight * (z - mean(z))^2)
	qpv <- pchisq(Q, df=k-1, lower.tail=FALSE);
	return(list("Q"=Q, "p.value"=qpv));
}