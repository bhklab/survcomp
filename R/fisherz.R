'fisherz' <-
function(x, inv=FALSE,  eps=1e-10) {
	if(is.na(x)) { return(NA) }
	if(!inv) {
		if((1 - abs(x)) < eps) { x <- ifelse(x < 0, -Inf, Inf) }
		else { x <- (log(1 + x) - log(1 - x)) / 2 }
	}
	else {
		if(is.infinite(x) || x > (1 / eps)) { x <- ifelse(x < 0, -1, 1) }
		else { x <- (exp(2 * x) - 1) / (exp(2 * x) + 1) }
	}
	return(x)
}