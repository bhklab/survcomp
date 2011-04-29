'fisherz' <-
function(x, inv=FALSE, eps=1e-16) {
	
	myfoo <- function(x, inv, eps) {
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
	
	if(is.matrix(x) || is.data.frame(x)) { return(apply(X=x, MARGIN=c(1, 2), FUN=myfoo, inv=inv, eps=eps)) }
	if(is.vector(x)) { return(sapply(X=x, FUN=myfoo, inv=inv, eps=eps)) }
	return(myfoo(x=x, inv=inv, eps=eps))
}