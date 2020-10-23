#' @name combine.test
#'
#' @title Function to combine probabilities
#'
#' @description
#' The function combines several p-value estimated from the same null hypothesis
#'   in different studies involving independent data.
#'
#' @usage combine.test(p, weight, method = c("fisher", "z.transform", "logit"),
#'   hetero = FALSE, na.rm = FALSE)
#'
#' @param p vector of p-values
#' @param weight vector of weights (e.g. sample size of each study)
#' @param hetero fisher for the Fisher's combined probability test, z.transform
#'   for the Z-transformed test, logit for the weighted Z-method
#' @param na.rm TRUE if the missing values should be removed from the data,
#'   FALSE otherwise
#'
#' @details The p-values must be one-sided and computed from the same null
#'   hypothesis.
#'
#' @return [`numeric`] p-value
#'
#' @references
#' Whitlock, M. C. (2005) "Combining probability from independent tests: the
#'   weighted Z-method is superior to Fisher's approach", J. Evol. Biol., 18,
#'   pages 1368–1373.
#'
#' @seealso
#' [survcomp::test.hetero.test]
#'
#' @examples
#' p <- c(0.01, 0.13, 0.07, 0.2)
#' w <- c(100, 50, 200, 30)
#'
#' #with equal weights
#' combine.test(p=p, method="z.transform")
#' #with p-values weighted by the sample size of the studies
#' combine.test(p=p, weight=w, method="z.transform")
#'
#' @md
#' @export
combine.test <-
function(p, weight, method=c("fisher", "z.transform", "logit"), hetero=FALSE, na.rm=FALSE) {
	if(hetero) { stop("function to deal with heterogeneity is not implemented yet!") }
	method <- match.arg(method)
	na.ix <- is.na(p)
	if(any(na.ix) && !na.rm) { stop("missing values are present!") }
	if(all(na.ix)) { return(NA) } ## all p-values are missing
	p <- p[!na.ix]
	k <- length(p)
	if(k == 1) { return(p) }
	if(missing(weight)) { weight <- rep(1, k); }
	switch(method,  
	"fisher"={
		cp <- pchisq(-2 * sum(log(p)), df=2*k, lower.tail=FALSE)
	}, 
	"z.transform"={
		z <- qnorm(p, lower.tail=FALSE)
		cp <- pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
	}, 
	"logit"={
		tt <- (- sum(log(p / (1 - p)))) / sqrt(k * pi^2 * (5 * k + 2) / (3 * (5 * k + 4)))
		cp <- pt(tt,df=5*k+4, lower.tail=FALSE)
	})
	return(cp)
}