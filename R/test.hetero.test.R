#' @title Function to test the heterogeneity of set of probabilities
#'
#' @description
#' The function tests whether a set of p-values are heterogeneous.
#'
#' @param p vector of p-values
#' @param weight vector of weights (e.g. sample size of each study)
#' @param na.rm `TRUE` if the missing values should be removed from the data,
#'   `FALSE` otherwise
#'
#' @details
#' The p-values should be one-sided and computed from the same null hypothesis.
#'
#' @return
#' A list with items:
#' - Q: Q statistic
#' - p.value: p-value of the heterogeneity test
#'
#' @references
#' Cochrane, W. G. (1954) "The combination of estimates from different
#'   experiments", Biometrics, 10, pages 101–129.
#'
#' Whitlock, M. C. (2005) "Combining probability from independent tests: the
#'   weighted Z-method is superior to Fisher's approach", J. Evol. Biol., 18,
#'   pages 1368–1373.
#'
#' @seealso
#' [survcomp::combie.est]
#'
#' @examples
#' p <- c(0.01, 0.13, 0.07, 0.2)
#' w <- c(100, 50, 200, 30)
#'
#' #with equal weights
#' test.hetero.test(p=p)
#' #with p-values weighted by the sample size of the studies
#' test.hetero.test(p=p, weight=w)
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