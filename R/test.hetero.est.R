#' @name test.hetero.est
#'
#' @title Function to test the heterogeneity of set of probabilities
#'
#' @description
#' The heterogeneity test is known to be very conservative. Consider a
#'   p-value < 0.1 as significant.
#'
#' @param x vector of estimates
#' @param x.se vector of standard errors of the corresponding estimates
#' @param na.rm `TRUE` if the missing values should be removed from the data,
#'   `FALSE` otherwise
#'
#' @details
#' The heterogeneity test is known to be very conservative. Consider a
#'   p-value < 0.1 as significant.
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
#' @seealso
#' [survcomp::combine.est]
#'
#' @examples
#' set.seed(12345)
#' x1 <- rnorm(100, 50, 10) + rnorm(100, 0, 2)
#' m1 <- mean(x1)
#' se1 <- sqrt(var(x1))
#' x2 <- rnorm(100, 75, 15) + rnorm(100, 0, 5)
#' m2 <- mean(x2)
#' se2 <- sqrt(var(x2))
#'
#' test.hetero.est(x=c(m1, m2), x.se=c(se1, se2))
#'
#' @md
#' @export
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