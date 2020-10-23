#' @name combine.est
#'
#' @title Function to combine estimates
#'
#' @description
#' The function combines several estimators using meta-analytical formula to
#'   compute a meta-estimate.
#'
#' @usage combine.est(x, x.se, hetero = FALSE, na.rm = FALSE)
#'
#' @param x vector of estimates
#' @param x.se vector of standard errors of the corresponding estimates
#' @param hetero `TRUE` is the heterogeneity should be taken into account
#'   (random effect model), `FALSE` otherwise (fixed effect model)
#' @param na.rm TRUE if the missing values should be removed from the data,
#'   `FALSE` otherwise.
#'
#' @return
#' A list with items:
#' - estimate: meta-estimate
#' - se: standard error of the meta-estimate
#'
#' @authors
#' Benjamin Haibe-Kains
#'
#' @references
#' Cochrane, W. G. (1954) "The combination of estimates from different
#'   experiments", Biometrics, 10, pages 101–129.
#'
#' @seealso
#' [survcomp::test.hetero.est]
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
#' #fixed effect model
#' combine.est(x=c(m1, m2), x.se=c(se1, se2), hetero=FALSE)
#' #random effect model
#' combine.est(x=c(m1, m2), x.se=c(se1, se2), hetero=TRUE)
#'
#' @md
#' @export
combine.est <-
function(x, x.se, hetero=FALSE, na.rm=FALSE) {
	cc.ix <- complete.cases(x, x.se)
	if(!all(cc.ix) && !na.rm) { stop("missing values are present!") }
	if(all(!cc.ix)) { return(list("estimate"=NA, "se"=NA)) } ## all estimates/standard errors are missing
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
