#' @name fisherz
#'
#' @title Function to compute Fisher z transformation
#'
#' @description
#' The function computes the Fisher z transformation useful to calculate the
#'   confidence interval of Pearson's correlation coefficient.
#'
#' @usage fisherz(x, inv = FALSE, eps = 1e-16)
#'
#' @param x value, e.g. Pearson's correlation coefficient
#' @param inv `TRUE` for inverse Fisher z transformation, `FALSE` otherwise
#' @param eps tolerance for extreme cases, i.e. \eqn{|x| ~ 1} when `inv=FALSE`
#'   and \eqn{|x| ~ inf} when `inv=TRUE`.
#'
#' @details
#' The sampling distribution of Pearson's ρ is not normally distributed.
#'   R. A. Fisher developed a transformation now called “Fisher's z
#'   transformation” that converts Pearson's ρ to the normally distributed
#'   variable z. The formula for the transformation is
#'   \eqn{z = 1 / 2 [ \log(1 + ρ) - \log(1 - ρ) ]}.
#'   Two attributes of the distribution of the z statistic: (1) It is normally
#'   distributed and (2) it has a known standard error of \eqn{σ_z = 1 / √{N - 3}}
#'   where \eqn{N} is the number of samples. Fisher's z is used for computing
#'   confidence intervals on Pearson's correlation and for confidence intervals
#'   on the difference between correlations.
#'
#' @return [`numeric`] Fisher's z statistic
#'
#' @authors
#' Benjamin Haibe-Kains
#'
#' @references
#' R. A. Fisher (1915) "Frequency distribution of the values of the correlation
#'   coefficient in samples of an indefinitely large population". Biometrika,
#'   10,pages 507–521.
#'
#' @seealso
#' [stats::cor]
#'
#' @examples
#' set.seed(12345)
#' x1 <- rnorm(100, 50, 10)
#' x2 <- runif(100,.5,2)
#' cc <- cor(x1, x2)
#' z <- fisherz(x=cc, inv=FALSE)
#' z.se <- 1 / sqrt(100 - 3)
#' fisherz(z, inv=TRUE)
#'
#' @md
#' @export
fisherz <-
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