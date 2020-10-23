#' @name bhr.comp
#'
#' @title Function to statistically compare two balanced hazard ratios
#'
#' @description
#' This function compares two balanced hazard ratios from their betas and
#'   standard errors as computed by a Cox model for instance. The statistical
#'   test is a Student t test for dependent samples. The two balanced hazard
#'   ratios must be computed from the same survival data.
#'
#' @details
#' The two balanced hazard ratios must be computed from the same samples
#'   (and corresponding survival data). The function uses a Student t test for
#'   dependent samples.
#'
#' @param bhr1 first balanced hazard ratio.
#' @param bhr2 second balanced hazard ratio.
#'
#' @return
#' A list with items:
#' - p.value:	p-value from the Student t test for the comparison beta1 >
#' beta2 (equivalently bhr1 > bhr2)
#' - bhr1: value of the first balanced hazard ratio
#' - bhr2: value of the second balanced hazard ratio
#'
#' @references
#' Student. (1908) "The Probable Error of a Mean", Biometrika, 6, 1, pages 1–25.
#' Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008)
#'   "A comparative study of survival models for breast cancer prognostication
#'   based on microarray data: does a single gene beat them all?",
#'   Bioinformatics, 24, 19, pages 2200–2208.
#' Branders, S. and Dupont, P. (2015) "A balanced hazard ratio for risk group
#'   evaluation from survival data", Statistics in Medicine, 34(17),
#'   pages 2528–2543.
#'
#' @seealso
#' [survcomp::balanced.hazard.ratio]
#' [survival::coxph]
#' [stats::t.test]
#'
#' @examples
#' set.seed(12345)
#' age <- as.numeric(rnorm(100, 50, 10) >= 50)
#' size <- as.numeric(rexp(100,1) > 1)
#' stime <- rexp(100)
#' cens <- runif(100,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' bhr1 <- balanced.hazard.ratio(x=age, surv.time=stime, surv.event=sevent)
#' bhr2 <- balanced.hazard.ratio(x=size, surv.time=stime, surv.event=sevent)
#' bhr.comp(bhr1=bhr1, bhr2=bhr2)
#'
#' @md
#' @export
bhr.comp <-
function(bhr1, bhr2) {
	if(bhr1$n != bhr2$n) { stop("the balanced hazard ratios are computed from different number of samples!") }
	n <- bhr1$n
	x1 <- bhr1$data$x
	x2 <- bhr2$data$x
	beta1 <- bhr1$coef
	beta2 <- bhr2$coef
	se1 <- bhr1$se
	se2 <- bhr2$se
	r <- cor(x1, x2, method="spearman", use="complete.obs")
	if(abs(r) < 1) {
		t.stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "bhr1"=exp(beta1), "bhr2"=exp(beta2)))
}

