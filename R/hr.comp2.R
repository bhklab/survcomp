#' @title Function to statistically compare two hazard ratios (alternative
#'   interface)
#'
#' @description
#' This function compares two hazard ratios from their betas and standard errors
#'   as computed by a Cox model for instance. The statistical test is a Student
#'   t test for dependent samples. The two hazard ratios must be computed from
#'   the same survival data.
#'
#' @usage hr.comp2(x1, beta1, se1, x2, beta2, se2, n)
#'
#' @param x1 risk score used to estimate the first hazard ratio.
#' @param beta1 beta estimate for the first hazard ratio.
#' @param se1 standard error of beta estimate for the first hazard ratio.
#' @param x2 risk score used to estimate the second hazard ratio.
#' @param beta2 beta estimate for the second hazard ratio.
#' @param se2 standard error of beta estimate for the first hazard ratio.
#' @param n number of samples from which the hazard ratios were estimated.
#'
#' @details
#' The two hazard ratios must be computed from the same samples (and
#'   corresponding survival data). The function uses a Student t test for
#'   dependent samples.
#'
#' @return
#' A list with items:
#' - p.value: p-value from the Student t test for the comparison beta1 > beta2
#' (equivalently hr1 > hr2)
#' - hr1: value of the first hazard ratio
#' - hr2: value of the second hazard ratio
#'
#' @references
#' Student (1908). "The Probable Error of a Mean", Biometrika, 6, 1, pages 1–25.
#'
#' Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A
#'   comparative study of survival models for breast cancer prognostication
#'   based on microarray data: does a single gene beat them all?",
#'   Bioinformatics, 24, 19, pages 2200–2208.
#'
#' @seealso
#' [survival::coxph], [stats::t.test]
#'
#' @examples
#' set.seed(12345)
#' age <- as.numeric(rnorm(100, 50, 10) >= 50)
#' size <- as.numeric(rexp(100,1) > 1)
#' stime <- rexp(100)
#' cens <- runif(100,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' coxm1 <- coxph(Surv(stime, sevent) ~ age)
#' coxm2 <- coxph(Surv(stime, sevent) ~ size)
#' hr.comp2(x1=age, beta1=coxm1$coefficients, se1=drop(sqrt(coxm1$var)),
#'   x2=size, beta2=coxm2$coefficients, se2=drop(sqrt(coxm2$var)), n=100)
#'
#' @md
#' @export
hr.comp2 <-
function(x1, beta1, se1, x2, beta2, se2, n) {
	r <- cor(x1, x2, method="spearman", use="complete.obs")
	if(abs(r) < 1) {
		t.stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "hr1"=exp(beta1), "hr2"=exp(beta2)))
}

