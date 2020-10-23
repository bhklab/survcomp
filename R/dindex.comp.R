#' @name dindex.comp
#'
#' @title Function to compare two D indices
#'
#' @description
#' This function compares two D indices from their betas and standard errors as
#'   computed by a Cox model for instance. The statistical test is a Student
#'   t test for dependent samples. The two D indices must be computed using
#'   the [survcomp::D.index] function from the same survival data.
#'
#' @usage dindex.comp(dindex1, dindex2)
#'
#' @param dindex1 first D index
#' @param dindex2 second D index
#'
#' @details
#' The two D indices must be computed using the D.index function from the same
#'   samples (and corresponding survival data). The function uses a Student
#'   t test for dependent samples.
#'
#' @return
#' A list with items:
#' - p.value: p-value from the Wilcoxon rank sum test for the comparison
#' dindex1 > dindex2
#' - dindex1: value of the first D index
#' - dindex2: value of the second D index
#'
#' @authors
#' Benjamin Haibe-Kains
#'
#' @references
#' Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008)
#'   "A comparative study of survival models for breast cancer prognostication
#'   based on microarray data: does a single gene beat them all?",
#'   Bioinformatics, 24, 19, pages 2200–2208.
#'
#' @seealso
#' [survcomp::D.index], [stats::t.test]
#'
#' @examples
#' set.seed(12345)
#' age <- rnorm(100, 50, 10)
#' size <- rexp(100,1)
#' stime <- rexp(100)
#' cens <- runif(100,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' d1 <- D.index(x=age, surv.time=stime, surv.event=sevent)
#' d2 <- D.index(x=size, surv.time=stime, surv.event=sevent)
#' dindex.comp(d1, d2)
#'
#' @md
#' @export
dindex.comp <-
function(dindex1, dindex2) {
	if(dindex1$n != dindex2$n) { stop("the D indices are computed from different number of samples!") }
	n <- dindex1$n
	x1 <- dindex1$data$z
	x2 <- dindex2$data$z
	beta1 <- dindex1$coef
	beta2 <- dindex2$coef
	se1 <- dindex1$se
	se2 <- dindex2$se
	r <- cor(x1, x2, method="spearman", use="complete.obs")
	if(abs(r) < 1) {
		t.stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "dindex1"=exp(beta1), "dindex2"=exp(beta2)))
}

