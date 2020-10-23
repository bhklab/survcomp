#' @name cindex.comp
#'
#' @title Function to compare two concordance indices
#'
#' @description
#' This function compares two concordance indices computed from the same data by
#'   using the function concordance.index. The statistical test is a Student
#'   t-test for dependent samples.
#'
#' @usage cindex.comp(cindex1, cindex2)
#'
#' @param cindex1 first concordance index as returned by the concordance.index function.
#' @param cindex2 second concordance index as returned by the concordance.index function.
#'
#' @return
#' A list with items:
#' - p.value: p-value from the Student t test for the comparison cindex1 >
#' - cindex1: value of the first concordance index.
#' - cindex2: value of the second concordance index.
#'
#' @references
#' Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008)
#'   "A comparative study of survival models for breast cancer prognostication
#'   based on microarray data: does a single gene beat them all?",
#'   Bioinformatics, 24(19), pages 2200–2208.
#'
#' @seealso
#' [survcomp::concordance.index]
#'
#' @examples
#' set.seed(12345)
#' age <- rnorm(100, 50, 10)
#' size <- rexp(100,1)
#' stime <- rexp(100)
#' cens <- runif(100,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' c1 <- concordance.index(x=age, surv.time=stime, surv.event=sevent,
#'   method="noether")
#' c2 <- concordance.index(x=size, surv.time=stime, surv.event=sevent,
#'   method="noether")
#' cindex.comp(c1, c2)
#'
#' @md
#' @export
cindex.comp <-
function(cindex1, cindex2) {

	if(cindex1$n != cindex2$n) { stop("the concordance indices are computed from different number of samples!") }
	if(is.na(cindex1$se) || is.na(cindex2$se)){stop("the concordance indices must be computed using method noether!")}
	eps <- 1E-15
	
	n <- cindex1$n
	r <- cor(cindex1$data$x, cindex2$data$x, use="complete.obs", method="spearman")
	if((1 - abs(r)) > eps) {
		t.stat <- (cindex1$c.index - cindex2$c.index) / sqrt(cindex1$se^2 + cindex2$se^2 - 2 * r * cindex1$se * cindex2$se)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "cindex1"=cindex1$c.index, "cindex2"=cindex2$c.index))
}

