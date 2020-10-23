#' @title Function to compare two concordance indices
#'
#' @description
#' This function compares two lists of concordance indices computed from the
#'   same survival data by using the function concordance.index. The statistical
#'   test is a Student t test for dependent samples.
#'
#' @param list.cindex1 first list of concordance indices as returned by the
#'   [survcomp::concordance.index] function.
#' @param list.cindex2 second list of concordance indices as returned by the
#'   [survcomp::concordance.index] function.
#' @param hetero if `TRUE`, a random effect model is use to compute the
#'   meta-estimators. Otherwise a fixed effect model is used.
#'
#' @details
#' In meta-analysis, we estimate the statistic of interest in several
#'   independent datasets. It results a list of estimates such as list of
#'   concordance indices. The two lists of concordance indices must be computed
#'   from the same samples (and corresponding survival data). The function
#'   computes a meta-estimator for the correlations between the two scores and
#'   uses a Student t test for dependent samples.
#'
#' @return
#' A list with items:
#' - p.value: p-value from the Student t test for the comparison cindex1 > cindex2.
#' - cindex1: meta-estimator of the first concordance index.
#' - cindex2: meta-estimator of the second concordance index.
#'
#'
#' @references
#' Cochrane, W. G. (1954) "The combination of estimates from different
#'   experiments", Biometrics, 10, pages 101–129.
#' Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008)
#'   "A comparative study of survival models for breast cancer prognostication
#'   based on microarray data: does a single gene beat them all?",
#'   Bioinformatics, 24, 19, pages 2200–2208.
#'
#' @seealso
#' [survcomp::concordance.index]
#'
#' @examples
#' #first dataset
#' set.seed(12345)
#' age <- rnorm(100, 50, 10)
#' size <- rexp(100,1)
#' stime <- rexp(100)
#' cens <- runif(100,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' c1.1 <- concordance.index(x=age, surv.time=stime, surv.event=sevent,
#'   method="noether")
#' c2.1 <- concordance.index(x=size, surv.time=stime, surv.event=sevent,
#'   method="noether")
#' #second dataset
#' set.seed(54321)
#' age <- rnorm(110, 53, 10)
#' size <- rexp(110,1.1)
#' stime <- rexp(110)
#' cens <- runif(110,.55,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' c1.2 <- concordance.index(x=age, surv.time=stime, surv.event=sevent,
#'   method="noether")
#' c2.2 <- concordance.index(x=size, surv.time=stime, surv.event=sevent,
#'   method="noether")
#' cindex.comp.meta(list.cindex1=list("cindex.age1"=c1.1, "cindex.age2"=c1.2),
#'   list.cindex2=list("cindex.size1"=c2.1, "cindex.size2"=c2.2))
#'
#' @md
#' @export
cindex.comp.meta <-
function(list.cindex1, list.cindex2, hetero=FALSE) {

	if(length(list.cindex1) != length(list.cindex2)) { stop("the number of concordance indices is not the same!") }
	eps <- 1E-15
	
	n <- 0
	x1 <- x1.se <- x2 <- x2.se <- corz <- corz.se <- NULL
	for(i in 1:length(list.cindex1)) {
		nn <- list.cindex1[[i]]$n
		if(nn != list.cindex2[[i]]$n) { stop("the number of samples to compute the concordance indices is not the same!") }
		if(nn > 3) {
			n <- n + nn
			x1 <- c(x1, list.cindex1[[i]]$c.index)
			x1.se <- c(x1.se, list.cindex1[[i]]$se)
			x2 <- c(x2, list.cindex2[[i]]$c.index)
			x2.se <- c(x2.se, list.cindex2[[i]]$se)
			cort <- cor(list.cindex1[[i]]$data$x, list.cindex2[[i]]$data$x, method="spearman", use="complete.obs")
			## since r is the spearman correlation coefficient and not the Pearson's one, we should apply a correction factor (see http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient for details)
			corz <- c(corz, sqrt((nn - 3) / 1.06) * fisherz(cort, inv=FALSE))
			corz.se <- c(corz.se, 1 / sqrt(nn - 3))
		} else {
			corz <- c(corz, NA)
			corz.se <- c(corz.se, NA)
		}
	}
	x1.meta <- combine.est(x=x1, x.se=x1.se, hetero=hetero, na.rm=TRUE)
	x2.meta <- combine.est(x=x2, x.se=x2.se, hetero=hetero, na.rm=TRUE)
	if(x1.meta$estimate == x2.meta$estimate && x1.meta$se == x2.meta$se) {
	## same concordance indices	
		return(list("p.value"=1, "cindex1"=x1.meta$estimate, "cindex2"=x2.meta$estimate))
	}
	rz <- combine.est(x=corz, x.se=corz.se, na.rm=TRUE, hetero=hetero)$estimate
	## since r is the spearman correlation coefficient and not the Pearson's one, we should apply a correction factor (see http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient for details)
	rz <- rz / (sqrt((n - 3) / 1.06))
	r <- fisherz(rz, inv=TRUE)
	if((1 - abs(r)) > eps) {
		t.stat <- (x1.meta$estimate - x2.meta$estimate) / sqrt(x1.meta$se^2 + x2.meta$se^2 - 2 * r * x1.meta$se * x2.meta$se)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "cindex1"=x1.meta$estimate, "cindex2"=x2.meta$estimate))
}
