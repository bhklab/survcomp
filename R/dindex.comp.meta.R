#' @name dindex.comp.meta
#'
#' @title Function to compare two D indices
#'
#' @description
#' This function compares two lists of D indices computed from the same
#'   survival data by using the function D.index. The statistical test is a
#'   Student t test for dependent samples.
#'
#' @param list.dindex1 first list of D indices as returned by the `D.index`
#'   function.
#' @param list.dindex2 second list of D indices as returned by the `D.index`
#'   function.
#' @param hetero if `TRUE`, a random effect model is use to compute the
#'   meta-estimators. Otherwise a fixed effect model is used.
#'
#' @details
#' In meta-analysis, we estimate the statistic of interest in several
#'   independent datasets. It results a list of estimates such as list of D
#'   indices. The two lists of D indices must be computed from the same samples
#'   (and corresponding survival data). The function computes a meta-estimator
#'   for the correlations between the two scores and uses a Student t test for
#'   dependent samples.
#'
#' @return
#' A list with items:
#' - p.value: p-value from the Student t test for the comparison dindex1 >
#' dindex2.
#' - dindex1: meta-estimator of the first D index.
#' - dindex2: meta-estimator of the second D index.
#'
#' @references
#' Cochrane, W. G. (1954) "The combination of estimates from different
#'   experiments", Biometrics, 10, pages 101–129.
#'
#' Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A
#'   comparative study of survival models for breast cancer prognostication
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
#' d1.1 <- D.index(x=age, surv.time=stime, surv.event=sevent)
#' d2.1 <- D.index(x=size, surv.time=stime, surv.event=sevent)
#' #second dataset
#' set.seed(54321)
#' age <- rnorm(110, 53, 10)
#' size <- rexp(110,1.1)
#' stime <- rexp(110)
#' cens <- runif(110,.55,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' d1.2 <- D.index(x=age, surv.time=stime, surv.event=sevent)
#' d2.2 <- D.index(x=size, surv.time=stime, surv.event=sevent)
#' dindex.comp.meta(list.dindex1=list("dindex.age1"=d1.1, "dindex.age2"=d1.2),
#'   list.dindex2=list("dindex.size1"=d2.1, "dindex.size2"=d2.2))
#'
#' @md
#' @export
dindex.comp.meta <-
function(list.dindex1, list.dindex2, hetero=FALSE) {

	if(length(list.dindex1) != length(list.dindex2)) { stop("the concordance indices are computed from different number of samples!") }

	n <- 0
	x1 <- x1.se <- x2 <- x2.se <- corz <- corz.se <- NULL
	for(i in 1:length(list.dindex1)) {
		nn <- list.dindex1[[i]]$n
		if(nn != list.dindex2[[i]]$n) { stop("the number of samples to compute the concordance indices is not the same!") }
		n <- n + nn
		x1 <- c(x1, list.dindex1[[i]]$coef)
		x1.se <- c(x1.se, list.dindex1[[i]]$se)
		x2 <- c(x2, list.dindex2[[i]]$coef)
		x2.se <- c(x2.se, list.dindex2[[i]]$se)
		cort <- cor(list.dindex1[[i]]$data$z, list.dindex2[[i]]$data$z, method="spearman", use="complete.obs")
		## since r is the spearman correlation coefficient and not the Pearson's one, we should apply a correction factor (see http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient for details)
		corz <- c(corz, sqrt((nn - 3) / 1.06) * fisherz(cort, inv=FALSE))
		if(nn > 3) { corz.se <- c(corz.se, 1 / sqrt(nn - 3)) } else { corz.se <- c(corz.se, NA) }
	}
	x1.meta <- combine.est(x=x1, x.se=x1.se, hetero=hetero, na.rm=TRUE)
	x2.meta <- combine.est(x=x2, x.se=x2.se, hetero=hetero, na.rm=TRUE) 
	if(x1.meta$estimate == x2.meta$estimate && x1.meta$se == x2.meta$se) {
	## same D indices	
		return(list("p.value"=1, "cindex1"=x1.meta$estimate, "cindex2"=x2.meta$estimate))
	}
	rz <- combine.est(x=corz, x.se=corz.se, na.rm=TRUE, hetero=hetero)$estimate
	## since r is the spearman correlation coefficient and not the Pearson's one, we should apply a correction factor (see http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient for details)
	rz <- rz / (sqrt((n - 3) / 1.06))
	r <- fisherz(rz, inv=TRUE)

	if(abs(r) < 1) {
		t.stat <- (x1.meta$estimate - x2.meta$estimate) / sqrt(x1.meta$se^2 + x2.meta$se^2 - 2 * r * x1.meta$se * x2.meta$se)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "dindex1"=exp(x1.meta$estimate), "dindex2"=exp(x2.meta$estimate)))
}
