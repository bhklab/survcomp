#' @title Function to compare two IAUCs through time-dependent ROC curves
#'
#' @description
#' This function compares two integrated areas under the curves (IAUC) through
#'   the results of time-dependent ROC curves at some points in time. The
#'   statistical test is a Wilcoxon rank sum test fordependent samples.
#'
#' @usage iauc.comp(auc1, auc2, time)
#'
#' @param auc1 vector of AUCs computed from the first time-dependent ROC curves
#'   for somepoints in time
#' @param auc2 c2vector  of  AUCs  computed  from  the  second  time-dependent
#'   ROC  curves  forsome points in time
#' @param time vector of points in time for which the AUCs are computed
#'
#' @details
#' The two vectors of AUCs must be computed from the same samples (and
#'   corresponding survivaldata) and for the same points in time. The function
#'   uses a Wilcoxon rank sum test for dependent samples.
#'
#' @return
#' A list with items:
#' - p.value: p-value from the Wilcoxon rank sum test for the comparison
#' iauc1 > iauc2
#' - iauc1: value of the IAUC for the first set of time-dependent ROC curves
#' - iauc2: value of the IAUC for the second set of time-dependent ROC curves
#'
#' @references
#' Wilcoxon, F. (1945) "Individual comparisons by ranking methods", Biometrics
#'   Bulletin,1, pages80–83
#'
#' Haibe-Kains, B. and Desmedt, C. and Sotiriou, C. and Bontempi, G. (2008) "A
#'   comparative studyof survival models for breast cancer prognostication
#'   based on microarray data:  does a single genebeat them all?",
#'   Bioinformatics, 24, 19, pages 2200–2208.
#'
#' @seealso
#' [survcomp::tdrocc], [stats::wilcox.test]
#'
#' @examples
#' set.seed(12345)
#' age <- rnorm(30, 50, 10)
#' size <- rexp(30,1)
#' stime <- rexp(30)
#' cens <- runif(30,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' ##time-dependent ROC curves
#' tt <- unique(sort(stime[sevent == 1]))
#' ##size
#' mytdroc1 <- NULL
#' for(i in 1:length(tt)) {
#' 	  rr <- tdrocc(x=size, surv.time=stime, surv.event=sevent, time=tt[i], na.rm=TRUE, verbose=FALSE)
#' 	  mytdroc1 <- c(mytdroc1, list(rr))
#' }
#' auc1 <- unlist(lapply(mytdroc1, function(x) { return(x$AUC) }))
#' ##age
#' mytdroc2 <- NULL
#' for (i in 1:length(tt)) {
#' 	  rr <- tdrocc(x=age, surv.time=stime, surv.event=sevent, time=tt[i], na.rm=TRUE, verbose=FALSE)
#'    mytdroc2 <- c(mytdroc2, list(rr))
#' }
#' auc2 <- unlist(lapply(mytdroc2, function(x) { return(x$AUC) }))
#' iauc.comp(auc1=auc1, auc2=auc2, time=tt)
#'
#' @md
#' @export
iauc.comp <-
function(auc1, auc2, time) {
	if((length(auc1) + length(auc2) + length(time)) != 3 * length(time)) { stop("auc1, auc2 and time must have the same length!") }
	cc.ix <- complete.cases(auc1, auc2, time)
	auc1 <- auc1[cc.ix]
	auc2 <- auc2[cc.ix]
	time <- time[cc.ix]
	diffs <- c(time[1], time[2:length(time)] - time[1:(length(time) - 1)])
	iauc1 <- sum(diffs * auc1) / max(time)
	iauc2 <- sum(diffs * auc2) / max(time)
	rr <- wilcox.test(x=auc1, y=auc2, alternative="greater", paired=TRUE, exact=FALSE)
	return(list("p.value"=rr$p.value, "iauc1"=iauc1, "iauc2"=iauc2))
}

