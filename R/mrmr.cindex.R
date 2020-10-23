#' @title
#' Function to compute the concordance index for survival or binary class
#'   prediction
#'
#' @description
#' Function to compute the minimum redundancy - maximum relevance (mRMR)
#'   ranking for a risk prediction or a binary classification task based on the
#'   concordance index. The mRMR feature selection has been adapted to use the
#'   concordance index to estimate the correlation between a variable and the
#'   output (binary or survival) data.
#'
#' @usage
#' mrmr.cindex(x, surv.time, surv.event, cl, weights, comppairs=10, strat,
#' alpha = 0.05, outx = TRUE, method = c("conservative", "noether", "nam"),
#' alternative = c("two.sided", "less", "greater"), na.rm = FALSE)
#'
#' @param x a vector of risk predictions.
#' @param surv.time a vector of event times.
#' @param surv.event a vector of event occurence indicators.
#' @param cl a vector of binary class indicators.
#' @param weights weight of each sample.
#' @param comppairs threshold for compairable patients.
#' @param strat stratification indicator.
#' @param alpha apha level to compute confidence interval.
#' @param outx set to `TRUE` to not count pairs of observations tied on x as a
#'   relevant pair. This results in a Goodman-Kruskal gamma type rank
#'   correlation.
#' @param method can take the value conservative, noether or name (see paper
#'   Pencina et al. for details).
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of "two.sided" (default), "greater" (concordance index is
#'   greater than 0.5) or "less" (concordance index is less than 0.5). You can
#'   specify just the initial letter.
#' @param na.rm TRUE if missing values should be removed.
#'
#' @return A mRMR ranking
#'
#' @section Note:
#' The "direction" of the concordance index (< 0.5 or > 0.5) is the opposite
#'   than the [Hmisc::rcorr.cens] function in the Hmisc package. So you can
#'   easily get the same results than [Hmisc::rcorr.cens] by changing the sign
#'   of x.
#'
#' @references
#' Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996) "Tutorial in
#'   biostatistics: multivariable prognostic models: issues in developing
#'   models, evaluating assumptions and adequacy, and measuring and reducing
#'   error", Statistics in Medicine, 15, pages 361–387.
#'
#' Pencina, M. J. and D'Agostino, R. B. (2004) "Overall C as a measure of
#'   discrimination in survival analysis: model specific population value and
#'   confidence interval estimation", Statistics in Medicine, 23, pages
#'   2109–2123, 2004.
#'
#' @seealso
#' [Hmisc::rcorr.cens], [CPE::phcpe], [clinfun::coxphCPE]
#'
#' @examples
#' set.seed(12345)
#' age <- rnorm(100, 50, 10)
#' sex <- sample(0:1, 100, replace=TRUE)
#' stime <- rexp(100)
#' cens   <- runif(100,.5,2)
#' sevent  <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' strat <- sample(1:3, 100, replace=TRUE)
#' weight <- runif(100, min=0, max=1)
#' comppairs <- 10
#' xx <- data.frame("age"=age, "sex"=sex)
#' cat("survival prediction:\n")
#' mrmr.cindex(x=xx, surv.time=stime, surv.event=sevent, strat=strat, weights=weight,
#' method="noether", comppairs=comppairs)
#'
#' @md
#' @export
mrmr.cindex <-
function(x, surv.time, surv.event, cl, weights, comppairs=10, strat, alpha=0.05, outx=TRUE, method=c("conservative", "noether", "nam"), alternative=c("two.sided", "less", "greater"), na.rm=FALSE) {

	nvar<-ncol(x)
	nsample<-nrow(x)
	if(!missing(weights)) {
		if(length(weights) != nsample) { stop("bad length for parameter weights!") }
		if(min(weights, na.rm=TRUE) < 0 && max(weights, na.rm=TRUE) > 1) { stop("weights must be a number between 0 and 1!") }
	} else { weights <- rep(1, nsample) }
	if(!missing(strat)) {
		if(length(strat) != nsample) { stop("bad length for parameter strat!") }
	} else { strat <- rep(1, nsample) }

	if(missing(cl) && (missing(surv.time) || missing(surv.event))) { stop("binary classes and survival data are missing!") }
	if(!missing(cl) && (!missing(surv.time) || !missing(surv.event))) { stop("choose binary classes or survival data but not both!") }


	res_cIndex<-rep(0,nvar)

		msurv <- FALSE
		if(missing(cl)) { ## survival data
			msurv <- TRUE
			cl <- rep(0, nsample)
		} else { surv.time <- surv.event <- rep(0, nsample) } ## binary classes


	#### get cIndex for each variable in dataset with surv.time, surv.event ####
	for(ivar in 1:nvar){
		is.correct <- TRUE
		xx <- x[ ,ivar]
		cc.ix <- complete.cases(xx, surv.time, surv.event, cl, weights, strat)

		if(sum(cc.ix) < 3) {
			## not enough observations
			if(msurv) { data <- list("x"=xx, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=xx, "cl"=cl) }
			return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=sum(cc.ix), "data"=data, "comppairs"=NA))
		}

		if(any(!cc.ix) & !na.rm) { stop("NA values are present!") }

		# remove samples whose the weight is equal to 0 to speed up the computation of the concordance index
		cc.ix <- cc.ix & weights != 0
		x2 <- xx[cc.ix]
		cl2 <- cl[cc.ix]
		st <- surv.time[cc.ix]
		se <- surv.event[cc.ix]
		if(msurv && sum(se) == 0) {
			warning(paste("\nno events, the concordance index cannot be computed for variable ", colnames(x)[ivar]," !"))
			res_cIndex[ivar]<-NA
			is.correct<-FALSE
		}
		if(!msurv && length(unique(cl2)) == 1) {
			warning(paste("\nonly one class, the concordance index cannot be computed for variable", colnames(x)[ivar]," !"))
			res_cIndex[ivar]<-NA
			is.correct<-FALSE
		}
		weights <- weights[cc.ix]
		strat <- strat[cc.ix]
		strat <- as.numeric(as.factor(strat))
		ustrat <- sort(unique(strat)) ## to check later
		N <- sum(weights) ##length(x2)
		if(N <= 1) {
			warning(paste("\nweights of observations are too small (sum should be > 1), the concordance index cannot be computed for variable", colnames(x)[ivar]," !"))
			res_cIndex[ivar]<-NA
			is.correct<-FALSE
		}

		lenS <- length(strat)
		lenU <- length(ustrat)
		if(is.correct){
			res_cIndex[ivar] <- .Call(.C_get_concordanceIndex_onevariable, as.integer(as.logical(msurv)), as.integer(ustrat), as.double(x2),as.integer(cl2), as.double(st), as.integer(se), as.double(weights), as.integer(strat),as.integer(N), as.integer(as.logical(outx)), as.integer(lenS), as.integer(lenU))
		}
	}

	res_mrmr_cIndex <- .Call(.C_mrmr_cIndex, data.matrix(x),as.integer(is.na(x)),as.double(res_cIndex),ncol(x), nrow(x), -1000)
	names(res_mrmr_cIndex)<-colnames(x)
	return(res_mrmr_cIndex)

}


