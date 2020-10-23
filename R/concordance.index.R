#' @title Function to compute the concordance index for survival or binary class prediction
#'
#' @details
#' Function to compute the concordance index for a risk prediction, i.e. the
#'   probability that, for a pair of randomly chosen comparable samples, the
#'   sample with the higher risk prediction will experience an event before
#'   the other sample or belongs to a higher binary class.
#'
#' @section Note:
#' The "direction" of the concordance index (< 0.5 or > 0.5) is the opposite
#'   than the [Hmisc::rcorr.cens] function in the Hmisc package. So you can
#'   easily get the same results than [Hmisc::rcorr.cens] by changing the sign
#'   of x.
#'
#' @param x a vector of risk predictions.
#' @param surv.time a vector of event times.
#' @param cl a vector of event occurence indicators.
#' @param weights weight of each sample.
#' @param comppairs threshold for comparable patients.
#' @param strat stratification indicator.
#' @param alpha apha level to compute confidence interval.
#' @param outx set to TRUE to not count pairs of observations tied on x as a
#'   relevant pair. This results in a Goodman-Kruskal gamma type rank
#'   correlation.
#' @param method can take the value conservative, noether or name
#'   (see paper Pencina et al. for details).
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of "two.sided" (default), "greater" (concordance index is
#'   greater than 0.5) or "less" (concordance index is less than 0.5). You can
#'   specify just the initial letter.
#' @param na.rm `TRUE` if missing values should be removed.
#'
#' @return A [`list`] containing the items:
#'   - c.index concordance index estimate.
#'   - se standard error of the estimate.
#'   - lower lower bound for the confidence interval.
#'   - upper upper bound for the confidence interval.
#'   - p.value p-value for the statistical test if the estimate if different from 0.5.
#'   - n number of samples used for the estimation.
#'   - data list of data used to compute the index (x, surv.time and surv.event, or cl).
#'   - comppairs number of compairable pairs.
#'
#' @references
#' Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996) "Tutorial in
#'   biostatistics: multivariable prognostic models: issues in developing
#'   models, evaluating assumptions and adequacy, and measuring and reducing
#'   error", Statistics in Medicine, 15, pages 361–387.
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
#' cat("survival prediction:\n")
#' concordance.index(x=age, surv.time=stime, surv.event=sevent, strat=strat,
#'   weights=weight, method="noether", comppairs=comppairs)
#' cat("binary class prediction:\n")
#' ## is age predictive of sex?
#' concordance.index(x=age, cl=sex, strat=strat, method="noether")
#'
#' @useDynLib
#' @md
#' @export concordance.index
concordance.index <-
function(x, surv.time, surv.event, cl, weights, comppairs=10, strat, alpha=0.05, outx=TRUE, method=c("conservative", "noether", "nam"), alternative=c("two.sided", "less", "greater"), na.rm=FALSE) {
	method <- match.arg(method)
	alternative <- match.arg(alternative)
	if(!missing(weights)) {
		if(length(weights) != length(x)) { stop("bad length for parameter weights!") }
		if(min(weights, na.rm=TRUE) < 0 && max(weights, na.rm=TRUE) > 1) { stop("weights must be a number between 0 and 1!") }
	} else { weights <- rep(1, length(x)) }
	if(!missing(strat)) {
		if(length(strat) != length(x)) { stop("bad length for parameter strat!") }
	} else { strat <- rep(1, length(x)) }

	if(missing(cl) && (missing(surv.time) || missing(surv.event))) { stop("binary classes and survival data are missing!") }
	if(!missing(cl) && (!missing(surv.time) || !missing(surv.event))) { stop("choose binary classes or survival data but not both!") }
	msurv <- FALSE
	if(missing(cl)) { ## survival data
		msurv <- TRUE
		cl <- rep(0, length(x))
	} else { surv.time <- surv.event <- rep(0, length(x)) } ## binary classes

	cc.ix <- complete.cases(x, surv.time, surv.event, cl, weights, strat)
	if(sum(cc.ix) < 3) {
	## not enough observations
		if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
		return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=sum(cc.ix), "data"=data, "comppairs"=NA))
	}
	if(any(!cc.ix) & !na.rm) { stop("NA values are present!") }
  ## remove samples whose the weight is equal to 0 to speed up the computation of the concordance index
	cc.ix <- cc.ix & weights != 0
	x2 <- x[cc.ix]
	cl2 <- cl[cc.ix]
	st <- surv.time[cc.ix]
	se <- surv.event[cc.ix]
	if(msurv && sum(se) == 0) {
		warning("\nno events, the concordance index cannot be computed!")
		data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event)
		return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=0, "data"=data, "comppairs"=NA))
	}
	if(!msurv && length(unique(cl2)) == 1) {
		warning("\nonly one class, the concordance index cannot be computed!")
		data  <- list("x"=x, "cl"=cl)
		return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=0, "data"=data, "comppairs"=NA))
	}
	weights <- weights[cc.ix]
	strat <- strat[cc.ix]
	strat <- as.numeric(as.factor(strat))
	ustrat <- sort(unique(strat)) ## to check later
	N <- sum(weights) ##length(x2)
	if(N <= 1) {
    warning("\nweights of observations are too small (sum should be > 1), the concordance index cannot be computed!")
    if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
    return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=length(x2), "data"=data, "comppairs"=NA))
  }

	ch <- dh <- uh <- rph <- rep(0, times=length(strat))
	lenS <- length(strat)
	lenU <- length(ustrat)
	out <- .C(.C_concordanceIndexC, as.integer(as.logical(msurv)), as.integer(ustrat), as.double(x2),
			as.integer(cl2), as.double(st), as.integer(se), as.double(weights), as.integer(strat),
			as.integer(N), as.integer(as.logical(outx)), ch = as.numeric(ch), dh = as.numeric(dh),
			uh = as.numeric(uh), rph = as.numeric(rph), as.integer(lenS), as.integer(lenU), PACKAGE="survcomp")
  ch <- out$ch
  dh <- out$dh
  uh <- out$uh
  rph <- out$rph
  cscount <- sum(ch + dh) ## comparable pairs
  if(sum(ch)==0 || sum(dh)==0 || sum(ch * (ch - 1))==0 || sum(dh * (dh - 1))==0 || sum(ch * dh)==0 || cscount < comppairs){
    if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
    return(list("c.index"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=length(x2), "data"=data, "comppairs"=cscount))
  }

	pc <- (1 / (N * (N - 1))) * sum(ch)
	pd  <- (1 / (N * (N - 1))) * sum(dh)
	cindex <- pc / (pc + pd)

	switch(method,
	"noether"={
    pcc <- (1 / (N * (N - 1) * (N - 2))) * sum(ch * (ch - 1))
    pdd <- (1 / (N * (N - 1) * (N - 2))) * sum(dh * (dh - 1))
    pcd <- (1 / (N * (N - 1) * (N - 2))) * sum(ch * dh)
    varp <- (4 / (pc + pd)^4) * (pd^2 * pcc - 2 * pc * pd * pcd + pc^2 * pdd)
    if((varp / N) > 0) {
      ci <- qnorm(p=alpha / 2, lower.tail=FALSE) * sqrt(varp / N)
      lower <- cindex - ci
      upper <- cindex + ci
      switch(alternative,
		   "two.sided"={ p <- pnorm((cindex - 0.5) / sqrt(varp / N), lower.tail=cindex < 0.5) * 2 },
		   "less"={ p <- pnorm((cindex - 0.5) / sqrt(varp / N), lower.tail=TRUE) },
		   "greater"={  p <- pnorm((cindex - 0.5) / sqrt(varp / N), lower.tail=FALSE) }
	    )
    } else { ci <- lower <- upper <- p <- NA }
  },
	"conservative"={
    C <- cindex
    ## pc and pd have been computed previously
    w <- (2 * qnorm(p=alpha / 2, lower.tail=FALSE)^2) / (N * (pc + pd))
    ci <- sqrt(w^2 + 4 * w * C * (1 - C)) / (2 * (1 + w))
    point <- (w + 2 * C) / (2 * (1 + w))
    lower <- point - ci
    upper <- point + ci
    cindex <- C
    p <- NA
    varp <- NA
  },
  "name"={
    stop("method not implemented!")
  })
  #bound the confidence interval
  lower <- ifelse(lower < 0, 0, lower)
  lower <- ifelse(lower > 1, 1, lower)
  upper <- ifelse(upper < 0, 0, upper)
  upper <- ifelse(upper > 1, 1, upper)
  if(msurv) { data <- list("x"=x, "surv.time"=surv.time, "surv.event"=surv.event) } else { data  <- list("x"=x, "cl"=cl) }
  if(!is.na(varp) && (varp / N) > 0) { se <- sqrt(varp / N) } else { se <- NA }
  return(list("c.index"=cindex, "se"=se, "lower"=lower, "upper"=upper, "p.value"=p, "n"=length(x2), "data"=data, "comppairs"=cscount))
}
