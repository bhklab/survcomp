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


