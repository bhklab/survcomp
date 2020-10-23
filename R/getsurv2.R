#' @name getsurv2
#'
#' @title Function to retrieve the survival probabilities at a specific point
#'   in time
#'
#' @description
#' The function retrieves the survival probabilities from a survfit object, for
#'   a specific point in time.
#'
#' @usage getsurv2(sf, time, which.est = c("point", "lower", "upper"))
#'
#' @param sf survfit object
#' @param time time at which the survival probabilities must be retrieved
#' @param which.est which estimation to be returned? point for the point
#'   estimate, lower for the lower bound and upper for the upper bound
#'
#' @details
#' The survival probabilities are estimated through the [survival::survfit]
#'   function.
#'
#' @return vector of survival probabilities
#'
#' @authors
#' Benjamin Haibe-Kains
#'
#' @seealso
#' [survival::survfit]
#'
#' @exampless
#' set.seed(12345)
#' age <- rnorm(30, 50, 10)
#' stime <- rexp(30)
#' cens <- runif(30,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' sf <- survfit(Surv(stime, sevent) ~ 1)
#' getsurv2(sf, time=1)
#'
#' @md
#' @export
getsurv2 <-
function(sf, time, which.est=c("point", "lower", "upper")) {
	which.est <- match.arg(which.est)
	switch(which.est,
	"lower"={
		survres <- as.matrix(sf$lower)	
	},
	"point"={
		survres <- as.matrix(sf$surv)
	},
	"upper"={
		survres <- as.matrix(sf$upper)
	})
	if(time >= max(sf$time)) {
		time.ix <- length(sf$time)
	} else { time.ix <- order(sf$time <= time)[1] - 1 }
	if(time.ix == 0) {
		res <- rep(1, ncol(sf$surv[time.ix, ]))
		names(res) <- dimnames(sf$surv)[[2]]
	} else { res <- survres[time.ix, ] }	
	return(res)
}

