#' @title Function to artificially censor survival data
#'
#' @description
#' The function censors the survival data at a specific point in time. This is
#'   is useful if you used datasets having different follow-up periods.
#'
#' @usage
#' censor.time(surv.time, surv.event, time.cens = 0)
#'
#' @param surv.time vector of times to event occurrence
#' @param surv.event vector of indicators for event occurrence
#' @param time.cens point in time at which the survival data must be censored
#'
#' @return
#' A list with the items:
#' - surv.time.cens: vector of censored times to event occurrence
#' - srv.event.cens: vector of censored indicators for event occurrence
#'
#'
#' @examples
#' set.seed(12345)
#' stime <- rexp(30)
#' cens <- runif(30,0.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' censor.time(surv.time=stime, surv.event=sevent, time.cens=1)
#'
#' @md
#' @export
censor.time <-
function(surv.time, surv.event, time.cens=0) {
	stc <- surv.time
   sec <- surv.event
   cc.ix <- complete.cases(stc, sec)
   if(time.cens != 0) { 
   	stc[cc.ix][surv.time[cc.ix] > time.cens] <- time.cens
   	sec[cc.ix][surv.time[cc.ix] > time.cens] <- 0
   }
   return(list("surv.time.cens"=stc, "surv.event.cens"=sec))
}

