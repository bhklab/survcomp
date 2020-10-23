#' @title Function to compute the number of individuals at risk
#'
#' @description
#' Function to compute the number of individuals at risk at certain time points,
#'   as used in the Kaplan-Meier estimator for instance, depending on
#'   stratification.
#'
#' @param formula.s formula composed of a Surv object and a strata variable
#'   (i.e. stratification).
#' @param data.s data frame composed of the variables used in the formula.
#' @param sub.s vector of booleans specifying if only a subset of the data
#'   should be considered.
#' @param t.step time step at which the number of individuals at risk is
#'   computed.
#' @param t.end maximum time to be considered.
#'
#' @details
#' The original version of this function was kindly provided by Dr Christos
#'   Hatzis (January, 17th 2006).
#'
#' @return
#' number of individuals at risk at each time step specified in t.step up to
#'   t.end.
#'
#' @seealso
#' [survival::survfit], [survcomp::km.coxph.plot]
#'
#' @examples
#' set.seed(12345)
#' stime <- rexp(100)
#' cens   <- runif(100,.5,2)
#' sevent  <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' strat <- sample(1:3, 100, replace=TRUE)
#' dd <- data.frame("surv.time"=stime, "surv.event"=sevent, "strat"=strat)
#' no.at.risk(formula.s=Surv(surv.time,surv.event) ~ strat, data.s=dd,
#'   sub.s="all", t.step=0.05, t.end=1)
#'
#' @md
#' @export
no.at.risk <-
function( formula.s, data.s, sub.s="all", t.step, t.end ) {
# Updated 6.6.11 to work from summary.surfvit

    if( length(sub.s)==1 && sub.s=="all" ) sub.s <- rep(TRUE, nrow(data.s))
    pos <- 1
    envir = as.environment(pos)
    assign("sub.s", sub.s, envir = envir)

    sf <- survfit( formula.s, data=data.s, subset=sub.s )
    if (is.null(sf$strata))
        sf$strata <- c("All" = length(sf$time))
    n.strata <- length(sf$strata)

    t.pts <- seq(0, t.end, t.step)
    sumsf <- summary(sf, times = t.pts, extend = TRUE)
    tms <- with(sumsf, split(time, strata))
    rsk <- with(sumsf, split(n.risk, strata))

    nar <- do.call("rbind", rsk)
    nar <- data.frame(names(sf$strata), nar)
    colnames(nar) <- c("risk.factor", as.character(tms[[1]]))

    remove("sub.s", envir=.GlobalEnv)
    nar
}
