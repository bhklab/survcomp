`no.at.risk` <- 
function( formula.s, data.s, sub.s="all", t.step, t.end ) {
# Updated 6.6.11 to work from summary.surfvit

    require(survival)

    if( length(sub.s)==1 && sub.s=="all" ) sub.s <- rep(TRUE, nrow(data.s))
    assign("sub.s", sub.s, envir=.GlobalEnv)

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