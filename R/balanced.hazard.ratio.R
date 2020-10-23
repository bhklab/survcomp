balanced.hazard.ratio <-
function(x, surv.time, surv.event, alpha=0.05, method.test = c("logrank", "likelihood.ratio", "wald"), ties=c("efron","breslow","exact"), weights, strat, ...)
{
    #Balanced Hazard ratio
    
    if (missing(method.test))
    {
        method.test = "logrank"
    }
    if (missing(ties))
    {
        ties = "breslow"
    }
    if(!missing(weights))
    {
        if(length(weights) != length(x))
        {
            stop("bad length for parameter weights!")
        }
    } else {
        weights <- rep(1,  length(x))
    }
    if(!missing(strat)) {
        if(length(strat) != length(x))
        {
            stop("bad length for parameter strat!")
        }
        ## remove weights=0 because the coxph function does not deal with them properly
        iix <- weights <= 0
        if(any(iix)) { warning("samples with weight<=0 are discarded") }
        weights[iix] <- NA
    } else { 
        strat <- rep(1,  length(x))
    }
    
    
    #Duplicating the patients
    cl = sort(unique(x))
    
    repTime = c()
    repEvent = c()
    repX = c()
    repWeights = c()
    repStrat = c()
    
    xOld = x
    for(i in 1:(length(cl)))
    {
        x[xOld==cl[i]] = i*2
    }
    
    for(i in 1:(length(cl)-1))
    {
        ind = c(which(x==i*2),which(x==(i+1)*2))
        repTime = c(repTime,surv.time[ind])
        repEvent = c(repEvent,surv.event[ind])
        repX = c(repX,rep(i*2+1,length(ind)))
        repWeights = c(repWeights,weights[ind])
        repStrat = c(repStrat,strat[ind])
    }
    
    #Compute the hazard ratio on the duplicated patients
    BHr = hazard.ratio( x=c(x,repX), surv.time=c(surv.time,repTime), surv.event=c(surv.event,repEvent), weights=c(weights,repWeights), strat=c(strat,repStrat), alpha=alpha, method.test=method.test, ties="breslow", ...)
    
    BHr$balanced.hazard.ratio = BHr$hazard.ratio
    BHr$hazard.ratio = NULL
    BHr$n = length(x)
    BHr$data$x = x
    BHr$data$surv.time = surv.time
    BHr$data$surv.event = surv.event
    
    return(BHr)
}
