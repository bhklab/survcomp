#' @name balanced.harzard.ratio
#'
#' @title Function to estimate the balanced hazard ratio through Cox regression
#'
#' @description
#' Function to compute the balanced hazard ratio for a risk group prediction.
#' The balanced hazard ratio is computed using the Cox model.
#'
#' @usage
#' balanced.hazard.ratio(x, surv.time, surv.event, alpha = 0.05,
#' method.test = c("logrank", "likelihood.ratio", "wald"),
#' ties = c("efron", "breslow", "exact"), weights, strat, ...)
#'
#' @param x a vector of risk group predictions.
#' @param surv.time a vector of event times.
#' @param surv.event a vector of event occurrence indicators.
#' @param alpha alpha level to compute confidence interval.
#' @param method.test ...
#' @param ties ...
#' @param weights ...
#' @param strat ...
#' @param ... Additional parmaeters to the hazard.ratio and coxph functions.
#'
#' @return
#' A list with the items:
#' - balanced.hazard.ratio: balanced hazard ratio estimate.
#' - coef: coefficient (beta) estimated in the cox regression model.
#' - se: standard error of the coefficient (beta) estimate.
#' - lower: lower bound for the confidence interval.
#' - upper: upper bound for the confidence interval.
#' - p.value: p-value computed using the score (logrank) test whether the
#' balanced hazard ratio is different from 1.
#' - n: number of samples used for the estimation.
#' - coxm: coxph.object fitted on the survival data and x.
#' - data: list of data used to compute the balanced hazard ratio (x,
#' surv.time and surv.event).
#'
#' @references
#' Branders, S. and Dupont, P. (2015) "A balanced hazard ratio for risk group
#'   evaluation from survival data", Statistics in Medicine, 34(17),
#'   pages 2528–2543.
#'
#' @seealso
#' [survcomp::hazard.ratio], [survival::coxph], [survival::coxph.object]
#'
#' @examples
#' set.seed(12345)
#' age <- rnorm(100, 50, 10)
#' stime <- rexp(100)
#' cens   <- runif(100,.5,2)
#' sevent  <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' strat <- sample(1:3, 100, replace=TRUE)
#' weight <- runif(100, min=0, max=1)
#' balanced.hazard.ratio(x=age, surv.time=stime, surv.event=sevent,
#'   weights=weight, strat=strat)
#'
#' @md
#' @export
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
