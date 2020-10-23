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

