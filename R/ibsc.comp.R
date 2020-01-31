`ibsc.comp` <-
function(bsc1, bsc2, time) {
	if((length(bsc1) + length(bsc2) + length(time)) != 3 * length(time)) { stop("bsc1, bsc2 and time must have the same length!") }
	cc.ix <- complete.cases(bsc1, bsc2, time) & !duplicated(time)
	bsc1 <- bsc1[cc.ix]
	bsc2 <- bsc2[cc.ix]
	time <- time[cc.ix]
	diffs <- c(time[1], time[2:length(time)] - time[1:(length(time) - 1)])
	ibsc1 <- sum(diffs * bsc1) / max(time)
	ibsc2 <- sum(diffs * bsc2) / max(time)
	rr <- wilcox.test(x=bsc1, y=bsc2, alternative="less", paired=TRUE, exact=FALSE)
	return(list("p.value"=rr$p.value, "ibsc1"=ibsc1, "ibsc2"=ibsc2))
}

