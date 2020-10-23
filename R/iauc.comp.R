iauc.comp <-
function(auc1, auc2, time) {
	if((length(auc1) + length(auc2) + length(time)) != 3 * length(time)) { stop("auc1, auc2 and time must have the same length!") }
	cc.ix <- complete.cases(auc1, auc2, time)
	auc1 <- auc1[cc.ix]
	auc2 <- auc2[cc.ix]
	time <- time[cc.ix]
	diffs <- c(time[1], time[2:length(time)] - time[1:(length(time) - 1)])
	iauc1 <- sum(diffs * auc1) / max(time)
	iauc2 <- sum(diffs * auc2) / max(time)
	rr <- wilcox.test(x=auc1, y=auc2, alternative="greater", paired=TRUE, exact=FALSE)
	return(list("p.value"=rr$p.value, "iauc1"=iauc1, "iauc2"=iauc2))
}

