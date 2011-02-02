`cindex.comp.meta` <-
function(list.cindex1, list.cindex2, hetero=FALSE) {

	if(length(list.cindex1) != length(list.cindex2)) { stop("the number of concordance indices is not the same!") }
	eps <- 1E-15
	
	n <- 0
	x1 <- x1.se <- x2 <- x2.se <- corz <- corz.se <- NULL
	for(i in 1:length(list.cindex1)) {
		nn <- list.cindex1[[i]]$n
		if(nn != list.cindex2[[i]]$n) { stop("the number of samples to compute the concordance indices is not the same!") }
		n <- n + nn
		x1 <- c(x1, list.cindex1[[i]]$c.index)
		x1.se <- c(x1.se, list.cindex1[[i]]$se)
		x2 <- c(x2, list.cindex2[[i]]$c.index)
		x2.se <- c(x2.se, list.cindex2[[i]]$se)
		cort <- cor(list.cindex1[[i]]$data$x, list.cindex2[[i]]$data$x, method="pearson", use="complete.obs")
		corz <- c(corz, fisherz(cort, inv=FALSE))
		if(nn > 3) { corz.se <- c(corz.se, 1 / sqrt(nn - 3)) } else { corz.se <- c(corz.se, NA) }
	}
	x1.meta <- combine.est(x=x1, x.se=x1.se, hetero=hetero, na.rm=TRUE)
	x2.meta <- combine.est(x=x2, x.se=x2.se, hetero=hetero, na.rm=TRUE) 
	r <- fisherz(combine.est(x=corz, x.se=corz.se, na.rm=TRUE, hetero=hetero)$estimate, inv=TRUE)

	if((1 - abs(r)) > eps) {
		t.stat <- (x1.meta$estimate - x2.meta$estimate) / sqrt(x1.meta$se^2 + x2.meta$se^2 - 2 * r * x1.meta$se * x2.meta$se)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "cindex1"=x1.meta$estimate, "cindex2"=x2.meta$estimate))
}
