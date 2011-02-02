`hr.comp.meta` <-
function(list.hr1, list.hr2, hetero=FALSE) {

	if(length(list.hr1) != length(list.hr2)) { stop("the concordance indices are computed from different number of samples!") }

	n <- 0
	x1 <- x1.se <- x2 <- x2.se <- corz <- corz.se <- NULL
	for(i in 1:length(list.hr1)) {
		nn <- list.hr1[[i]]$n
		if(nn != list.hr2[[i]]$n) { stop("the number of samples to compute the concordance indices is not the same!") }
		n <- n + nn
		x1 <- c(x1, list.hr1[[i]]$coef)
		x1.se <- c(x1.se, list.hr1[[i]]$se)
		x2 <- c(x2, list.hr2[[i]]$coef)
		x2.se <- c(x2.se, list.hr2[[i]]$se)
		cort <- cor(list.hr1[[i]]$data$x, list.hr2[[i]]$data$x, method="pearson", use="complete.obs")
		corz <- c(corz, fisherz(cort, inv=FALSE))
		if(nn > 3) { corz.se <- c(corz.se, 1 / sqrt(nn - 3)) } else { corz.se <- c(corz.se, NA) }
	}
	x1.meta <- combine.est(x=x1, x.se=x1.se, hetero=hetero, na.rm=TRUE)
	x2.meta <- combine.est(x=x2, x.se=x2.se, hetero=hetero, na.rm=TRUE) 
	r <- fisherz(combine.est(x=corz, x.se=corz.se, na.rm=TRUE, hetero=hetero)$estimate, inv=TRUE)

	if(abs(r) < 1) {
		t.stat <- (x1.meta$estimate - x2.meta$estimate) / sqrt(x1.meta$se^2 + x2.meta$se^2 - 2 * r * x1.meta$se * x2.meta$se)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "hr1"=exp(x1.meta$estimate), "hr2"=exp(x2.meta$estimate)))
}
