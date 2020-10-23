hr.comp <-
function(hr1, hr2) {
	if(hr1$n != hr2$n) { stop("the hazard ratios are computed from different number of samples!") }
	n <- hr1$n
	x1 <- hr1$data$x
	x2 <- hr2$data$x
	beta1 <- hr1$coef
	beta2 <- hr2$coef
	se1 <- hr1$se
	se2 <- hr2$se
	r <- cor(x1, x2, method="spearman", use="complete.obs")
	if(abs(r) < 1) {
		t.stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "hr1"=exp(beta1), "hr2"=exp(beta2)))
}

