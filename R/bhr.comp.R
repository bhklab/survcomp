`bhr.comp` <-
function(bhr1, bhr2) {
	if(bhr1$n != bhr2$n) { stop("the balanced hazard ratios are computed from different number of samples!") }
	n <- bhr1$n
	x1 <- bhr1$data$x
	x2 <- bhr2$data$x
	beta1 <- bhr1$coef
	beta2 <- bhr2$coef
	se1 <- bhr1$se
	se2 <- bhr2$se
	r <- cor(x1, x2, method="spearman", use="complete.obs")
	if(abs(r) < 1) {
		t.stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "bhr1"=exp(beta1), "bhr2"=exp(beta2)))
}

