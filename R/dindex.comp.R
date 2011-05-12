`dindex.comp` <-
function(dindex1, dindex2) {
	if(dindex1$n != dindex2$n) { stop("the D indices are computed from different number of samples!") }
	n <- dindex1$n
	x1 <- dindex1$data$z
	x2 <- dindex2$data$z
	beta1 <- dindex1$coef
	beta2 <- dindex2$coef
	se1 <- dindex1$se
	se2 <- dindex2$se
	r <- cor(x1, x2, method="spearman", use="complete.obs")
	if(abs(r) < 1) {
		t.stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "dindex1"=exp(beta1), "dindex2"=exp(beta2)))
}

