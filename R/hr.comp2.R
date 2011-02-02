`hr.comp2` <-
function(x1, beta1, se1, x2, beta2, se2, n) {
	r <- cor(x1, x2, use="complete.obs")
	if(abs(r) < 1) {
		t.stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "hr1"=exp(beta1), "hr2"=exp(beta2)))
}

