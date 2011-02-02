`cindex.comp` <-
function(cindex1, cindex2) {

	if(cindex1$n != cindex2$n) { stop("the concordance indices are computed from different number of samples!") }
	if(is.na(cindex1$se) || is.na(cindex2$se)){stop("the concordance indices must be computed using method noether!")}
	eps <- 1E-15
	
	n <- cindex1$n
	r <- cor(cindex1$data$x, cindex2$data$x, use="complete.obs")
	if((1 - abs(r)) > eps) {
		t.stat <- (cindex1$c.index - cindex2$c.index) / sqrt(cindex1$se^2 + cindex2$se^2 - 2 * r * cindex1$se * cindex2$se)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "cindex1"=cindex1$c.index, "cindex2"=cindex2$c.index))
}

