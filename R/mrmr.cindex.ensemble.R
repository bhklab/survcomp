`mrmr.cindex.ensemble` <-
function(x, surv.time, surv.event, cl, weights, comppairs=10, strat, alpha=0.05, outx=TRUE, method=c("conservative", "noether", "nam"), alternative=c("two.sided", "less", "greater"), maxparents, maxnsol, nboot=200, na.rm=FALSE) {

	nvar<-ncol(x)
	nsample<-nrow(x)
	
	vec_ensemble<-.Call("mrmr_cIndex_ensemble_remove",data.matrix(x),as.integer(is.na(x)),maxparents,ncol(x),nrow(x),1,1,nboot,maxnsol,-1000,as.integer(as.logical(TRUE)),as.integer(sort(unique(strat))),as.integer(cl ),as.double(surv.time),as.integer(surv.event),as.double(weights),as.integer(strat),as.integer(sum(weights)),as.integer(as.logical(outx)),as.integer(length(strat)),as.integer(length(sort(unique(strat)))))
	
	vec_ensemble[2:vec_ensemble[1]+1]<-vec_ensemble[2:vec_ensemble[1]+1]-1
	
	models.equiv <- .extract.all.parents(x,vec_ensemble,maxparents,1)
	models.equiv[1,]<-rep("T",ncol(models.equiv))
	
	return(models.equiv)
}


