`mrmr.surv.ensemble` <-
function(mydata, surv.time, surv.event, cl, weights, comppairs=10, strat, alpha=0.05, outx=TRUE, na.rm=FALSE, maxparents, maxnsol,nboot=200) {

	nvar<-ncol(mydata)
	nsample<-nrow(mydata)
	
	
	vec_ensemble<-.Call("mrmr_cIndex_ensemble_remove",data.matrix(mydata),as.integer(is.na(mydata)),maxparents,ncol(mydata),nrow(mydata),1,1,nboot,maxnsol,-1000,as.integer(as.logical(TRUE)),as.integer(sort(unique(strat))),as.integer(cl ),as.double(surv.time),as.integer(surv.event),as.double(weights),as.integer(strat),as.integer(sum(weight)),as.integer(as.logical(outx)),as.integer(length(strat)),as.integer(length(sort(unique(strat)))))
	
	vec_ensemble[2:vec_ensemble[1]+1]<-vec_ensemble[2:vec_ensemble[1]+1]-1
	models.equiv <- .extract.all.parents(mydata,vec_ensemble,maxparents,1)
																																																			
	models.equiv[1,]<-rep("T",ncol(models.equiv))
																																																										
	return(models.equiv)																																																										

}


