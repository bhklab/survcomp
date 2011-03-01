##########
## chunk 1: source biocLite, install and load survcomp
##########
source("http://bioconductor.org/biocLite.R")
biocLite("survcomp")
library("survcomp")

##########
## chunk 2: show survcomp help page
##########
library(help=survcomp)

##########
## chunk 3: load Suggests: packages
##########
library(Biobase)

##########
## chunk 4: load sample dataset
##########
data(sampleData)

##########
## chunk 5: specify 7 gene signature, spaces and censored time point
##########
gsList <- c("ESR1", "ERBB2", "AURKA", "PLAU", "VEGF", "STAT1", "CASP3")
gidList <- c(2099, 2064, 6790, 5328, 7422, 6772, 836)
probesNKI <- c("NM_000125", "NM_004448", "NM_003600", "NM_002658", "NM_003376", "NM_007315", "NM_004346")
probesAffy <- c("205225_at", "216836_s_at", "208079_s_at", "211668_s_at", "211527_x_at", "209969_s_at", "202763_at")
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI","ALL")
myspace <- " "
mybigspace <- "    "
tc <- 10 * 365

##########
## chunk 6: compute concordance indices for all datasets
##########
cindexall.mainz.small <- t(apply(X=exprs(mainzSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainzSample)[ ,"t.dmfs"], z=pData(mainzSample)[ ,"e.dmfs"]))

cindexall.transbig.small <- t(apply(X=exprs(transbigSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbigSample)[ ,"t.dmfs"], z=pData(transbigSample)[ ,"e.dmfs"]))

cindexall.vdx.small <- t(apply(X=exprs(vdxSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdxSample)[ ,"t.dmfs"], z=pData(vdxSample)[ ,"e.dmfs"])) 

cindexall.upp.small <- t(apply(X=exprs(uppSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(uppSample)[ ,"t.rfs"], z=pData(uppSample)[ ,"e.rfs"]))

cindexall.unt.small <- t(apply(X=exprs(untSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(untSample)[ ,"t.dmfs"], z=pData(untSample)[ ,"e.dmfs"]))

cindexall.nki.small <- t(apply(X=exprs(nkiSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nkiSample)[ ,"t.dmfs"], z=pData(nkiSample)[ ,"e.dmfs"]))


##########
## chunk 7: compute D.indices for all datasets
##########
dindexall.mainz.small <- t(apply(X=exprs(mainzSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainzSample)[ ,"t.dmfs"], z=pData(mainzSample)[ ,"e.dmfs"]))

dindexall.transbig.small <- t(apply(X=exprs(transbigSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbigSample)[ ,"t.dmfs"], z=pData(transbigSample)[ ,"e.dmfs"]))

dindexall.upp.small <- t(apply(X=exprs(uppSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(uppSample)[ ,"t.rfs"], z=pData(uppSample)[ ,"e.rfs"]))

dindexall.unt.small <- t(apply(X=exprs(untSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(untSample)[ ,"t.dmfs"], z=pData(untSample)[ ,"e.dmfs"]))

dindexall.vdx.small <- t(apply(X=exprs(vdxSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdxSample)[ ,"t.dmfs"], z=pData(vdxSample)[ ,"e.dmfs"]))

dindexall.nki.small <- t(apply(X=exprs(nkiSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nkiSample)[ ,"t.dmfs"], z=pData(nkiSample)[ ,"e.dmfs"]))


##########
## chunk pre 8: define rescaling function
##########
rescale <- function(x, na.rm=FALSE, q=0) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return(x)
}


##########
## chunk 8: compute hazard.ratios for all datasets
##########
hratio.mainz.small <- t(apply(X=((rescale(exprs(mainzSample) , q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainzSample)[ ,"t.dmfs"], z=pData(mainzSample)[ ,"e.dmfs"]))

hratio.transbig.small <- t(apply(X=((rescale(exprs(transbigSample), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbigSample)[ ,"t.dmfs"], z=pData(transbigSample)[ ,"e.dmfs"]))

hratio.upp.small <- t(apply(X=((rescale(exprs(uppSample), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(uppSample)[ ,"t.rfs"], z=pData(uppSample)[ ,"e.rfs"]))

hratio.unt.small <- t(apply(X=((rescale(exprs(untSample), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(untSample)[ ,"t.dmfs"], z=pData(untSample)[ ,"e.dmfs"]))

hratio.vdx.small <- t(apply(X=((rescale(exprs(vdxSample), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdxSample)[ ,"t.dmfs"], z=pData(vdxSample)[ ,"e.dmfs"]))

hratio.nki.small <- t(apply(X=((rescale(exprs(nkiSample), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nkiSample)[ ,"t.dmfs"], z=pData(nkiSample)[ ,"e.dmfs"]))


##########
## chunk 9: combining the cindices for each gene for all datasets
##########
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(x=cbind(   cindexall.mainz.small[i,"cindex"],
                              cindexall.transbig.small[i,"cindex"],
                              cindexall.upp.small[i,"cindex"],
                              cindexall.unt.small[i,"cindex"],
                              cindexall.vdx.small[i,"cindex"],
                              cindexall.nki.small[i,"cindex"]),
                   x.se=cbind(cindexall.mainz.small[i,"cindex"],
                              cindexall.transbig.small[i,"cindex.se"],
                              cindexall.upp.small[i,"cindex.se"],
                              cindexall.unt.small[i,"cindex.se"],
                              cindexall.vdx.small[i,"cindex.se"],
                              cindexall.nki.small[i,"cindex.se"]),)
              )
}
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
rownames(tt) <- gsList
colnames(tt) <- c("cindex","cindex.se","lower","upper")
ccindex <- tt


##########
## chunk 10: combining the D.indices for each gene for all datasets
##########
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(x=cbind(   dindexall.mainz.small[i,"dindex"],
                              dindexall.transbig.small[i,"dindex"],
                              dindexall.upp.small[i,"dindex"],
                              dindexall.unt.small[i,"dindex"],
                              dindexall.vdx.small[i,"dindex"],
                              dindexall.nki.small[i,"dindex"]),
                   x.se=cbind(dindexall.mainz.small[i,"dindex.se"],
                              dindexall.transbig.small[i,"dindex.se"],
                              dindexall.upp.small[i,"dindex.se"],
                              dindexall.unt.small[i,"dindex.se"],
                              dindexall.vdx.small[i,"dindex.se"],
                              dindexall.nki.small[i,"dindex.se"]),)
              )
}
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
rownames(tt) <- gsList
colnames(tt) <- c("dindex","dindex.se","lower","upper")
cdindex <- tt


##########
## chunk 11: combining all hazard ratios for each gene for all datasets
##########
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(x=cbind(   hratio.mainz.small[i,"hratio"],
                              hratio.transbig.small[i,"hratio"],
                              hratio.upp.small[i,"hratio"],
                              hratio.unt.small[i,"hratio"],
                              hratio.vdx.small[i,"hratio"],
                              hratio.nki.small[i,"hratio"]),
                   x.se=cbind(hratio.mainz.small[i,"hratio.se"],
                              hratio.transbig.small[i,"hratio.se"],
                              hratio.upp.small[i,"hratio.se"],
                              hratio.unt.small[i,"hratio.se"],
                              hratio.vdx.small[i,"hratio.se"],
                              hratio.nki.small[i,"hratio.se"]),)
              )
}
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
rownames(tt) <- gsList
colnames(tt) <- c("hratio","hratio.se","lower","upper")
chratio <- tt


##########
## chunk 12: forestplot for each gene combining all datasets
##########
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
r.mean <- c(NA,ccindex$cindex)
r.lower <- c(NA,ccindex$cindex + qnorm(0.025, lower.tail=TRUE) * ccindex$cindex.se)
r.upper <- c(NA,ccindex$cindex + qnorm(0.025, lower.tail=FALSE) * ccindex$cindex.se)

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.05), xlab=paste("Concordance Index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.4,1))
##title(paste("cindex forestplot, each gene, ccindex"))

##########
## chunk 13: forestplot for each gene combining all datasets
##########
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
tt <- log2(cdindex)

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1.3,0.1), xlab=paste("log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,1.3))
##title(paste("log2 dindex forestplot, each gene, cdindex"))


##########
## chunk 14: forestplot for each gene combining all datasets
##########
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(mybigspace,8)))
bs <- rep(0.5, nrow(labeltext))   
tt <- log2(chratio)

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$hratio), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,3.5,0.5), xlab=paste("log2 Hazard Ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.75,3.5))
title(paste("log2 hratio, each gene, all datasets"))


##########
## chunk 15: forestplot for the cindices from AURKA for all datasets
##########
tt <- rbind(cindexall.mainz.small[3,],
            cindexall.transbig.small[3,],
            cindexall.upp.small[3,],
            cindexall.unt.small[3,],
            cindexall.vdx.small[3,],
            cindexall.nki.small[3,],
            as.numeric(ccindex[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(mybigspace,length(datasetList))))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA Concordance Index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.5,1), is.summary=(c(rep(FALSE,7),TRUE)))
##title(paste("cindex forestplot, AURKA"))


#######################
## chunk 15.1: forestplot for the cindices from VEGF for all datasets
#######################
tt <- rbind(cindexall.mainz.small[5,],
            cindexall.transbig.small[5,],
            cindexall.upp.small[5,],
            cindexall.unt.small[5,],
            cindexall.vdx.small[5,],
            cindexall.nki.small[5,],
            as.numeric(ccindex[5,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetList),c(mybigspace,rep(mybigspace,length(datasetList))))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.75,0.05), xlab=paste("VEGF cindex", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.4,0.75), is.summary=(c(rep(FALSE,7),TRUE)))
##title(paste("cindex forestplot, VEGF"))


#######################
## chunk 15.2: forestplot for the cindices from AURKA and VEGF for all datasets
#######################
tt <- rbind(cindexall.mainz.small[3,],
            cindexall.mainz.small[5,],
            NA,
            cindexall.transbig.small[3,],            
            cindexall.transbig.small[5,],
            NA,
            cindexall.upp.small[3,],            
            cindexall.upp.small[5,],
            NA,
            cindexall.unt.small[3,],            
            cindexall.unt.small[5,],
            NA,
            cindexall.vdx.small[3,],            
            cindexall.vdx.small[5,],
            NA,
            cindexall.nki.small[3,],
            cindexall.nki.small[5,],
            NA,
            as.numeric(ccindex[3,]),            
            as.numeric(ccindex[5,]))

rownames(tt) <- c("MAINZa","MAINZv","a","TRANSBIGa","TRANSBIGv","b","UPPa","UPPv","c","UNTa","UNTv","d","VDXa","VDXv","e","NKIa","NKIv","f","ALLa","ALLv")
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset","MAINZ",NA,NA,"TRANSBIG",NA,NA,"UPP",NA,NA,"UNT",NA,NA,"VDX",NA,NA,"NKI",NA,NA,"ALL",NA),
                   c("Gene",rep(c("aurka","vegf",NA),length(datasetList)-1),c("aurka","vegf")),
                   c(rep(mybigspace,21)))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA and VEGF cindex", myspace, sep=""), col=meta.colors(line=c(rep(c(NA,"darkblue","darkred"),7)),zero="darkred",box=c(rep(c(NA,"royalblue","red"),7))), box.size=bs, clip=c(0.4,1), is.summary=(c(rep(FALSE,19),TRUE,TRUE)))
##title(paste("cindex forestplot, AURKA and VEGF"))    



##########
## chunk 16: forestplot for the D.indices from AURKA for all datasets
##########
tt <- rbind(dindexall.mainz.small[3,],
            dindexall.transbig.small[3,],
            dindexall.upp.small[3,],
            dindexall.unt.small[3,],
            dindexall.vdx.small[3,],
            dindexall.nki.small[3,],
            as.numeric(cdindex[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(mybigspace,8)))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,2,0.5), xlab=paste("AURKA log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.25,2), is.summary=(c(rep(FALSE,7),TRUE)))
##title(paste("log2 dindex forestplot, AURKA"))

##########
## chunk 16: forestplot for the D.indices from VEGF for all datasets
##########
tt <- rbind(dindexall.mainz.small[5,],
            dindexall.transbig.small[5,],
            dindexall.upp.small[5,],
            dindexall.unt.small[5,],
            dindexall.vdx.small[5,],
            dindexall.nki.small[5,],
            as.numeric(cdindex[5,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",datasetList),c(mybigspace,rep(mybigspace,length(datasetList))))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1.5,0.5), xlab=paste("VEGF log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,1.5), is.summary=(c(rep(FALSE,7),TRUE)))
##title(paste("log2 D.index forestplot, VEGF"))



##########
## chunk 17: forestplot for the hazard ratios from AURKA for all datasets
##########
tt <- rbind(hratio.mainz.small[3,],
            hratio.transbig.small[3,],
            hratio.upp.small[3,],
            hratio.unt.small[3,],
            hratio.vdx.small[3,],
            hratio.nki.small[3,],
            as.numeric(chratio[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))
   
forestplot.surv(labeltext=labeltext, mean=c(NA,tt$hratio), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,3.5,0.5), xlab=paste("AURKA log2 Hazard Ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,3.5),is.summary=(c(rep(FALSE,7),TRUE)))
##title(paste("log2 hratio forestplot, AURKA"))


##########
## chunk 18: forestplot for every single gene over all datasets
##########
for(i in 1:length(gsList)) {
##  par(mfrow=c(3,3))
  myspace <- " "
  tt <- rbind(cindexall.mainz.small[i,],
              cindexall.transbig.small[i,],
              cindexall.upp.small[i,],
              cindexall.unt.small[i,],
              cindexall.vdx.small[i,],
              cindexall.nki.small[i,],
              as.numeric(ccindex[i,]))

  rownames(tt) <- datasetList
  tt <- as.data.frame(tt)
  labeltext <- cbind(c("Dataset",datasetList),c(rep(myspace,8)))
  bs <- rep(0.5, nrow(labeltext))

  r.mean <- c(NA,tt$cindex)
  r.lower <- c(NA,tt$cindex + qnorm(0.025, lower.tail=TRUE) * tt$cindex.se)
  r.upper <- c(NA,tt$cindex + qnorm(0.025, lower.tail=FALSE) * tt$cindex.se)
  
  x.ticks.lower <- (floor((min(r.mean,na.rm=TRUE) - 0.1) * 10)/10)
  x.ticks.upper <- (floor((max(r.mean,na.rm=TRUE) + 0.2) * 10)/10)
   
  forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(x.ticks.lower,x.ticks.upper,0.05), xlab=paste(gsList[i], myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.3,0.8),is.summary=(c(rep(FALSE,7),TRUE)))
##  title(paste("cindex forestplot, ", gsList[i]))
}

##########
## chunk 19: kaplan meier survival curve for all datasets
##########
surv.data <- censor.time(surv.time=c(pData(mainzSample)[ ,"t.dmfs"], pData(transbigSample)[ ,"t.dmfs"], pData(untSample)[ ,"t.dmfs"], pData(vdxSample)[ ,"t.dmfs"], pData(uppSample)[ ,"t.rfs"], pData(nkiSample)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainzSample)[ ,"e.dmfs"], pData(transbigSample)[ ,"e.dmfs"], pData(untSample)[ ,"e.dmfs"], pData(vdxSample)[ ,"e.dmfs"], pData(uppSample)[ ,"e.rfs"], pData(nkiSample)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainzSample))), rep("transbig", nrow(pData(transbigSample))), rep("unt", nrow(pData(untSample))), rep("vdx", nrow(pData(vdxSample))), rep("upp", nrow(pData(uppSample))), rep("nki", nrow(pData(nkiSample)))), levels=c("mainz", "transbig", "unt", "vdx", "upp", "nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange", "red", "darkblue", "darkgreen", "black", "brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
##title(paste("km survival curve, all datasets"))


##########
## chunk 20: kaplan meier survival curve for all six datasets
##########
##par(mfrow=c(2,3))
surv.data <- censor.time(surv.time=c(pData(mainzSample)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainzSample)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainzSample)))), levels=c("mainz"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
##title(paste("km survival curve, MAINZ"))

surv.data <- censor.time(surv.time=c(pData(transbigSample)[ ,"t.dmfs"]) / 365, surv.event=c(pData(transbigSample)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("transbig", nrow(pData(transbigSample)))), levels=c("transbig"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("red"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
##title(paste("km survival curve, TRANSBIG"))

surv.data <- censor.time(surv.time=c(pData(uppSample)[ ,"t.rfs"]) / 365, surv.event=c(pData(uppSample)[ ,"e.rfs"]), time.cens=tc / 365)
gg <- factor(c(rep("upp", nrow(pData(uppSample)))), levels=c("upp"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of RFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE) 
##title(paste("km survival curve, UPP"))

surv.data <- censor.time(surv.time=c(pData(untSample)[ ,"t.dmfs"]) / 365, surv.event=c(pData(untSample)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("unt", nrow(pData(untSample)))), levels=c("unt"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("black"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
##title(paste("km survival curve, UNT"))

surv.data <- censor.time(surv.time=c(pData(vdxSample)[ ,"t.dmfs"]) / 365, surv.event=c(pData(vdxSample)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("vdx", nrow(pData(vdxSample)))), levels=c("vdx"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkblue"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
##title(paste("km survival curve, VDX"))

surv.data <- censor.time(surv.time=c(pData(nkiSample)[ ,"t.dmfs"]) / 365, surv.event=c(pData(nkiSample)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("nki", nrow(pData(nkiSample)))), levels=c("nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkgreen"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)      
##title(paste("km survival curve, NKI"))


##########
## chunk 21: kaplan meier survival plot for AURKA combining 6 datasets
##########
aurkaGs <- "AURKA"
aurkaGid <- 6790
aurkaPaf <- "208079_s_at"
aurkaPagi <- "NM_003600"

surv.time.all <- c(pData(mainzSample)[ ,"t.dmfs"], pData(transbigSample)[ ,"t.dmfs"], pData(untSample)[ ,"t.dmfs"], pData(uppSample)[ ,"t.rfs"], pData(vdxSample)[ ,"t.dmfs"], pData(nkiSample)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainzSample)[ ,"e.dmfs"], pData(transbigSample)[ ,"e.dmfs"], pData(untSample)[ ,"e.dmfs"], pData(uppSample)[ ,"e.rfs"], pData(vdxSample)[ ,"e.dmfs"], pData(nkiSample)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainzSample)[aurkaPaf,], exprs(transbigSample)[aurkaPaf,], exprs(untSample)[aurkaPaf,], exprs(uppSample)[aurkaPaf,], exprs(vdxSample)[aurkaPaf,], exprs(nkiSample)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainzSample)[aurkaPaf,]), length(exprs(transbigSample)[aurkaPaf,]), length(exprs(untSample)[aurkaPaf,]), length(exprs(uppSample)[aurkaPaf,]), length(exprs(vdxSample)[aurkaPaf,]), length(exprs(nkiSample)[aurkaPagi,]))

pos <- 0
mygroup <- NULL
for(i in aurka.exprs.length){
  qq <- aurka.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[aurka.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[aurka.exprs[(pos+1):(pos+i)] >= myq[1] & aurka.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[aurka.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainzSample))), rep("transbig", nrow(pData(transbigSample))), rep("unt", nrow(pData(untSample))), rep("upp", nrow(pData(uppSample))), rep("vdx", nrow(pData(vdxSample))), rep("nki", nrow(pData(nkiSample)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="", sub.title=NULL, leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE)


##########
## chunk 22: SBRIER
##########
dd.tr <- data.frame("time"=pData(mainzSample)[,"t.dmfs"], "event"=pData(mainzSample)[,"e.dmfs"], "score"=pData(mainzSample)[,"age"])
##sbrier.score2proba(data.tr=dd, data.ts=dd, method="cox")

#brier
mysbrier <- sbrier.score2proba(data.tr=dd.tr, data.ts=dd.tr, method="cox")
#time-dependent ROC curves
tt <- unique(sort(dd.tr$time[dd.tr$event == 1]))
if(max(tt) < tc / 365) { tt <- c(tt, tc / 365) }
mytdroc <- NULL
for(i in 1:length(tt)) {
	rr <- tdrocc(x=dd.tr$score, surv.time=dd.tr$time, surv.event=dd.tr$event, time=tt[i], na.rm=TRUE, verbose=FALSE)
	mytdroc <- c(mytdroc, list(rr))
}
names(mytdroc) <- paste("years", tt, sep=".")


myperf <- mysbrier##sbrier.score2proba.tr.all

mycol <- c("#000000FF", rainbow(length(myperf))[-1])
mylty <- 1:(length(myperf))
mylwd <- c(3, rep(2, length(myperf)-1))

plot(myperf[[1]]$time, myperf[[1]]$bsc, typ="l", xlab="Time (years)", ylab="Brier score", col=mycol[i], lty=mylty[i], ylim=c(0,0.4), lwd=mylwd[i], main="Brier Score")

for(i in 1:length(myperf)) {
	if(i == 1) {
		plot(myperf[[i]]$time, myperf[[i]]$bsc, typ="l", xlab="Time (years)", ylab="Brier score", col=mycol[i], lty=mylty[i], ylim=c(0,0.4), lwd=mylwd[i], main="Brier Score")
	} else {  lines(myperf[[i]]$time, myperf[[i]]$bsc, col=mycol[i], lty=mylty[i], lwd=mylwd[i]) }
	mysbrier.int <- unlist(lapply(myperf, function(x) { return(x$bsc.integrated)}))
}
smartlegend(x="left", y="top", legend=sprintf("%s, IBSC = %.3g", names(myperf), mysbrier.int), col=mycol, lty=mylty, lwd=mylwd)




##########
## chunk 23: cindex.comp.meta example for AURKA and VEGF in MAINZ dataset
##########
#cindex.meta.mainz <- t(apply(X=exprs(mainzSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(tt); }, y=pData(mainzSample)[ ,"t.dmfs"], z=pData(mainzSample)[ ,"e.dmfs"]))
#
#ccmMainz <- tt <- rr <- NULL
#ccmNames <- colnames(cindex.meta.mainz)
#
#for(i in 1:length(cindex.meta.mainz)){
#  tt <- NULL
#  for(j in 1:length(cindex.meta.mainz)){
#    rr <- cindex.comp.meta(list.cindex1=list("cindex.meta.mainz"=cindex.meta.mainz[[i]]), list.cindex2=list("cindex.meta.mainz"=cindex.meta.mainz[[j]]))
#    tt <- cbind(tt, list(round(rr$p.value,4)))
#  }
#  ccmMainz <- rbind(ccmMainz, tt)
#}
#ccmMainz <- as.data.frame(ccmMainz)
#colnames(ccmMainz) <- gsList
#rownames(ccmMainz) <- gsList
#
#latex(format.df(ccmMainz), file="")



##########
## chunk 23: cindex.comp.meta example for AURKA and VEGF in MAINZ dataset
##########
cindexMetaMainz <- t(apply(X=exprs(mainzSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(tt); }, y=pData(mainzSample)[ ,"t.dmfs"], z=pData(mainzSample)[ ,"e.dmfs"]))

cindexMetaTransbig <- t(apply(X=exprs(transbigSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(tt); }, y=pData(transbigSample)[ ,"t.dmfs"], z=pData(transbigSample)[ ,"e.dmfs"]))

cindexMetaUpp <- t(apply(X=exprs(uppSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(tt); }, y=pData(uppSample)[ ,"t.rfs"], z=pData(uppSample)[ ,"e.rfs"]))

cindexMetaUnt <- t(apply(X=exprs(untSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(tt); }, y=pData(untSample)[ ,"t.dmfs"], z=pData(untSample)[ ,"e.dmfs"]))

cindexMetaVdx <- t(apply(X=exprs(vdxSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(tt); }, y=pData(vdxSample)[ ,"t.dmfs"], z=pData(vdxSample)[ ,"e.dmfs"]))


ccNki <- complete.cases(exprs(nkiSample)[1,],exprs(nkiSample)[2,],exprs(nkiSample)[3,],exprs(nkiSample)[4,],exprs(nkiSample)[5,],exprs(nkiSample)[6,],exprs(nkiSample)[7,],pData(nkiSample)[,"e.dmfs"],pData(nkiSample)[,"e.dmfs"])

cindexMetaNki <- t(apply(X=exprs(nkiSample)[,ccNki], MARGIN=1, function(x, y, z) {tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(tt); }, y=pData(nkiSample)[ccNki ,"t.dmfs"], z=pData(nkiSample)[ccNki ,"e.dmfs"]))


ccmData <- tt <- rr <- NULL
ccmNames <- colnames(cindexMetaMainz)       

for(i in 1:7){
  tt <- NULL
  firstOne <- list("mainz"=cindexMetaMainz[[i]],
                   "transbig"=cindexMetaTransbig[[i]],
                   "upp"=cindexMetaUpp[[i]],
                   "unt"=cindexMetaUnt[[i]],
                   "vdx"=cindexMetaVdx[[i]],
                   "nki"=cindexMetaNki[[i]])
                   
  for(j in 1:7){
    listTwo <- list("mainz"=cindexMetaMainz[[j]],
                    "transbig"=cindexMetaTransbig[[j]],
                    "upp"=cindexMetaUpp[[j]],
                    "unt"=cindexMetaUnt[[j]],
                    "vdx"=cindexMetaVdx[[j]],
                    "nki"=cindexMetaNki[[j]])
                    
    paste("iteration [i,j]: ", i, ", ", j)
    rr <- cindex.comp.meta(list.cindex1=listOne, list.cindex2=listTwo)
    tt <- cbind(tt, list(round(rr$p.value,4)))
  }
  ccmData <- rbind(ccmData, tt)
}
ccmData <- as.data.frame(ccmData)
colnames(ccmData) <- gsList
rownames(ccmData) <- gsList

latex(format.df(ccmData),file="")





##########
## 
##########



