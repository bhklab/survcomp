###################################################
### chunk number 1: setup
###################################################
## #line 72 "survcomp.Rnw"
## library(pgfSweave)
## setCacheDir("cache")
## options(keep.source=TRUE)


###################################################
### chunk number 2: install-pkg eval=FALSE
###################################################
## #line 97 "survcomp.Rnw"
## source("http://bioconductor.org/biocLite.R")
## biocLite("survcomp")


###################################################
### chunk number 3: loadlib
###################################################
#line 103 "survcomp.Rnw"
library(survcomp)


###################################################
### chunk number 4: survcomphelp eval=FALSE
###################################################
## #line 113 "survcomp.Rnw"
## library(help=survcomp)


###################################################
### chunk number 5: loadDepends
###################################################
#line 137 "survcomp.Rnw"
library(Biobase)
library(xtable)


###################################################
### chunk number 6: loadbreastCancerData
###################################################
#line 144 "survcomp.Rnw"
data(breastCancerData)
mainz7g


###################################################
### chunk number 7: createVars
###################################################
#line 151 "survcomp.Rnw"
gsList <- tolower(fData(mainz7g)[,"Gene.symbol"])
gidList <- fData(mainz7g)[,"Gene.ID"]
probesNKI <- as.character(fData(nki7g)[,"probe"])
probesAffy <- fData(mainz7g)[,"probe"]
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI","","Overall")
myspace <- " "
mybigspace <- "    "
tc <- 10 * 365


###################################################
### chunk number 8: showVars
###################################################
#line 164 "survcomp.Rnw"
overview <- as.data.frame(cbind("Gene Symbol"=gsList,"Gene ID"=gidList,"Probes Agilent"=probesNKI,"Probes Affy"=probesAffy))
print(xtable(overview), floating=FALSE)


###################################################
### chunk number 9: computeCindexSample
###################################################
#line 178 "survcomp.Rnw"
cindexall.mainz.small <- t(apply(X=exprs(mainz7g), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))


###################################################
### chunk number 10: computeCindex
###################################################
#line 182 "survcomp.Rnw"
cindexall.transbig.small <- t(apply(X=exprs(transbig7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))

cindexall.vdx.small <- t(apply(X=exprs(vdx7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"])) 

cindexall.upp.small <- t(apply(X=exprs(upp7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))

cindexall.unt.small <- t(apply(X=exprs(unt7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))

cindexall.nki.small <- t(apply(X=exprs(nki7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki7g)[ ,"t.dmfs"], z=pData(nki7g)[ ,"e.dmfs"]))


###################################################
### chunk number 11: dindexSample
###################################################
#line 211 "survcomp.Rnw"
dindexall.mainz.small <- t(apply(X=exprs(mainz7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))


###################################################
### chunk number 12: dindex
###################################################
#line 218 "survcomp.Rnw"
dindexall.transbig.small <- t(apply(X=exprs(transbig7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))

dindexall.upp.small <- t(apply(X=exprs(upp7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))

dindexall.unt.small <- t(apply(X=exprs(unt7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))

dindexall.vdx.small <- t(apply(X=exprs(vdx7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"]))

dindexall.nki.small <- t(apply(X=exprs(nki7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki7g)[ ,"t.dmfs"], z=pData(nki7g)[ ,"e.dmfs"]))


###################################################
### chunk number 13: showGeneExpressionRange
###################################################
#line 247 "survcomp.Rnw"
tt <- list(mainz7g, transbig7g, upp7g, unt7g, vdx7g, nki7g)
ttNames <- c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI")
dataRange <- NULL
for(i in tt) {
  rr <- range(exprs(i), na.rm=TRUE)
  dataRange$Min <- rbind(dataRange$Min, rr[1])
  dataRange$Max <- rbind(dataRange$Max, rr[2])
}
dataRange <- as.data.frame(dataRange)
rownames(dataRange) <- ttNames


###################################################
### chunk number 14: dataRangeTable
###################################################
#line 262 "survcomp.Rnw"
print(xtable(dataRange), floating=FALSE)


###################################################
### chunk number 15: rescaleExpressionData
###################################################
#line 272 "survcomp.Rnw"
rescale <- function(x, na.rm=FALSE, q=0.05) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return((x - 0.5) * 2)
}


###################################################
### chunk number 16: hratioSample
###################################################
#line 283 "survcomp.Rnw"
hratio.mainz.small <- t(apply(X=rescale(exprs(mainz7g) , q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))


###################################################
### chunk number 17: hratio
###################################################
#line 290 "survcomp.Rnw"
hratio.transbig.small <- t(apply(X=rescale(exprs(transbig7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))

hratio.upp.small <- t(apply(X=rescale(exprs(upp7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))

hratio.unt.small <- t(apply(X=rescale(exprs(unt7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))

hratio.vdx.small <- t(apply(X=rescale(exprs(vdx7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"]))

hratio.nki.small <- t(apply(X=rescale(exprs(nki7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki7g)[ ,"t.dmfs"], z=pData(nki7g)[ ,"e.dmfs"]))


###################################################
### chunk number 18: combineCindices
###################################################
#line 322 "survcomp.Rnw"
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(
    x=cbind(   cindexall.mainz.small[i,"cindex"],
        cindexall.transbig.small[i,"cindex"],
        cindexall.upp.small[i,"cindex"],
        cindexall.unt.small[i,"cindex"],
        cindexall.vdx.small[i,"cindex"],
        cindexall.nki.small[i,"cindex"]),
    x.se=cbind(cindexall.mainz.small[i,"cindex.se"],
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


###################################################
### chunk number 19: combineCindicesOutput
###################################################
#line 352 "survcomp.Rnw"
print(xtable(ccindex), floating=FALSE)


###################################################
### chunk number 20: combineDindices
###################################################
#line 364 "survcomp.Rnw"
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(
    x=cbind( dindexall.mainz.small[i,"dindex"],
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
print(xtable(log2(cdindex)), floating=FALSE)


###################################################
### chunk number 21: combineHratio
###################################################
#line 399 "survcomp.Rnw"
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(
    x=cbind( hratio.mainz.small[i,"hratio"],
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
print(xtable(log2(chratio)), floating=FALSE)


###################################################
### chunk number 22: allDatasetsForestPlotCindex
###################################################
#line 435 "survcomp.Rnw"
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
r.mean <- c(NA,ccindex$cindex)
r.lower <- c(NA,ccindex$lower)
r.upper <- c(NA,ccindex$upper)

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.05), xlab=paste("Concordance Index", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(0.4,1))


###################################################
### chunk number 23: allDatasetsForestPlotDindex
###################################################
#line 449 "survcomp.Rnw"
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
r.mean <- c(NA,log2(cdindex$dindex))
r.lower <- c(NA,log2(cdindex$lower))
r.upper <- c(NA,log2(cdindex$upper))

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1,0.1), xlab=paste("log2 D Index", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(-0.5,1.25))


###################################################
### chunk number 24: allDatasetsForestPlotHratio
###################################################
#line 463 "survcomp.Rnw"
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(mybigspace,8)))
bs <- rep(0.5, nrow(labeltext))   
r.mean <- c(NA,log2(chratio$hratio))
r.lower <- c(NA,log2(chratio$lower))
r.upper <- c(NA,log2(chratio$upper))

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,3.5,0.5), xlab=paste("log2 Hazard Ratio", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(-0.75,3.5))


###################################################
### chunk number 25: AURKAForestPlotCindices
###################################################
#line 477 "survcomp.Rnw"
tt <- rbind(cindexall.mainz.small[3,],
            cindexall.transbig.small[3,],
            cindexall.upp.small[3,],
            cindexall.unt.small[3,],
            cindexall.vdx.small[3,],
            cindexall.nki.small[3,],
            NA,
            as.numeric(ccindex[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(mybigspace,length(datasetList)+1)))
bs <- rep(0.5, nrow(labeltext))
r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                               

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA Concordance Index", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(0.5,1), is.summary=(c(rep(FALSE,8),TRUE)))


###################################################
### chunk number 26: VEGFForestPlotCindices
###################################################
#line 500 "survcomp.Rnw"
tt <- rbind(cindexall.mainz.small[5,],
            cindexall.transbig.small[5,],
            cindexall.upp.small[5,],
            cindexall.unt.small[5,],
            cindexall.vdx.small[5,],
            cindexall.nki.small[5,],
            NA,
            as.numeric(ccindex[5,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(mybigspace,length(datasetList)+1)))
bs <- rep(0.5, nrow(labeltext))
r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.75,0.05), xlab=paste("VEGF Concordance Index", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(0.3,0.75), is.summary=(c(rep(FALSE,8),TRUE)))


###################################################
### chunk number 27: AURKAandVEGFForestPlotCindices
###################################################
#line 525 "survcomp.Rnw"
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

rownames(tt) <- c("MAINZa", "MAINZv", "a", "TRANSBIGa", "TRANSBIGv", "b", "UPPa", "UPPv", "c", "UNTa", "UNTv", "d", "VDXa", "VDXv", "e", "NKIa", "NKIv", "f", "ALLa", "ALLv")
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset", "MAINZ", NA, NA, "TRANSBIG", NA, NA, "UPP", NA, NA, "UNT", NA, NA, "VDX", NA, NA, "NKI", NA, NA, "Overall", NA),
                   c("Gene", rep(c("aurka","vegf",NA), length(datasetList)-2), c("aurka","vegf")), c(rep(mybigspace,21)))
bs <- rep(0.5, nrow(labeltext))   
r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA and VEGF Concordance Index", myspace, sep=""),
    col=meta.colors(line=c(rep(c(NA, "darkblue", "seagreen"),7)), zero="firebrick", box=c(rep(c(NA," royalblue", "forestgreen"),7))), box.size=bs,
    clip=c(0.3,1), is.summary=(c(rep(FALSE,19), TRUE, TRUE)))


###################################################
### chunk number 28: AURKAForestPlotDindices
###################################################
#line 564 "survcomp.Rnw"
tt <- rbind(dindexall.mainz.small[3,],
            dindexall.transbig.small[3,],
            dindexall.upp.small[3,],
            dindexall.unt.small[3,],
            dindexall.vdx.small[3,],
            dindexall.nki.small[3,],
            NA,
            as.numeric(cdindex[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(mybigspace,length(datasetList)+1)))
bs <- rep(0.5, nrow(labeltext))                              
r.mean <- c(NA,log2(tt$dindex))
r.lower <- c(NA,log2(tt$lower))
r.upper <- c(NA,log2(tt$upper))

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,2,0.5), xlab=paste("AURKA log2 D Index", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(-0.25,2), is.summary=(c(rep(FALSE,8), TRUE)))


###################################################
### chunk number 29: VEGFForestPlotDindices
###################################################
#line 587 "survcomp.Rnw"
tt <- rbind(dindexall.mainz.small[5,],
            dindexall.transbig.small[5,],
            dindexall.upp.small[5,],
            dindexall.unt.small[5,],
            dindexall.vdx.small[5,],
            dindexall.nki.small[5,],
            NA,
            as.numeric(cdindex[5,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(mybigspace,length(datasetList)+1)))
bs <- rep(0.5, nrow(labeltext))  
r.mean <- c(NA,log2(tt$dindex))
r.lower <- c(NA,log2(tt$lower))
r.upper <- c(NA,log2(tt$upper))                            

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-1.25,1.5,0.25), xlab=paste("VEGF log2 D Index", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(-1.5,1.75), is.summary=(c(rep(FALSE,8), TRUE)))


###################################################
### chunk number 30: AURKAForestPlotHratio
###################################################
#line 612 "survcomp.Rnw"
tt <- rbind(hratio.mainz.small[3,],
            hratio.transbig.small[3,],
            hratio.upp.small[3,],
            hratio.unt.small[3,],
            hratio.vdx.small[3,],
            hratio.nki.small[3,],
            NA,
            as.numeric(chratio[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(myspace,length(datasetList)+1)))
bs <- rep(0.5, nrow(labeltext))
r.mean <- c(NA,log2(tt$hratio))
r.lower <- c(NA,log2(tt$lower))
r.upper <- c(NA,log2(tt$upper))  
   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,3.5,0.5), xlab=paste("AURKA log2 Hazard Ratio", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(-0.5,3.5),is.summary=(c(rep(FALSE,8), TRUE)))


###################################################
### chunk number 31: singleGenesForestplotCindex eval=FALSE
###################################################
## #line 637 "survcomp.Rnw"
## for(i in 1:length(gsList)) {
##   tt <- rbind(cindexall.mainz.small[i,],
##               cindexall.transbig.small[i,],
##               cindexall.upp.small[i,],
##               cindexall.unt.small[i,],
##               cindexall.vdx.small[i,],
##               cindexall.nki.small[i,],
##               NA,
##               as.numeric(ccindex[i,]))
## 
##   rownames(tt) <- datasetList
##   tt <- as.data.frame(tt)
##   labeltext <- cbind(c("Dataset",datasetList), c(rep(myspace,length(datasetList)+1)))
##   bs <- rep(0.5, nrow(labeltext))
##   r.mean <- c(NA,tt$cindex)
##   r.lower <- c(NA,tt$lower)
##   r.upper <- c(NA,tt$upper)
##   
##   x.ticks.lower <- (floor((min(r.mean,na.rm=TRUE) - 0.1) * 10)/10)
##   x.ticks.upper <- (floor((max(r.mean,na.rm=TRUE) + 0.2) * 10)/10)
##    
##   forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
##       align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(x.ticks.lower,x.ticks.upper,0.05), xlab=paste(gsList[i], " Concordance Index", myspace, sep=""),
##       col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(0.3,0.8),is.summary=(c(rep(FALSE,8), TRUE)))
## }


###################################################
### chunk number 32: survivalCurveAllDatasets
###################################################
#line 670 "survcomp.Rnw"
surv.data <- censor.time(surv.time=c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(vdx7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(nki7g)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(vdx7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(nki7g)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("vdx", nrow(pData(vdx7g))), rep("upp", nrow(pData(upp7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "vdx", "upp", "nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="", sub.title="", leg.pos="bottomright", leg.inset=0.05, o.text=FALSE, v.line=FALSE, h.line=FALSE, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange", "red", "darkblue", "darkgreen", "black", "brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE, ylim=c(0.1,1))


###################################################
### chunk number 33: survivalCurveAURKA
###################################################
#line 679 "survcomp.Rnw"
aurkaGs <- "AURKA"
aurkaGid <- 6790
aurkaPaf <- "208079_s_at"
aurkaPagi <- "NM_003600"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz7g)[aurkaPaf,]), length(exprs(transbig7g)[aurkaPaf,]), length(exprs(unt7g)[aurkaPaf,]), length(exprs(upp7g)[aurkaPaf,]), length(exprs(vdx7g)[aurkaPaf,]), length(exprs(nki7g)[aurkaPagi,]))

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
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))


###################################################
### chunk number 34: cindexcompmetaAnalysis
###################################################
#line 714 "survcomp.Rnw"
cindexMetaMainz <- t(apply(X=exprs(mainz7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); }, y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))

cindexMetaTransbig <- t(apply(X=exprs(transbig7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); }, y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))

cindexMetaUpp <- t(apply(X=exprs(upp7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); }, y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))

cindexMetaUnt <- t(apply(X=exprs(unt7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); }, y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))

cindexMetaVdx <- t(apply(X=exprs(vdx7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); }, y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"]))

ccNki <- complete.cases(exprs(nki7g)[1,], exprs(nki7g)[2,], exprs(nki7g)[3,], exprs(nki7g)[4,], exprs(nki7g)[5,], exprs(nki7g)[6,], exprs(nki7g)[7,], pData(nki7g)[,"e.dmfs"], pData(nki7g)[,"e.dmfs"])

cindexMetaNki <- t(apply(X=exprs(nki7g)[,ccNki], MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); }, y=pData(nki7g)[ccNki ,"t.dmfs"], z=pData(nki7g)[ccNki ,"e.dmfs"]))

ccmData <- tt <- rr <- NULL
for(i in 1:7){
  tt <- NULL
  listOne <- list("mainz"=cindexMetaMainz[[i]],
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
                    
    rr <- cindex.comp.meta(list.cindex1=listOne, list.cindex2=listTwo)
    tt <- cbind(tt, rr$p.value) ##list(round(rr$p.value,5)))
  }
  ccmData <- rbind(ccmData, tt)
}
ccmData <- as.data.frame(ccmData)
colnames(ccmData) <- gsList
rownames(ccmData) <- gsList


###################################################
### chunk number 35: cindexcompmetaAnalysisTable
###################################################
#line 771 "survcomp.Rnw"
print(xtable(ccmData, digits=5), floating=FALSE)


###################################################
### chunk number 36: sessionInfo
###################################################
#line 791 "survcomp.Rnw"
toLatex(sessionInfo())


