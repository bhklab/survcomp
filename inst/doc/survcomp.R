################
## chunk 1: setup
################
library(pgfSweave)
setCacheDir("cache")
options(keep.source=TRUE, device=pdf)


################
## chunk 2: install package
################
source("http://bioconductor.org/biocLite.R")
biocLite("survcomp")

################
## chunk 3: load package
################
library(survcomp)


################
## chunk 4: call package help
################
library(help=survcomp)


################
## chunk 5: load dependencies
################
library(Biobase)
library(xtable)


################
## chunk 6: load dataset
################
data(breastCancerData)
mainz7g


################
## chunk 7: create info variables
################
gsList <- c("ESR1", "ERBB2", "AURKA", "PLAU", "VEGF", "STAT1", "CASP3")
gidList <- c(2099, 2064, 6790, 5328, 7422, 6772, 836)
probesNKI <- c("NM_000125", "NM_004448", "NM_003600", "NM_002658", "NM_003376", "NM_007315", "NM_004346")
probesAffy <- c("205225_at", "216836_s_at", "208079_s_at", "211668_s_at", "211527_x_at", "209969_s_at", "202763_at")
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI","","Overall")
myspace <- " "
mybigspace <- "    "
tc <- 10 * 365


################
## chunk 8: compute cindex sample
################
cindexall.mainz.small <- t(apply(X=exprs(mainz7g), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))


################
## chunk 9: compute cindices
################
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


################
## chunk 10: compute D index sample
################
dindexall.mainz.small <- t(apply(X=exprs(mainz7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))


################
## chunk 11: compute D indices
################
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


################
## chunk 12: create table that shows the range of the gene expression values in the six datasets
################
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


################
## chunk 13: display table
################
xtable(dataRange)


################
## chunk 14: create rescaling function
################
rescale <- function(x, na.rm=FALSE, q=0) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return(x)
}


################
## chunk 15: compute hazard ratio sample
################
hratio.mainz.small <- t(apply(X=((rescale(exprs(mainz7g) , q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))


################
## chunk 16: compute hazard ratios
################
hratio.transbig.small <- t(apply(X=((rescale(exprs(transbig7g), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))

hratio.upp.small <- t(apply(X=((rescale(exprs(upp7g), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))

hratio.unt.small <- t(apply(X=((rescale(exprs(unt7g), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))

hratio.vdx.small <- t(apply(X=((rescale(exprs(vdx7g), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"]))

hratio.nki.small <- t(apply(X=((rescale(exprs(nki7g), q=0.05, na.rm=TRUE) - 0.5) * 2), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki7g)[ ,"t.dmfs"], z=pData(nki7g)[ ,"e.dmfs"]))


################
## chunk 17: combine the concordance indices with combine.est()
################
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


################
## chunk 18: display table with the combined concordance indices
################
xtable(ccindex)


################
## chunk 19: combine the D indices with combine.est() and display table with the combined log2 D indices
################
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
xtable(log2(cdindex))


################
## chunk 20: combine the hazard ratios with combine.est() and display table with the combined log2 hazard ratios
################
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
xtable(log2(chratio))


################
## chunk 21: create forestplot for all seven genes showing the combined concordance indices
################
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
r.mean <- c(NA,ccindex$cindex)
r.lower <- c(NA,ccindex$lower)
r.upper <- c(NA,ccindex$upper)

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.05), xlab=paste("Concordance Index", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(0.4,1))


################
## chunk 22: create forestplot for all seven genes showing the combined log2 D indices
################
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
r.mean <- c(NA,log2(cdindex$dindex))
r.lower <- c(NA,log2(cdindex$lower))
r.upper <- c(NA,log2(cdindex$upper))

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1,0.1), xlab=paste("log2 D Index", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(-0.5,1.25))


################
## chunk 23: create forestplot for all seven genes showing the combined log2 hazard ratios
################
labeltext <- cbind(c("Gene Symbol",gsList),c(rep(mybigspace,8)))
bs <- rep(0.5, nrow(labeltext))   
r.mean <- c(NA,log2(chratio$hratio))
r.lower <- c(NA,log2(chratio$lower))
r.upper <- c(NA,log2(chratio$upper))

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,3.5,0.5), xlab=paste("log2 Hazard Ratio", myspace, sep=""),
    col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(-0.75,3.5))


################
## chunk 24: create forestplot for AURKA showing the concordance indices of all six datasets
################
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


################
## chunk 25: create forestplot for VEGF showing the concordance indices of all six datasets
################
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


################
## chunk 26: create advanced forestplot for AURKA and VEGF showing the concordance indices of all six datasets and the combined one
################
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


################
## chunk 27: create forestplot for AURKA showing the D indices of all six datasets and the combined one
################
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


################
## chunk 28: create forestplot for VEGF showing the D indices of all six datasets and the combined one
################
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


################
## chunk 29: create forestplot for AURKA showing the hazard ratios of all six datasets and the combined one
################
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


################
## chunk 30: for loop that creates forestplots for each gene showing the concordance indices of all six datasets and the combined one
################
for(i in 1:length(gsList)) {
  tt <- rbind(cindexall.mainz.small[i,],
              cindexall.transbig.small[i,],
              cindexall.upp.small[i,],
              cindexall.unt.small[i,],
              cindexall.vdx.small[i,],
              cindexall.nki.small[i,],
              NA,
              as.numeric(ccindex[i,]))

  rownames(tt) <- datasetList
  tt <- as.data.frame(tt)
  labeltext <- cbind(c("Dataset",datasetList), c(rep(myspace,length(datasetList)+1)))
  bs <- rep(0.5, nrow(labeltext))
  r.mean <- c(NA,tt$cindex)
  r.lower <- c(NA,tt$lower)
  r.upper <- c(NA,tt$upper)
  
  x.ticks.lower <- (floor((min(r.mean,na.rm=TRUE) - 0.1) * 10)/10)
  x.ticks.upper <- (floor((max(r.mean,na.rm=TRUE) + 0.2) * 10)/10)
   
  forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
      align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(x.ticks.lower,x.ticks.upper,0.05), xlab=paste(gsList[i], " Concordance Index", myspace, sep=""),
      col=meta.colors(box="royalblue", line="darkblue", zero="darkred"), box.size=bs, clip=c(0.3,0.8),is.summary=(c(rep(FALSE,8), TRUE)))
}


################
## chunk 31: create Kaplan Meier survival curve showing the overall survival in percent for the patients from all six datasets
################
surv.data <- censor.time(surv.time=c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(vdx7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(nki7g)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(vdx7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(nki7g)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("vdx", nrow(pData(vdx7g))), rep("upp", nrow(pData(upp7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "vdx", "upp", "nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="", sub.title="", leg.pos="bottomright", leg.inset=0.05, o.text=FALSE, v.line=FALSE, h.line=FALSE, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange", "red", "darkblue", "darkgreen", "black", "brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE, ylim=c(0.3,1))


################
## chunk 32: create Kaplan Meier survival curve showing the overall survival for the gene AURKA, split in three parts with quantile() at 33% and 66%
################
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
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.4,1))


################
## chunk 33: create table with p-values that show how significant the difference between two concordance indices of two genes is.
################
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
xtable(ccmData, digits=4)

################
## chunk 34: session info
################
toLatex(sessionInfo())