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
library(genefu)
library(Biobase)

##########
## chunk 4: load sample dataset
##########
load(sampleSurvivalData)

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
cindexall.mainz.small <- t(apply(X=exprs(mainz.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz.small)[ ,"t.dmfs"], z=pData(mainz.small)[ ,"e.dmfs"]))

cindexall.transbig.small <- t(apply(X=exprs(transbig.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbig.small)[ ,"t.dmfs"], z=pData(transbig.small)[ ,"e.dmfs"]))

cindexall.vdx.small <- t(apply(X=exprs(vdx.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdx.small)[ ,"t.dmfs"], z=pData(vdx.small)[ ,"e.dmfs"])) 

cindexall.upp.small <- t(apply(X=exprs(upp.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(upp.small)[ ,"t.rfs"], z=pData(upp.small)[ ,"e.rfs"]))

cindexall.unt.small <- t(apply(X=exprs(unt.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(unt.small)[ ,"t.dmfs"], z=pData(unt.small)[ ,"e.dmfs"]))

cindexall.nki.small <- t(apply(X=exprs(nki.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nki.small)[ ,"t.dmfs"], z=pData(nki.small)[ ,"e.dmfs"]))


##########
## chunk 7: compute D.indices for all datasets
##########
dindexall.mainz.small <- t(apply(X=exprs(mainz.small), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz.small)[ ,"t.dmfs"], z=pData(mainz.small)[ ,"e.dmfs"]))

dindexall.transbig.small <- t(apply(X=exprs(transbig.small), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbig.small)[ ,"t.dmfs"], z=pData(transbig.small)[ ,"e.dmfs"]))

dindexall.upp.small <- t(apply(X=exprs(upp.small), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(upp.small)[ ,"t.rfs"], z=pData(upp.small)[ ,"e.rfs"]))

dindexall.unt.small <- t(apply(X=exprs(unt.small), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(unt.small)[ ,"t.dmfs"], z=pData(unt.small)[ ,"e.dmfs"]))

dindexall.vdx.small <- t(apply(X=exprs(vdx.small), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdx.small)[ ,"t.dmfs"], z=pData(vdx.small)[ ,"e.dmfs"]))

dindexall.nki.small <- t(apply(X=exprs(nki.small), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nki.small)[ ,"t.dmfs"], z=pData(nki.small)[ ,"e.dmfs"]))


##########
## chunk 8: compute hazard.ratios for all datasets
##########
hratio.mainz.small <- t(apply(X=exprs(mainz.small), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz.small)[ ,"t.dmfs"], z=pData(mainz.small)[ ,"e.dmfs"]))

hratio.transbig.small <- t(apply(X=exprs(transbig.small), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbig.small)[ ,"t.dmfs"], z=pData(transbig.small)[ ,"e.dmfs"]))

hratio.upp.small <- t(apply(X=exprs(upp.small), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(upp.small)[ ,"t.rfs"], z=pData(upp.small)[ ,"e.rfs"]))

hratio.unt.small <- t(apply(X=exprs(unt.small), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(unt.small)[ ,"t.dmfs"], z=pData(unt.small)[ ,"e.dmfs"]))

hratio.vdx.small <- t(apply(X=exprs(vdx.small), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdx.small)[ ,"t.dmfs"], z=pData(vdx.small)[ ,"e.dmfs"]))

hratio.nki.small <- t(apply(X=exprs(nki.small), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nki.small)[ ,"t.dmfs"], z=pData(nki.small)[ ,"e.dmfs"]))


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
rownames(tt) <- gs.list
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
rownames(tt) <- gs.list
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
rownames(tt) <- gs.list
colnames(tt) <- c("hratio","hratio.se","lower","upper")
chratio <- tt


##########
## chunk 12: forestplot for each gene combining all datasets
##########
labeltext <- cbind(c("Gene Symbol",gs.list),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
r.mean <- c(NA,ccindex$cindex)
r.lower <- c(NA,ccindex$cindex + qnorm(0.025, lower.tail=TRUE) * ccindex$cindex.se)
r.upper <- c(NA,ccindex$cindex + qnorm(0.025, lower.tail=FALSE) * ccindex$cindex.se)

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.3,0.7,0.05), xlab=paste("concordance index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.3,1))
title(paste("cindex forestplot, each gene, ccindex"))

##########
## chunk 13: forestplot for each gene vombining all datasets
##########
labeltext <- cbind(c("Gene Symbol",gs.list),c(rep(mybigspace,8)),c(rep(mybigspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
tt <- log2(cdindex)

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1.5,0.1), xlab=paste("log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,1.5))
title(paste("log2 dindex forestplot, each gene, cdindex"))


##########
## chunk 14: forestplot for each gene combining all datasets
##########
labeltext <- cbind(c("Gene Symbol",gs.list),c(rep(mybigspace,8)),c(rep(mybigspace,8)))
bs <- rep(0.5, nrow(labeltext))                              
tt <- log2(chratio)

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$hratio), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1,0.1), xlab=paste("log2 hazard.ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,1))
title(paste("log2 hratio forestplot, each gene, chratio"))


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

rownames(tt) <- dataset.list
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",dataset.list),c(rep(mybigspace,8)))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA cindex", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.5,1), is.summary=(c(rep(FALSE,7),TRUE)))
title(paste("cindex forestplot, AURKA"))


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

rownames(tt) <- dataset.list
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",dataset.list),c(rep(mybigspace,8)))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,2.5,0.5), xlab=paste("AURKA log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0,2.5), is.summary=(c(rep(FALSE,7),TRUE)))
title(paste("log2 dindex forestplot, AURKA"))


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

rownames(tt) <- dataset.list
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",dataset.list),c(rep(myspace,8)))
bs <- rep(0.5, nrow(labeltext))
   
forestplot.surv(labeltext=labeltext, mean=c(NA,tt$hratio), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,3.5,0.5), xlab=paste("AURKA log2 hazard.ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,3.5),is.summary=(c(rep(FALSE,7),TRUE)))
title(paste("log2 hratio forestplot, AURKA"))


##########
## chunk 18: forestplot for every single gene over all datasets
##########
for(i in 1:length(gs.list)) {
  par(mfrow=c(3,3))
  myspace <- " "
  tt <- rbind(cindexall.mainz.small[i,],
              cindexall.transbig.small[i,],
              cindexall.upp.small[i,],
              cindexall.unt.small[i,],
              cindexall.vdx.small[i,],
              cindexall.nki.small[i,],
              as.numeric(ccindex[i,]))

  rownames(tt) <- dataset.list
  tt <- as.data.frame(tt)
  labeltext <- cbind(c("Dataset",dataset.list),c(rep(myspace,8)))
  bs <- rep(0.5, nrow(labeltext))

  r.mean <- c(NA,tt$cindex)
  r.lower <- c(NA,tt$cindex + qnorm(0.025, lower.tail=TRUE) * tt$cindex.se)
  r.upper <- c(NA,tt$cindex + qnorm(0.025, lower.tail=FALSE) * tt$cindex.se)
  
  x.ticks.lower <- (floor((min(r.mean,na.rm=TRUE) - 0.1) * 10)/10)
  x.ticks.upper <- (floor((max(r.mean,na.rm=TRUE) + 0.2) * 10)/10)
   
  forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(x.ticks.lower,x.ticks.upper,0.05), xlab=paste(gs.list[i], myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.3,0.8),is.summary=(c(rep(FALSE,7),TRUE)))
  title(paste("cindex forestplot, ", gs.list[i]))
}

##########
## chunk 19: kaplan meier survival curve for all datasets
##########
surv.data <- censor.time(surv.time=c(pData(mainz.small)[ ,"t.dmfs"], pData(transbig.small)[ ,"t.dmfs"], pData(unt.small)[ ,"t.dmfs"], pData(vdx.small)[ ,"t.dmfs"], pData(upp.small)[ ,"t.rfs"], pData(nki.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz.small)[ ,"e.dmfs"], pData(transbig.small)[ ,"e.dmfs"], pData(unt.small)[ ,"e.dmfs"], pData(vdx.small)[ ,"e.dmfs"], pData(upp.small)[ ,"e.rfs"], pData(nki.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainz.small))), rep("transbig", nrow(pData(transbig.small))), rep("unt", nrow(pData(unt.small))), rep("vdx", nrow(pData(vdx.small))), rep("upp", nrow(pData(upp.small))), rep("nki", nrow(pData(nki.small)))), levels=c("mainz", "transbig", "unt", "vdx", "upp", "nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS(upp)", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange", "red", "darkblue", "darkgreen", "black", "brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, all datasets"))


##########
## chunk 20: kaplan meier survival curve for all six datasets
##########
par(mfrow=c(3,3))
surv.data <- censor.time(surv.time=c(pData(mainz.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainz.small)))), levels=c("mainz"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, MAINZ"))

surv.data <- censor.time(surv.time=c(pData(transbig.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(transbig.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("transbig", nrow(pData(transbig.small)))), levels=c("transbig"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("red"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, TRANSBIG"))

surv.data <- censor.time(surv.time=c(pData(upp.small)[ ,"t.rfs"]) / 365, surv.event=c(pData(upp.small)[ ,"e.rfs"]), time.cens=tc / 365)
gg <- factor(c(rep("upp", nrow(pData(upp.small)))), levels=c("upp"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of RFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE) 
title(paste("km survival curve, UPP"))

surv.data <- censor.time(surv.time=c(pData(unt.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(unt.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("unt", nrow(pData(unt.small)))), levels=c("unt"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("black"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, UNT"))

surv.data <- censor.time(surv.time=c(pData(vdx.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(vdx.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("vdx", nrow(pData(vdx.small)))), levels=c("vdx"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkblue"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, VDX"))

surv.data <- censor.time(surv.time=c(pData(nki.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(nki.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("nki", nrow(pData(nki.small)))), levels=c("nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkgreen"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)      
title(paste("km survival curve, NKI"))


##########
## chunk 21: kaplan meier survival plot for AURKA combining 6 datasets
##########
aurkaGs <- "AURKA"
aurkaGid <- 6790
aurkaPaf <- "208079_s_at"
aurkaPagi <- "NM_003600"

surv.time.all <- c(pData(mainz.small)[ ,"t.dmfs"], pData(transbig.small)[ ,"t.dmfs"], pData(unt.small)[ ,"t.dmfs"], pData(upp.small)[ ,"t.rfs"], pData(vdx.small)[ ,"t.dmfs"], pData(nki.small)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz.small)[ ,"e.dmfs"], pData(transbig.small)[ ,"e.dmfs"], pData(unt.small)[ ,"e.dmfs"], pData(upp.small)[ ,"e.rfs"], pData(vdx.small)[ ,"e.dmfs"], pData(nki.small)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz)[aurkaPaf,], exprs(transbig)[aurkaPaf,], exprs(unt)[aurkaPaf,], exprs(upp)[aurkaPaf,], exprs(vdx)[aurkaPaf,], exprs(nki)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz)[aurkaPaf,]), length(exprs(transbig)[aurkaPaf,]), length(exprs(unt)[aurkaPaf,]), length(exprs(upp)[aurkaPaf,]), length(exprs(vdx)[aurkaPaf,]), length(exprs(nki)[aurkaPagi,]))

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
gg <- factor(c(rep("mainz", nrow(pData(mainz.small))), rep("transbig", nrow(pData(transbig.small))), rep("unt", nrow(pData(unt.small))), rep("upp", nrow(pData(upp.small))), rep("vdx", nrow(pData(vdx.small))), rep("nki", nrow(pData(nki.small)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Survival", main.title="", sub.title=NULL, leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE)


##########
## 
##########


##########
## 
##########



##########
## 
##########


##########
## 
##########



