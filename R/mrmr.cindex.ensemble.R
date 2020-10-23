#' @name mrmr.cindex.ensemble
#'
#' @title Function to compute the concordance index for survival or binary class
#'   prediction
#'
#' @description
#' Function to compute the minimum redundancy - maximum relevance (mRMR) ranking
#'   for a risk prediction or a binary classification task based on the
#'   concordance index. The mRMR feature selection has been adapted to use the
#'   concordance index to estimate the correlation between a variable and the
#'   output (binary or survival) data.
#'
#' @usage
#' mrmr.cindex.ensemble(x, surv.time, surv.event, cl, weights, comppairs=10,
#' strat, alpha = 0.05, outx = TRUE, method = c("conservative", "noether",
#' "nam"), alternative = c("two.sided", "less", "greater"), maxparents,
#' maxnsol, nboot = 200, na.rm = FALSE)
#'
#' @param x
#' @param surv.time
#' @param surv.event
#' @param cl
#' @param weights
#' @param comppairs
#' @param strat
#' @param alpha
#' @param outx
#' @param method
#' @param alternative
#' @param maxparents
#' @param maxnsol
#' @param nboot
#' @param na.rm
#'
#' @return A mRMR ranking
#'
#' @section Note:
#' The "direction" of the concordance index (< 0.5 or > 0.5) is the opposite
#'   than the rcorr.cens function in the Hmisc package. So you can easily get
#'   the same results than rcorr.cens by changing the sign of x.
#'
#' @references
#' Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996) "Tutorial in
#'   biostatistics: multivariable prognostic models: issues in developing
#'   models, evaluating assumptions and adequacy, and measuring and reducing
#'   error", Statistics in Medicine, 15, pages 361–387.
#'
#' Pencina, M. J. and D'Agostino, R. B. (2004) "Overall C as a measure of
#'   discrimination in survival analysis: model specific population value and
#'   confidence interval estimation", Statistics in Medicine, 23, pages
#'   2109–2123, 2004.
#'
#' @seealso
#' [Hmisc::rcorr.cens], [CPE::phcpe], [clinfun::coxphCPE]
#'
#' @md
#' @export
mrmr.cindex.ensemble <-
    function(x, surv.time, surv.event, cl, weights, comppairs=10, strat, alpha=0.05, outx=TRUE, method=c("conservative", "noether", "nam"), alternative=c("two.sided", "less", "greater"), maxparents, maxnsol, nboot=200, na.rm=FALSE) {

        nvar<-ncol(x)
        nsample<-nrow(x)

        vec_ensemble<-.Call(.C_mrmr_cIndex_ensemble_remove,data.matrix(x),as.integer(is.na(x)),maxparents,ncol(x),nrow(x),1,1,nboot,maxnsol,-1000,as.integer(as.logical(TRUE)),as.integer(sort(unique(strat))),as.integer(cl ),as.double(surv.time),as.integer(surv.event),as.double(weights),as.integer(strat),as.integer(sum(weights)),as.integer(as.logical(outx)),as.integer(length(strat)),as.integer(length(sort(unique(strat)))))

        vec_ensemble[2:vec_ensemble[1]+1]<-vec_ensemble[2:vec_ensemble[1]+1]-1

        models.equiv <- .extract.all.parents(x,vec_ensemble,maxparents,1)
        models.equiv[1,]<-rep("T",ncol(models.equiv))

        return(models.equiv)
    }


