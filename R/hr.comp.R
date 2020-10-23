#' @title Function to statistically compare two hazard ratios
#'
#' @description
#' This function compares two hazard ratios from their betas and standard
#'   errors as computed by a Cox model for instance. The statistical test is
#'   a Student t test for dependent samples. The two hazard ratios must be
#'   computed from the same survival data.
#'
#' @usage hr.comp(hr1, hr2)
#'
#' @param hr1 first hazard ratio.
#' @param hr2 second hazard ratio.
#'
#' @details
#' The two hazard ratios must be computed from the same samples (and
#'   corresponding survival data). The function uses a Student t test for
#'   dependent samples.
#'
#' @return
#' A list with items:
#' - p.value: p-value from the Student t test for the comparison beta1 > beta2
#' (equivalently hr1 > hr2)
#' - hr1: value of the first hazard ratio
#' - hr2: value of the second hazard ratio
#'
#' @references
#' Student (1908). "The Probable Error of a Mean", Biometrika, 6, 1, pages 1–25.
#'
#' @seealso
#' [survival::coxph], [stats::t.test]
#'
#' @examples
#' set.seed(12345)
#' age <- as.numeric(rnorm(100, 50, 10) >= 50)
#' size <- as.numeric(rexp(100,1) > 1)
#' stime <- rexp(100)
#' cens <- runif(100,.5,2)
#' sevent <- as.numeric(stime <= cens)
#' stime <- pmin(stime, cens)
#' hr1 <- hazard.ratio(x=age, surv.time=stime, surv.event=sevent)
#' hr2 <- hazard.ratio(x=size, surv.time=stime, surv.event=sevent)
#' hr.comp(hr1=hr1, hr2=hr2)
#'
#' @md
#' @export
hr.comp <-
    function(hr1, hr2) {
        if(hr1$n != hr2$n) { stop("the hazard ratios are computed from different number of samples!") }
        n <- hr1$n
        x1 <- hr1$data$x
        x2 <- hr2$data$x
        beta1 <- hr1$coef
        beta2 <- hr2$coef
        se1 <- hr1$se
        se2 <- hr2$se
        r <- cor(x1, x2, method="spearman", use="complete.obs")
        if(abs(r) < 1) {
            t.stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2 - 2 * r * se1 * se2)
            diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
        } else { diff.ci.p <- 1 }
        return(list("p.value"=diff.ci.p, "hr1"=exp(beta1), "hr2"=exp(beta2)))
    }

