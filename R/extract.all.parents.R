#' .extract.all.parents
#'
#' Function taking the output of the regrnet.ensemble method and returns a matrix
#' containing one equivalent model in each column. The target variable is in the first row.
#'
#' @param res.main output of regrnet.ensemble
#' @param maxparents maxparents parameter of the netinf method
#' @param predn list of target variables for which ensemble method was run.
#'
#' @md
#' @keywords internal
#' @noRds
.extract.all.parents  <-
function(data,res.main,maxparents,predn) {


	final <- NULL
	cnt_main <- 1
	for(imain in 1:length(predn)){
		res.vec <- res.main[cnt_main:(cnt_main+2*res.main[cnt_main])]
		if(length(res.vec)>3){
			cnt_main <- cnt_main+2*res.main[cnt_main]+1
			nsol <- sum(res.vec==0)
			res <- matrix(0,ncol=nsol,nrow=(maxparents+1))

			val <- res.vec[2:(res.vec[1]+1)]
			ind <- res.vec[(res.vec[1]+2):(2*res.vec[1]+1)]
			res[1,1:ind[1]] <- rep(val[1],ind[1])
			nvar <- length(val)
			level <- 2
			cnt <- 1
			nelem <- 0
			sum_old <- sum(ind[1])
			last_level <- FALSE
			i <- 2
			ind2 <- 1

			while(i<=nvar && !last_level){
				if(ind[i]!=0){
					if(ind[i]>1){
						tmp <- res[,(ind2+1):nsol]
						res[level,ind2] <- val[i]
						for(j in 1:level){
							res[j,(ind2+1):(ind2+ind[i]-1)] <- rep(res[j,ind2],(ind[i]-1))
						}
						if((nsol-(ind2+ind[i]-1))>0){
							res <- cbind(res[,1:(ind2+ind[i]-1)],tmp[,1:(nsol-(ind2+ind[i]-1))])
						}else{
							res <- res[,1:(ind2+ind[i]-1)]
						}
					}
					res[level,ind2:(ind2-1+ind[i])] <- rep(val[i],ind[i])
					ind2 <- ind2+ind[i]
				}else{
					res[level,ind2] <- val[i]
					ind2 <- ind2+1
				}
				nelem <- nelem+ind[i]
				if(cnt==sum_old){
					sum_old <- nelem
					if(nelem==0){
						last_level <- TRUE
					}
					nelem <- 0
					cnt <- 1
					level <- level+1
					ind2 <- 1
				}
				else if(cnt< sum_old){
					cnt <- cnt+1
				}
				i <- i+1
			}
			final <- cbind(final,res)
		}
	}
	dimension <- dim(final)
	final <- colnames(data)[final]
	dim(final) <- dimension
	return(final)
}
