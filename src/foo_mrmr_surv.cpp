#include "mrmr_surv.h"

double get_correlation(double data [],int namat[], int ind_x, int ind_y, int size){
	//compute correlation of two variables;
	//data: contains all data in a vector; variable-wise appended
	//ind_x: starting index of first variable in data
	//ind_y: starting index of second variable in data
	//size: number of samples for both variables
	double mean_data_x=0.0,mean_data_y=0.0;
	double correlation_nom=0.0,correlation_den_x=0.0,correlation_den_y=0.0;
	
	for( unsigned int i=0; i< size; ++i ) {
		if (namat[ind_x+i]==0 && namat[ind_y+i]==0 ) {
			mean_data_x+=data[ind_x+i];
			mean_data_y+=data[ind_y+i];
		}
	}
	
	mean_data_x=mean_data_x/size;
	mean_data_y=mean_data_y/size;
	
	for( unsigned int i=0; i< size; ++i ) {		
		if(namat[ind_x+i]==0 && namat[ind_y+i]==0){
		correlation_nom+=(data[ind_x+i]-mean_data_x)*(data[ind_y+i]-mean_data_y);
		correlation_den_x+=(data[ind_x+i]-mean_data_x)*(data[ind_x+i]-mean_data_x);
		correlation_den_y+=(data[ind_y+i]-mean_data_y)*(data[ind_y+i]-mean_data_y);
		}
	}
	return correlation_nom/(sqrt(correlation_den_x*correlation_den_y));
}


void build_mim_subset(double mim[],double data[], int namat [],int nvar,int nsample, int subset [],int size_subset){
	//compute mutual information matrix
	//mim:			matrix (stored as vector) in which the mi values will be stored
	//data:			contains all data in a vector; variable-wise appended
	//nvar:			number of variables
	//nsample:		number of samples in dataset
	//subset:		indices of samples to be included in the bootstrapping data
	//size_subset:	number of variables in the bootstrapped dataset
	
	double tmp;
	double *data_x;
	int *namat_x;
	
	namat_x = (int*) R_alloc(nvar*size_subset, sizeof(int));
	data_x = (double *) R_alloc(nvar*size_subset, sizeof(double));
	
	for(unsigned int i=0; i< size_subset; ++i){
		for(unsigned int j=0; j< nvar; ++j){
			data_x[size_subset*j+i]=data[(subset[i])+nsample*j];
			namat_x[size_subset*j+i]=namat[(subset[i])+nsample*j];
		}
	}
	
	for(unsigned int i=0; i< nvar; ++i){
		mim[i*nvar+i]=0;
		for(unsigned int j=i+1; j< nvar; ++j){
			tmp=get_correlation(data_x,namat_x,i*size_subset,j*size_subset,size_subset);
			tmp=tmp*tmp;
			if(tmp>0.999999){
				tmp=0.999999;
			}
			mim[j*nvar+i]= -0.5* log (1-tmp);
			mim[i*nvar+j]=mim[j*nvar+i];
		}
	}
	//delete [] namat_x;
	//delete [] data_x;

}

double returnConcordanceIndexC(int *msurv, int *ustrat, double *x2, int *cl2,
					   double *st, int *se, double *weights, int *strat, int *N, int *outx, int *lenS, int *lenU)
{
	
	int lenUstrat = *lenU;
	int lenStrat = *lenS;
	
	double res_ch[lenStrat];
	double res_dh[lenStrat];

	double res_cIndex=0;
	
	int Ns_old = 0;
	int Ns = 0;
	for(int s=0; s < lenUstrat; s++) {
		int ixs[lenStrat];
		for(int i =0; i < lenStrat; i++){
			ixs[i] = 0;
			if(strat[i] == ustrat[s]){
				ixs[i] = 1;
			} else {
				ixs[i] = 0;
			}
		}
		Ns_old += Ns;
		Ns = 0;
		for(int i=0; i < lenStrat; i++){
			if(ixs[i] == 1){
				Ns++;
			}
		}
		double xs[Ns];
		int c = 0;
		for(int i=0; i < lenStrat; i++){
			if(ixs[i] == 1){
				xs[c] = x2[i];
				c++;
			}
		}
		int cls[Ns];
		c = 0;
		for(int i=0; i < lenStrat; i++){
			if(ixs[i] == 1){
				cls[c] = cl2[i];
				c++;
			}
		}
		double sts[Ns];
		c = 0;
		for(int i=0; i < lenStrat; i++){
			if(ixs[i] == 1){
				sts[c] = st[i];
				c++;
			}
		}        
		int ses[Ns];
		c = 0;
		for(int i=0; i < lenStrat; i++){
			if(ixs[i] == 1){
				ses[c] = se[i];
				c++;
			}
		}
		double weightss[Ns];
		c = 0;
		for(int i=0; i < lenStrat; i++){
			if(ixs[i] == 1){
				weightss[c] = weights[i];
				c++;
			}
		}
		double chs[Ns];
		double dhs[Ns];
		double uhs[Ns];
		double rphs[Ns];
		for (int h=0; h < Ns; h++) {       	
			double chsj, dhsj, uhsj, rphsj = 0;
			for (int j=0; j < Ns; j++) {
				double whj = weightss[h] * weightss[j];
				if((*msurv == 1 && (sts[h] < sts[j] && ses[h] == 1)) || (*msurv == 0 && cls[h] > cls[j])){
					rphsj = rphsj + whj;
					if (xs[h] > xs[j]) {  
						chsj = chsj + whj;
					} else if (xs[h] < xs[j]) {
						dhsj = dhsj + whj;
					}
					else {
						if (*outx == 1) {
							uhsj = uhsj + whj;
						} else {
							dhsj = dhsj + whj;
						}
					}
				}    
				if((*msurv == 1 && (sts[h] > sts[j] && ses[j] == 1)) || (*msurv == 0 && cls[h] < cls[j])){
					rphsj = rphsj + whj;
					if (xs[h] < xs[j]) {
						chsj = chsj + whj;
					}
					else if (xs[h] > xs[j]) {
						dhsj = dhsj + whj;
					}
					else {
						if (*outx == 1) {
							uhsj = uhsj + whj;
						} else {
							dhsj = dhsj + whj;
						}
					}
				}
			}
			chs[h] = chsj;
			dhs[h] = dhsj;
			uhs[h] = uhsj;
			rphs[h] = rphsj;
			chsj = 0;
			dhsj = 0;
			uhsj = 0;
			rphsj = 0;
		}
		for(int i = 0; i < Ns; i++){
			int pos = i + Ns_old;			
			res_ch [pos] =chs[i];
			res_dh [pos] =dhs[i];
		}
			
	}

	double tmp_ch=0, tmp_dh=0;
	for(int s=0; s < lenStrat; s++) {
		tmp_ch+=res_ch[s];
		tmp_dh+=res_dh[s];
	}
	
	double n=*N;
	tmp_ch=(1/(n *(n - 1))) * tmp_ch;
	tmp_dh=(1/(n *(n - 1))) * tmp_dh;
	res_cIndex=tmp_ch/ (tmp_ch+tmp_dh);
	
	/// scale value to be in the intervall [-1,1] as correlation and square to be on same scale as mutual information [0,1]
	
	res_cIndex=2*res_cIndex-1;
	res_cIndex=res_cIndex*res_cIndex;
	//cout<<"res cIndex "<<res_cIndex<<endl;
	return res_cIndex;
}


SEXP get_concordanceIndex_onevariable(SEXP Rmsurv, SEXP Rustrat, SEXP Rx2,SEXP Rcl2, SEXP Rst, SEXP Rse, SEXP Rweights, SEXP Rstrat, SEXP RN, SEXP Routx, SEXP RlenS, SEXP RlenU){
	double *res_cI , res_cIndex;
	
	double *x2, *st, *weights;
	int *msurv, *ustrat, *cl2, *se, *strat, *N, *outx, *lenS, *lenU;
	
	SEXP Rres;
	
	msurv=INTEGER_POINTER(Rmsurv);
	ustrat=INTEGER_POINTER(Rustrat);
	x2=NUMERIC_POINTER(Rx2);
	cl2=INTEGER_POINTER(Rcl2);
	st=NUMERIC_POINTER(Rst);
	se =INTEGER_POINTER(Rse);
	weights=NUMERIC_POINTER(Rweights);
	strat=INTEGER_POINTER(Rstrat);
	N =INTEGER_POINTER(RN);
	outx=INTEGER_POINTER(Routx);

	lenS=INTEGER_POINTER(RlenS);
	lenU=INTEGER_POINTER(RlenU);
		

	PROTECT(Rres = NEW_NUMERIC(1));
	res_cI=NUMERIC_POINTER(Rres);
	

	res_cIndex=returnConcordanceIndexC(msurv, ustrat, x2, cl2,st, se, weights,strat, N, outx, lenS, lenU) ;

	res_cI[0]=res_cIndex;

	UNPROTECT(1);
	
	return Rres;
}

SEXP mrmr_cIndex(SEXP Rdata, SEXP Rnamat, SEXP RcIndex, SEXP Rnvar, SEXP Rnsample, SEXP Rthreshold){     
	double *data, *cIndex;
	double *rel, *red, *res, *mim, score=1,*threshold, *res_final;
	
	const int *nvar, *nsample;
	int *ind, *namat;
	
	unsigned int n, jmax=0;
	PROTECT(Rdata = AS_NUMERIC(Rdata));
	PROTECT(Rnamat = AS_INTEGER(Rnamat));
	PROTECT(RcIndex = AS_NUMERIC(RcIndex));
	PROTECT(Rnvar= AS_INTEGER(Rnvar));
	PROTECT(Rnsample= AS_INTEGER(Rnsample));
	PROTECT(Rthreshold= AS_NUMERIC(Rthreshold));
	
	
	data = NUMERIC_POINTER(Rdata);
	namat=INTEGER_POINTER(Rnamat);
	cIndex= NUMERIC_POINTER(RcIndex);
	nvar = INTEGER_POINTER(Rnvar);
	nsample = INTEGER_POINTER(Rnsample);
	threshold = NUMERIC_POINTER(Rthreshold);

	n = *nvar;
	
	//new variables
	SEXP Rmim, Rres,Rred,Rrel,Rind,Rres_final;
		
	PROTECT(Rmim = NEW_NUMERIC(n*n));
	PROTECT(Rres = NEW_NUMERIC(n));
	PROTECT(Rres_final = NEW_NUMERIC(n));
	PROTECT(Rrel = NEW_NUMERIC(n));
	PROTECT(Rred = NEW_NUMERIC(n)); 
	PROTECT(Rind = NEW_INTEGER(*nsample));
	
	ind = INTEGER_POINTER(Rind);	
	mim = NUMERIC_POINTER(Rmim);
	res = NUMERIC_POINTER(Rres);
	rel = NUMERIC_POINTER(Rrel);
	red = NUMERIC_POINTER(Rred);
	res_final = NUMERIC_POINTER(Rres_final);
	

	for(unsigned int i=0;i < *nsample; ++i){
		ind[i]=i;
	}

	build_mim_subset(mim, data, namat, n, *nsample, ind, *nsample);
	
	
	for( unsigned int i=0; i< n; ++i ){
			res[i]=*threshold;
			res_final[i]=*threshold;
	}
	
		//init rel and red and select first
		for( unsigned int j=0; j<n; ++j ) {
			rel[j]=cIndex[j];
			red[j]=0;
			if( rel[j] > rel[jmax])
				jmax = j;
		}

		score = rel[jmax];
		if( res[jmax] < score ) {
			res[jmax] = score;
		}

		//select others
		for(unsigned int k=1; k < n+1; k++ ) { 
			jmax = 0;
	
			for(unsigned int j=1; j < n; ++j ) {
				if( (rel[j] - red[j]/k) > (rel[jmax] - red[jmax]/k) ) 
					jmax = j;
			} 
			score = (rel[jmax] - (red[jmax]/k));
			if( res[jmax] < score ) {
				res[jmax] = score;
			}
			
			//update rel and red
			rel[jmax]=-1000; 
			for( int l=0; l<n; ++l )  
				red[l] += mim[l*n+jmax];
			
			// stop criterion
			if( score < *threshold ) k=n;

		}


	for(unsigned int i=0; i< n; ++i ) {
		res_final[i]=res[i];
	}
	UNPROTECT(12);
	
	return Rres_final;
}
