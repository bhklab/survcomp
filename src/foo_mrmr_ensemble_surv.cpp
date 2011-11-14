#include "foo_mrmr_ensemble_surv.h"

double get_correlation_ensemble(double data [],int namat[], int ind_x, int ind_y, int size){
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


double returnConcordanceIndexC(int *msurv, int *ustrat, double *x2, int *cl2,
					   double *st, int *se, double *weights, int *strat, int *N, int *outx, int lenS, int *lenU)
{
	
	int lenUstrat = *lenU;
	int lenStrat = lenS;
	
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
//	cout<<"tmp_ch "<<tmp_ch<<endl;
//	cout<<"n "<<n<<endl;
//	cout<<"data [0] "<<x2[0]<<endl;
//	cout<<"msurv "<<*msurv<<endl;
	res_cIndex=2*res_cIndex-1;
	res_cIndex=res_cIndex*res_cIndex;
//	cout<<"res cIndex function "<<res_cIndex<<endl;
	return res_cIndex;
}


void build_mim_cIndex_subset(double mim[],double data[], int namat [],int nvar,int nsample, int subset [],int size_subset,int *msurv, int *ustrat, int *cl2, double *st, int *se, double *weights, int *strat, int *N, int *outx, int *lenU){
	//compute mutual information matrix
	//mim:			matrix (stored as vector) in which the mi values will be stored
	//data:			contains all data in a vector; variable-wise appended
	//nvar:			number of variables
	//nsample:		number of samples in dataset
	//subset:		indices of samples to be included in the bootstrapping data
	//size_subset:	number of variables in the bootstrapped dataset
	
	double tmp;
	double *data_x, *st_x, *weights_x;
	int *namat_x, *msurv_x, *ustrat_x, *cl2_x, *se_x, *strat_x;
	
	namat_x = (int*) R_alloc(nvar*size_subset, sizeof(int));
	cl2_x = (int*) R_alloc(size_subset, sizeof(int));
	se_x = (int*) R_alloc(size_subset, sizeof(int));
	strat_x = (int*) R_alloc(size_subset, sizeof(int));

	data_x = (double *) R_alloc(nvar*size_subset, sizeof(double));
	st_x = (double *) R_alloc(size_subset, sizeof(double));
	weights_x = (double *) R_alloc(size_subset, sizeof(double));

	for(unsigned int i=0; i< size_subset; ++i){
		for(unsigned int j=0; j< (nvar-1); ++j){
			data_x[size_subset*j+i]=data[(subset[i])+nsample*j];
			namat_x[size_subset*j+i]=namat[(subset[i])+nsample*j];
		}

		cl2_x[i]=cl2[subset[i]];
		se_x[i]=se[subset[i]];
		strat_x[i]=strat[subset[i]];
		st_x[i]=st[subset[i]];
		weights_x[i]=weights[subset[i]];
	}
	for(unsigned int i=0; i< nvar-1; ++i){
		mim[(i+1)*(nvar)+(i+1)]=0;
		for(unsigned int j=i+1; j< nvar-1; ++j){
			tmp=get_correlation_ensemble(data_x,namat_x,i*size_subset,j*size_subset,size_subset);
			tmp=tmp*tmp;
			if(tmp>0.999999){
				tmp=0.999999;
			}
			mim[(j+1)*(nvar)+i+1]= -0.5* log (1-tmp);
			mim[(i+1)*(nvar)+j+1]=mim[(j+1)*(nvar)+i+1];
		}
	}
	
	double *data_small;
	data_small =(double*) R_alloc(size_subset, sizeof(double));
	
	for(int j=0;j< nvar-1 ;++j){
		for(int i=0;i< size_subset;++i){
			data_small[i]=data_x[j* (nvar-1)+i];
		}
		mim[j+1]=returnConcordanceIndexC(msurv, ustrat, data_small, cl2_x,st_x, se_x, weights_x,strat_x, N, outx, size_subset, lenU) ;	
		mim[(nvar)*(j+1)] = mim[j+1];
	}

}

void remove_childless_nodes( tree<int>& res, tree<double>&res_mean, int max_elements_tmp){
	tree<int>::pre_order_iterator it_tmp,it_back, it=res.begin();
	tree<double>::pre_order_iterator  it_mean_tmp,it_mean_back, it_mean=res_mean.begin();
	tree<int>::leaf_iterator li;
	tree<double>::leaf_iterator li_mean;
	bool found_child, multiple;
	
	int depth_max=0;
	//determine max depth
	while(it!=res.end()){
		if(depth_max<res.depth(it)){
			depth_max=res.depth(it);
		}
		it++;
	}
	it=res.begin();
	while(it!=res.end() && max_elements_tmp<= (depth_max+1)) {
		if(res.depth(it)<=(max_elements_tmp-2) && res.number_of_children(it)==0){ //advance through the tree
			it_tmp=res.parent(it);it_mean_tmp=res_mean.parent(it_mean);
			found_child=false; multiple=false;
			while(!found_child && (it_tmp!=res.begin() || res.number_of_children(it_tmp)>1)){ //end loop if there is a node with more than one child or if back at top and number of children==1 for top node
				if(res.number_of_children(it_tmp)==1){
					it_back=it_tmp;it_mean_back=it_mean_tmp; //if this was the last level for which there is only one child
					multiple=true;
					if(it_tmp!=res.begin()){
						it_tmp=res.parent(it_tmp); it_mean_tmp=res_mean.parent(it_mean_tmp);
					}else{//in case of having arrived at the top node
						res.erase(it_back); res_mean.erase(it_mean_back);
						found_child=true;
					}		
				}else{
					if(multiple){
						res.erase(it_back); res_mean.erase(it_mean_back);
					}else{
						res.erase(it); res_mean.erase(it_mean);
					}
					found_child=true;
				}
			}
			it=it_tmp; it_mean=it_mean_tmp;
		}else{
			++it; ++it_mean;
		}
	}	
}
void bootstrap_tree(tree<int>& res,tree<double>& res_mrmr, double data[],int namat[], int nsamples,int n, int rep_boot, int *msurv, int *ustrat, double *x2, int *cl2,double *st, int *se, double *weights, int *strat, int *N, int *outx, int lenS, int *lenU){
	int  nsub, *prev_sel,nsamples_boot=nsamples,*to_remove;
	tree<int>::iterator li=res.begin_leaf(),li2;
	tree<double>::iterator li_mrmr=res_mrmr.begin_leaf(),li2_mrmr;
	double *mean, *sd;
	int cnt_leafs=0;
	int max_depth=res.depth(li),index;
	while (li!=res.end()) {
		if(res.depth(li)==max_depth){
			cnt_leafs++;
		}
		li++;
	}
	
	li=res.begin_leaf();	
	
	mean =(double*) R_alloc(cnt_leafs, sizeof(double));
	sd =(double*) R_alloc(cnt_leafs, sizeof(double));
	to_remove=(int*) R_alloc(cnt_leafs, sizeof(int));
	
	for(int k=0;k<cnt_leafs;k++){
		mean[k]=0;sd[k]=0;
	}
	
	int target=*res.begin();
	int nto_remove=0;
	
	prev_sel=(int*) R_alloc(max_depth, sizeof(int));
	
	int k=0;
	while (li!=res.end()) {
		if(res.depth(li)==max_depth){
			li2=li;
			prev_sel[max_depth-1]=*(li);
			li2=res.parent(li2);
			index=max_depth-2;
			while (li2!=res.begin()) {
				prev_sel[index]=*(li2);
				index--;
				li2=res.parent(li2);
			}
			bootstrap_mrmr(mean[k], sd[k], data,namat,n, rep_boot,nsamples_boot,nsamples, target, prev_sel[max_depth-1], max_depth-1,prev_sel, msurv, ustrat, x2, cl2,st, se, weights,strat, N, outx, lenS, lenU);
			k++;
		}
		li++;
	}
	double max_mrmr=-1000;
	int max_mrmr_ind=-1;
	for(int k=0;k<cnt_leafs;k++){
		if(mean[k]>max_mrmr){
			max_mrmr=mean[k];
			max_mrmr_ind=k;
		}
	}
	for(int k=0;k<cnt_leafs;k++){
		if(k!=max_mrmr_ind && (mean[k] < max_mrmr-sd[max_mrmr_ind])){
			to_remove[nto_remove]=k;
			nto_remove++;
		}
	}
	
	int cnt2=nto_remove;
	if(cnt2>0){
		li=res.begin_leaf(res.end());
		sort(to_remove,to_remove+cnt2);
		li_mrmr=res_mrmr.begin_leaf(res_mrmr.end());
		int cnt_back=cnt2;
		while (cnt_leafs>=0 && cnt2>0) {
			li2=li;
			li2_mrmr=li_mrmr;
			li--;li_mrmr--;
			
			while(res.depth(li)<max_depth && li2!=res.begin_leaf(res.begin())) {
				li--;li_mrmr--;
			}
			
			if(to_remove[cnt2-1]==cnt_leafs){
				res.erase(li2);res_mrmr.erase(li2_mrmr);
				cnt2--;
			}
			cnt_leafs--;
		}
	}
	remove_childless_nodes(res, res_mrmr,max_depth+1);
	
}
void bootstrap_mrmr(double &mean, double &sd, double data[],int namat[],int size, int rep_boot, int size_boot,int nsamples, int var_target, int var_interest, int nprev_sel,int* var_ind, int *msurv, int *ustrat, double *x2, int *cl2,
					double *st, int *se, double *weights, int *strat, int *N, int *outx, int lenS, int *lenU)
{
	//mean
	//sd
	//data
	//size
	//rep_boot
	//size_boot
	//nsamples
	//var_target
	//var_interest
	//nprev_sel
	//var_ind
	
	int *ind;
	double *mim, *boot_val, *mat_info;
	
	ind=(int*) R_alloc(size_boot, sizeof(int));
	boot_val =(double*) R_alloc(rep_boot, sizeof(double));
	mat_info =(double*) R_alloc(((size)*(size )), sizeof(double));

	for(unsigned int k=0; k< rep_boot; ++k){
		//in total there will be rep_boot times the mrmr sampled
		//determine the subset of samples that should be used (in total size_boot samples will be selected)
		for(unsigned int i=1;i<= size_boot;++i){
			ind[i-1]=rand () %nsamples;
		}
		// compute mi matrix for the subset
		for( unsigned int i=0; i< size ; ++i ){
			for( unsigned int j=0; j< size ; ++j ){
				mat_info[i+(size )*j]=0;
			}
		}
		
		build_mim_cIndex_subset(mat_info, data, namat, (size) , nsamples, ind, size_boot ,msurv, ustrat, cl2,st, se, weights,strat, N, outx,  lenU);
		boot_val[k]=mrnet_onegene(mat_info, size, nprev_sel, var_ind, var_target, var_interest);	

	}
	
	// determine mean and variance of bootstrapped values
	for(unsigned int i=0;i< rep_boot;++i){
		if(boot_val[i]==boot_val[i]){
			mean+=boot_val[i];
		}
	}
	mean=mean/rep_boot;
	
	for(unsigned int i=0;i< rep_boot;++i){
		if(boot_val[i]==boot_val[i]){
			sd+=(boot_val[i]-mean) * (boot_val[i]-mean);
		}
	}
	sd=sqrt(sd/rep_boot);
	
}	


double mrnet_onegene(double mim [], int size, int nbvar,int *var_ind,int var_target, int var_interest){
	// mim:			mutual information matrix
	// size:		total number of variables
	// nbvar:		number of previously selected variables (not the target)
	// var_ind:		the indices of the previously selected variables as vector
	// var_target:	the index of the target gene
	// var_interest: the variable for which the mrmr score with the target has to be computed; will be used for bootstrapping it
	
	unsigned int jmax;
	double rel, red,res;
	double max_val=-1000;
	
	jmax=var_target-1;
	//initialize the remaining entries to zero
	red=0;
	// the relevance for variable of interest with the target is simply its mutual information with it
	rel=mim[(var_target-1)*size+var_interest-1];
	if(nbvar > 0){
		// in case other variables have been previously selected; compute their redundancy with the variable of interest
		for(unsigned int j=0;j< nbvar; j++){
			red+=mim[(var_ind[j]-1)*size+var_interest-1];
		}		
		res=rel-red/nbvar ;
	}else{
		res=rel;
	}
	return res;
}

int power(int a, int b)
{
	int c=a;
	for (int n=b; n>1; n--) c*=a;
	return c;
}
int verify_equivalentset_nparents (tree<int>& tr, tree<int>::pre_order_iterator it, tree<int>::pre_order_iterator end,tree<double>& tr_mrmr, int maxnsol){
	if(!tr.is_valid(it)) return 0;
	
	bool found=false;
	int number_elements_to_remove=0, cnt=1, index=0;
	tree<int>::leaf_iterator li=tr.begin_leaf(it), li_tmp=tr.begin_leaf(tr.begin()), li_tmp2=li_tmp;
	tree<double>::leaf_iterator li_mrmr=tr_mrmr.begin_leaf(tr_mrmr.begin()), li_mrmr2=li_mrmr;
	int depth=tr.depth(li);
	tree<int>::pre_order_iterator it2;
	int vec_old[depth+1];
	int mat_res [power((maxnsol+1),(depth))][depth+2];
	int number_leafs=power((maxnsol+1),(depth));
	int to_remove[number_leafs];
	
	for (int k=0; k< number_leafs ; k++) {
		to_remove[k]=0;
	}
	int cnt2=0,cnt_leafs=0;
	
	while( li!=tr.end_leaf(it) ){
		vec_old[0]=*(li);
		it2=li;
		
		while(it2!=tr.begin()){
			it2=tr.parent(it2);
			vec_old[cnt]=*(it2);
			cnt++;
		}
		
		sort(vec_old,vec_old+depth+1);
		mat_res[index][depth+1]=0;
		
		for(int k=0;k<=depth;k++){
			mat_res[index][k]=vec_old[k];
			mat_res[index][depth+1]+=vec_old[k]+power(2,k);
		}
		index++;
		cnt=1;
		li++;
		cnt_leafs++;
	}
	
	
	index=0;
	bool found1=false,found2;
	
	for(int k=0;k<(cnt_leafs-1) && !found1;k++){
		for(int j=k+1;j<(cnt_leafs);j++){
			found2=false;
			if(mat_res[k][depth+1]==mat_res[j][depth+1]){
				for(int i=0;i<=depth && !found2;i++){
					if(mat_res[j][i]!=mat_res[k][i]){
						found2=true;
					}
				}
			}else{
				found2=true;
			}
			if(!found2){
				number_elements_to_remove++;
				int tmp=j;
				while(tmp>0){
					li_tmp++,li_mrmr++;
					tmp--;
				}
				
				
				tmp=k;
				while(tmp>0){
					li_tmp2++,li_mrmr2++;
					tmp--;
				}
				
				if(*(li_mrmr)< *(li_mrmr2)){
					to_remove[cnt2]=j;
				}else {
					to_remove[cnt2]=k;
				}
				cnt2++;
			}
		}
	}
	
	if(cnt2>0){
		li=tr.begin_leaf(tr.end());
		sort(to_remove,to_remove+cnt2);
		li_mrmr=tr_mrmr.begin_leaf(tr_mrmr.end());
		int cnt_back=cnt2;
		while (cnt_leafs>=0 && cnt2>0) {
			it2=li;li_mrmr2=li_mrmr;
			li--;li_mrmr--;
			while(to_remove[cnt2-1]==to_remove[cnt2] && cnt2!=cnt_back){
				cnt2--;
			}
			if(to_remove[cnt2-1]==cnt_leafs){
				tr.erase(it2);
				tr_mrmr.erase(li_mrmr2);
				cnt2--;
			}
			cnt_leafs--;
		}
	}
	return number_elements_to_remove;
}


void mrmr_ensemble_one_gene_remove (tree<int>& res, tree<int>::pre_order_iterator one, double data[], int namat[], int nsamples,int n , int max_elements, int predn , int rep_boot, int maxnsol, double threshold, int *msurv, int *ustrat, double *x2, int *cl2,
			double *st, int *se, double *weights, int *strat, int *N, int *outx, int lenS, int *lenU){
	//n					number of variables
	//predn:			index of target node
	
	// nsub: the variables which have been previously selected + target; prev_sel=nsub-target
	// number of samples to use for bootstrapping is equal to total number of samples

	int  *nsub, *prev_sel,nsamples_boot=nsamples, tmp_val_max_ind, *prev_sel_tmp,*vec_sol_local,ndelete; 
	double *vec_mean, *vec_sort, *vec_sd,  *vec_local_max_mean, *vec_local_max_sd,tmp_val_max, *mrmr_vec_sort,*vec_sol_local_mrmr;
	double *mat_info;
	
	
	int cnt=0, max_elements_tmp=1; //current depth in the tree
	int *ind;
	vec_mean =(double*) R_alloc(n, sizeof(double));
	vec_sd =(double*) R_alloc(n, sizeof(double));
	mrmr_vec_sort =(double*) R_alloc(n, sizeof(double));
	vec_local_max_mean =(double*) R_alloc(max_elements, sizeof(double));
	vec_local_max_sd =(double*) R_alloc(max_elements, sizeof(double));
	mat_info =(double*) R_alloc((n*n), sizeof(double));
	ind=(int*) R_alloc(nsamples, sizeof(int));

	for(unsigned int i=1;i<= nsamples;++i){
		ind[i-1]=i-1;
	}
	for( unsigned int i=0; i< n; ++i ){
		for( unsigned int j=0; j< n; ++j ){
			mat_info[i+(n)*j]=0;
		}
	}

	build_mim_cIndex_subset(mat_info, data, namat, n, nsamples, ind, nsamples ,msurv, ustrat, cl2,st, se, weights,strat, N, outx,  lenU);

	for(unsigned int k=0;k< max_elements ;++k){
		vec_local_max_mean[k]=-1000;
	}
	
	prev_sel=(int*) R_alloc(max_elements, sizeof(int));
	nsub=(int*) R_alloc(max_elements, sizeof(int));
	
	tree<int> res_tmp_new=res ;
	tree<int>::iterator it_local=res_tmp_new.begin(),it_local2=it_local;
	
	vec_sol_local=(int*) R_alloc(maxnsol, sizeof(int));
	vec_sol_local_mrmr=(double*) R_alloc(maxnsol, sizeof(double));
	//mrmr score should not be predicted for the target node
	vec_mean[predn-1]=-1000; vec_sd[predn-1]=-1000;
	prev_sel[0]=0; nsub[0]=predn;
	
	tree<double> res_mrmr;
	tree<double>::iterator top_mrmr;
	
	
	top_mrmr=res_mrmr.begin();
	res_mrmr.insert(top_mrmr, predn);
	tree<double>::iterator it_mrmr_local=res_mrmr.begin(),it_mrmr_local2=it_mrmr_local;
	int target_depth=max_elements, max_depth=0;
	int max_depth_local=2;
	
	while (res_tmp_new.depth(it_local)<target_depth && it_local!=res_tmp_new.end()) {
		max_depth=res_tmp_new.depth(it_local);
		while(it_local!=res_tmp_new.end()) {
			if(cnt!=0){
				it_local2=it_local; it_mrmr_local2=it_mrmr_local;
				while(res_tmp_new.depth(it_local2)<max_depth){
					it_local2++;it_mrmr_local2++;
					
				}
				while(it_local2!=res_tmp_new.begin()){
					nsub[res_tmp_new.depth(it_local2)]=*(it_local2);
					
					it_local2=res_tmp_new.parent(it_local2);
					it_mrmr_local2=res_mrmr.parent(it_mrmr_local2);
				}
			}
			for (unsigned int i=0;i<=max_depth;++i){	
				prev_sel[i]=nsub[i+1];		
			}
			
			////////////
			// initialize vec_mean and vec_sd for bootstrapping to -1000 if variable is not supposed to be tested (target or prev selected) otherwise 0
			////////////
			
			for(unsigned int k=0;k< n;++k){
				vec_mean[k]=0;vec_sd[k]=0;
			}
			for(unsigned int k=0;k<=max(res_tmp_new.depth(it_local),max_depth) ;++k){
				vec_mean[nsub[k]-1]=-1000;	vec_sd[nsub[k]-1]=-1000;
			}

			for(unsigned int k=0;k< n;++k){
				if(vec_mean[k]!= (-1000)){
					vec_mean[k]=mrnet_onegene( mat_info, n,min(cnt,max_elements_tmp),prev_sel, nsub[0], (k+1)); vec_sd[k]=0;
				
				}
				mrmr_vec_sort[k]=vec_mean[k];
			}
			sort(mrmr_vec_sort,mrmr_vec_sort+n);
	
			tmp_val_max=mrmr_vec_sort[n-maxnsol-1];
			int cnt_loop_max=0;
			
			while (res_tmp_new.depth(it_local)<max_depth) {
				it_local++;it_mrmr_local++;
			}
			it_local2=it_local;it_mrmr_local2=it_mrmr_local;
			it_local2++;it_mrmr_local2++;

			
			for(int k=0;k<n;k++){
				if(vec_mean[k]>tmp_val_max){
					vec_sol_local[cnt_loop_max]=k+1;
					vec_sol_local_mrmr[cnt_loop_max]=vec_mean[k];
					cnt_loop_max++;
				}
			}

			for(int k=maxnsol-1;k>=0;k--){
				res_tmp_new.append_child(it_local,vec_sol_local[k]);
				res_mrmr.append_child(it_mrmr_local,vec_sol_local_mrmr[k]);
			}
			
			if(res_tmp_new.depth(it_local)>0){
				it_local=it_local2;it_mrmr_local=it_mrmr_local2;
			}else{
				it_local++;it_mrmr_local++;
			}
			cnt++;
		}
		cnt++; max_elements_tmp++; 	
		ndelete= -1;
		while (ndelete!=0 ) {
			ndelete=verify_equivalentset_nparents (res_tmp_new, res_tmp_new.begin(),res_tmp_new.end(),res_mrmr, maxnsol);
		}
		
		remove_childless_nodes(res_tmp_new, res_mrmr,max_depth_local+1);
		
		it_local=res_tmp_new.begin_leaf();it_mrmr_local=res_mrmr.begin_leaf();
		max_depth_local++;
	}
	res=res_tmp_new;

	bootstrap_tree(res,res_mrmr, data, namat,  nsamples, n, rep_boot, msurv, ustrat, x2, cl2,st, se, weights,strat, N, outx, lenS, lenU);

}


SEXP mrmr_cIndex_ensemble_remove( SEXP Rdata, SEXP Rnamat, SEXP Rmaxparents, SEXP Rnvar, SEXP Rnsample, SEXP Rpredn, SEXP Rnpredn, SEXP Rrep_boot, SEXP Rmaxnsol, SEXP Rthreshold, SEXP Rmsurv, SEXP Rustrat, SEXP Rcl2, SEXP Rst, SEXP Rse, SEXP Rweights, SEXP Rstrat, SEXP RN, SEXP Routx, SEXP RlenS, SEXP RlenU){
	// Rdata:		data should be passed as vector, variable-wise appended
	// Rmaxparents:	number of maximum number of parents
	// Rnvar:		number of variables in the dataset
	// Rnsample:	number of samples in the dataset
	// Rpredn:		vector of target genes to consider
	// Rnpredn:		number of target genes (number of elements in Rpredn)
	// Rrep_boot:	how many bootstrap iterations
	// Rmaxnsol:	maximum number of children for each node at each step
	
	double *data, *threshold;
	const int* maxparents, * nvar, *nsample, *maxnsol;
	
	int *predn, *rep_boot,*res,*res_all,*res_all2, *namat;
	int vec_tmp;
	const int *npredn;
	
	
	double *x2, *st, *weights;
	int *msurv, *ustrat, *cl2, *se, *strat, *N, *outx, *lenS, *lenU;
	
	SEXP Rres;
	
	srand (time(NULL));
	PROTECT(Rdata = AS_NUMERIC(Rdata));
	PROTECT(Rnamat = AS_INTEGER(Rnamat));
	PROTECT(Rmaxparents= AS_INTEGER(Rmaxparents));
	PROTECT(Rnvar= AS_INTEGER(Rnvar));	
	PROTECT(Rnsample= AS_INTEGER(Rnsample));
	PROTECT(Rpredn = AS_INTEGER(Rpredn));
	PROTECT(Rnpredn = AS_INTEGER(Rnpredn));
	PROTECT(Rrep_boot = AS_INTEGER(Rrep_boot));
	PROTECT(Rmaxnsol= AS_INTEGER(Rmaxnsol));
	PROTECT(Rthreshold= AS_NUMERIC(Rthreshold));
	
	data=NUMERIC_POINTER(Rdata);
	namat=INTEGER_POINTER(Rnamat);
	maxparents = INTEGER_POINTER(Rmaxparents);
	nvar= INTEGER_POINTER(Rnvar);
	nsample= INTEGER_POINTER(Rnsample);
	predn= INTEGER_POINTER(Rpredn);
	npredn= INTEGER_POINTER(Rnpredn);
	rep_boot= INTEGER_POINTER(Rrep_boot);	
	maxnsol= INTEGER_POINTER(Rmaxnsol);
	threshold = NUMERIC_POINTER(Rthreshold);
	
	msurv=INTEGER_POINTER(Rmsurv);
	ustrat=INTEGER_POINTER(Rustrat);
	cl2=INTEGER_POINTER(Rcl2);
	st=NUMERIC_POINTER(Rst);
	se =INTEGER_POINTER(Rse);
	weights=NUMERIC_POINTER(Rweights);
	strat=INTEGER_POINTER(Rstrat);
	N =INTEGER_POINTER(RN);
	outx=INTEGER_POINTER(Routx);

	lenS=INTEGER_POINTER(RlenS);
	lenU=INTEGER_POINTER(RlenU);
	
	tree<int> res_tree;
	tree<int>::iterator top,one;
	tree<int>::breadth_first_queued_iterator it_final;
	
	top=res_tree.begin();
	int length_res=0;
	int length_res_old;
	for(unsigned int i=0;i< *npredn;++i){
	//	std::cout<<"model for node "<<predn[i]<< " is being built!"<<std::endl;
		one=res_tree.insert(top, predn[i]);
		
		//build ensemble tree
		mrmr_ensemble_one_gene_remove(res_tree, one, data,namat,*nsample,(*nvar+1),*maxparents,predn[i],*rep_boot, *maxnsol, *threshold, msurv, ustrat, x2, cl2,st, se, weights,strat, N, outx, *lenS, lenU);
		
		////////////////////////
		//convert tree to vector
		////////////////////////
		int *tmp_nchildren,*res_tmp;
		res_tmp=new int [2*(res_tree.size())+1];
		tmp_nchildren= new int [(res_tree.size())];
		
		it_final=res_tree.begin_breadth_first();
		int cnt=1,cnt2=0;
		
		res_tmp[0]=res_tree.size();
		int rootdepth=res_tree.depth(it_final);
		while(it_final!=res_tree.end_breadth_first()) {
			res_tmp[cnt]=*it_final;
			tmp_nchildren[cnt-1]=(res_tree.number_of_children(it_final));
			cnt++;
			++it_final;		
		}
		////////////////
		//save in final result vector
		////////////////
		length_res_old=length_res;
		length_res+=2*(res_tree.size())+1;
		int *res_all, *res_old;
		int ind=0;
		res_all=new int[length_res];
		if(length_res_old>0){
			for(unsigned int k=0;k<length_res_old;k++){
				res_all[k]=res_old[k];
			}
		}
		
		for(unsigned int k=0;k<=res_tree.size();k++){
			res_all[length_res_old+k]=res_tmp[k];
		}
		for(unsigned int k=0;k<res_tree.size();k++){
			res_all[length_res_old+k+res_tree.size()+1]=tmp_nchildren[k];
		}
		
		delete [] res_old;
		res_old=new int[length_res];
		for(unsigned int k=0;k<length_res;k++){
			res_old[k]=res_all[k];
		}
		
		delete [] res_all;
		if(i==(*npredn-1)){
			
			PROTECT(Rres = NEW_INTEGER(length_res));
			res = INTEGER_POINTER(Rres);
			for(unsigned int k=0;k<length_res;k++){
				res[k]=res_old[k];
			}
			delete [] res_old;
		}
		
		////////////////
		//erase old tree
		////////////////
		
		delete [] tmp_nchildren;
		delete [] res_tmp;
		res_tree.erase(res_tree.begin());
	}
	UNPROTECT(11);
	
	return Rres;
}
